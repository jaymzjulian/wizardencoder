#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fftw3.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>


#ifdef USE_OPENMP
#include <omp.h>
#else
int omp_get_max_threads() { return 1; }
int omp_get_thread_num() { return 0; }
#endif
#include "sidc.h"
#include "loud.h"
#include "ini.h"

#define DUMBARRAYSIZE 131072

#define DEDUPEOP
//#define NOSR

#define MAXTHREADS 32
#define MAXOPERATORS 32
//#define BENCHMARK (100/frameRange)

#define USE_FASTSID
#undef USE_RESID
#undef USE_RESIDFP

#ifdef USE_RESID
#include "resid/sid.h"
using namespace reSID;
#define SAMPLETYPE SAMPLE_FAST
#endif

#ifdef USE_RESIDFP
#include "resid-fp/sid.h"
#define SID SIDFP
#define SAMPLETYPE SAMPLE_INTERPOLATE
#endif

#ifdef USE_FASTSID
#include "fastsid/wrapper.h"
#define SID fastsid
#define SAMPLETYPE 0
#endif

// sensble defaults
int numChannels=3;
int maxIter=2;
double inputAmp=1.0;
int frameRange=1;
int popSize=1024;
int eliteSize=64;
int numOperators=4*frameRange;
int preWindowFrames=4;
int updateTime=1;
int simpleChannels=0;

signed short *outputBuffer;
signed short *throwawayBuffer;
int outBufferOffset;
int workBufferOffset;
signed short *inputBuffer;
signed short *workBuffer[MAXTHREADS];
double *loudMult;

// Benchtime:
// 13 seconds for '4'
// 20 seconds for '2'
// 30 seconds for '1'
int speedScale=1;




#define THROWAWAY (popSize-eliteSize)

#define GSIZE (numOperators*sizeof(struct so))


#define PALFRAME (19656*updateTime)
#define SINGLEFS (879*updateTime)
#define SHORTFS (870*updateTime)

#define FRAMESAMPLES (SINGLEFS*frameRange)
#define SHORTSAMPLES (SHORTFS*frameRange)
// scale this down by speedScale
#define FULLFFTSAMPLES (1024*frameRange*preWindowFrames)
#define FFTSAMPLES (FULLFFTSAMPLES/speedScale)
#define FFTSCALERATE ((65536*SHORTFS*preWindowFrames*frameRange)/FULLFFTSAMPLES)

int palFrame;
fftw_complex *srcIn, *srcOut, *dstIn[MAXTHREADS], *dstOut[MAXTHREADS];
fftw_plan srcPlan, dstPlan[MAXTHREADS];
FILE *mmmm;

//#define NOSR
//#define SRONLY

enum commands {
#ifndef SRONLY
	set_ad,
#endif
#ifndef NOSR
	set_sr,
#endif
	set_freq,
	set_ctrl,
	set_pw,
	last_command
};

#define NUMCOMMANDS (last_command)

const char *commandToString(int command) {
	switch(command) {
#ifndef SRONLY
		case set_ad: return "set_ad";
#endif
#ifndef NOSR
		case set_sr: return "set_sr";
#endif
		case set_freq: return "set_freq";
		case set_ctrl: return "set_ctrl";
		case set_pw: return "set_pw";
	}
	return "ERROR";
}	


typedef unsigned char byte;
struct so {
	byte frame;
	byte channel;
	byte command;
	byte value;
};

// we're going to pack the shit out of this instead
// FIXME: obviously, this should not be statically defined,
// but instead an object that takes a new or whatever - but 
// for now, this will do...
union sidOutput {
	struct so s[MAXOPERATORS];
	byte gval[MAXOPERATORS*sizeof(struct so)];
};

struct gpop {
	sidOutput population;
	double fitness;
};

byte simpleCtrl[]={
	0x10, 0x20, 0x40, 0x80,
	0x11, 0x21, 0x41, 0x81
};

byte possibleCtrl[]={
	0x10, 0x20, 0x40, 0x80,
	0x11, 0x21, 0x41, 0x81,
	0x12, 0x22, 0x42, 0x82,
	0x13, 0x23, 0x43, 0x83,
	0x14, 0x24, 0x44, 0x84,
	0x15, 0x25, 0x45, 0x85
};

void writeAsm(FILE *outasm, sidOutput o)
{
	for(int cframe=0;cframe<frameRange;cframe++) {
		for(int c=0;c<numOperators;c++)
		{
			if(o.s[c].frame==cframe){
				int baseChannel=o.s[c].channel;
				int command=o.s[c].command;
				int val=o.s[c].value;

				// FIXME: double val if gate is currently 1x on this channel!
				fprintf(outasm, ".byte %d, %d, %d\n", baseChannel, command, val);
			}
		}
		fprintf(outasm, ".byte $ff\n");
	}
}

void pushSid(sidOutput o, SID *testSid, int f)
{
	for(int c=0;c<numOperators;c++)
	{
		if(o.s[c].frame==f) {
		int baseChannel=o.s[c].channel;
		int command=o.s[c].command;
		int val=o.s[c].value;
		int igate;
		int cfreq;
		switch(command) {
#ifdef NOSR
			case set_ad:
				testSid->write(5+baseChannel*7, val&0xf0);
				testSid->write(6+baseChannel*7, (val&0x0f)|0xf0);
				break;
#elif defined(SRONLY)
			case set_sr:
				testSid->write(5+baseChannel*7, 0x00);
				testSid->write(6+baseChannel*7, val);
				break;
#else
			case set_ad:
				testSid->write(5+baseChannel*7, val);
				break;
			case set_sr:
				testSid->write(6+baseChannel*7, val);
				break;
#endif
			case set_pw:
				val=128;
				testSid->write(3+baseChannel*7, val>>4);
				testSid->write(2+baseChannel*7, (val<<4)&255);
				break;
			case set_freq:
#ifndef USE_FASTSID
				igate=testSid->voice[baseChannel].wave.waveform;
#else
				igate=testSid->me.psid.v[baseChannel].d[4]>>4;
#endif
				cfreq=sidFreq[val];
				if(igate==1)
					cfreq*=2;
				if(cfreq>65535)
					cfreq=65535;

				testSid->write(0+baseChannel*7, cfreq&255);
				testSid->write(1+baseChannel*7, cfreq>>8);
				break;
			case set_ctrl:
				// only allow modulation in channel 2!
				int realGate;
				if(baseChannel<simpleChannels)
					realGate=simpleCtrl[val%sizeof(simpleCtrl)];
				else
					realGate=possibleCtrl[val%sizeof(possibleCtrl)];
				testSid->write(4+baseChannel*7, realGate);
				break;
			default:
				printf("INVALID SID COMMAND: %d\n", command);
				abort();
		}
		}
	}
}

bool opSort(struct so a, struct so b)
{
	if(a.frame==b.frame)
	{
		return a.channel < b.channel;
	}
	return a.frame < b.frame;
}

void filterSingle(sidOutput &currentSid, int c) {
		if(currentSid.s[c].command==set_ctrl)
		{
			if(currentSid.s[c].channel<simpleChannels)
				currentSid.s[c].value%=sizeof(simpleCtrl);
			else
				currentSid.s[c].value%=sizeof(possibleCtrl);
		}
		currentSid.s[c].frame%=frameRange;
		currentSid.s[c].command%=NUMCOMMANDS;
		currentSid.s[c].channel%=numChannels;
#ifdef DEDUPEOP
		for(int d=0;d<numOperators;d++)
		{
			if(c!=d) {
				if(currentSid.s[c].frame == currentSid.s[d].frame && currentSid.s[c].command == currentSid.s[d].command && currentSid.s[c].channel==currentSid.s[d].channel) {
					currentSid.s[c].command=rand()%NUMCOMMANDS;
					currentSid.s[c].channel=rand()%numChannels;
					currentSid.s[c].frame=rand()%frameRange;
					filterSingle(currentSid, c);
				}
			}
		}
#endif

}

void filterValues(sidOutput &currentSid) {
	for(int c=0;c<numOperators;c++) {
		filterSingle(currentSid, c);
	}	
}

void singleRandom(sidOutput &currentSid, int c) {
	currentSid.s[c].frame=rand()%frameRange;
	currentSid.s[c].command=rand()%NUMCOMMANDS;
	currentSid.s[c].channel=rand()%numChannels;
	currentSid.s[c].value=rand()%256;
	filterSingle(currentSid, c);
}

void fullyRandomPop(sidOutput &currentSid) {
	for(int c=0;c<numOperators;c++) {
		currentSid.s[c].frame=rand()%frameRange;
		currentSid.s[c].command=rand()%NUMCOMMANDS;
		currentSid.s[c].channel=rand()%numChannels;
		currentSid.s[c].value=rand()%256;
	}
	filterValues(currentSid);
	// help it along a little!
	std::sort(currentSid.s, currentSid.s+numOperators, opSort);

}

// Do the src FFT once
void setupSrcGlobal() {
	// input is 44100, output is 44100/speedScale - so, really, we need to multiply by 
	// speedScale here.... if we got all of our #define's right, it all works out!
	double srcMax=0-32768.0;
	double srcMin= 32768.0;
	for(int c=0;c<FFTSAMPLES;c++) {
		srcIn[c][0]=((double)inputBuffer[(c*(FFTSCALERATE*speedScale))>>16])*inputAmp;
		srcIn[c][1]=0;
		if(srcMax<srcIn[c][0])
			srcMax=srcIn[c][0];
		if(srcMin>srcIn[c][0])
			srcMin=srcIn[c][0];
	}
	double srcCentre=(srcMax+srcMin)/2.0;
	for(int c=0;c<FFTSAMPLES;c++) 
		srcIn[c][0]-=srcCentre;
	printf("Samples per frame: %d/%d\n", FRAMESAMPLES, preWindowFrames);
	printf("Samples per fft: %d/%d\n", FFTSAMPLES, preWindowFrames);
	printf("Used %d samples\n", (FFTSAMPLES*(FFTSCALERATE*speedScale))>>16);
	printf("Cmp %d samples\n", (FFTSAMPLES*FFTSCALERATE)>>16);
	printf("Centre: %f\n", srcCentre);
	fftw_execute(srcPlan);
	
	// we need to copy the pre-window from the rendered output buffer
	// workBuffer is 44100, but of course we want this divided by speedscale if
	// speedscale!=1 - but at 1, memcpy is way faster!
	workBufferOffset=outBufferOffset/speedScale;
	for(int t=0;t<omp_get_max_threads();t++)
	{
		if(speedScale==1)
			memcpy(workBuffer[t], outputBuffer, outBufferOffset*sizeof(short));
		else
			for(int c=0;c<workBufferOffset;c++)
				workBuffer[t][c]=outputBuffer[c*speedScale];
	}
	fwrite(outputBuffer, sizeof(short), outBufferOffset, mmmm);
	fflush(mmmm);
}

double testFitness(SID *testSid, sidOutput currentSid, int baseCycle) {
	int t=omp_get_thread_num();
	int genSmpl=workBufferOffset;
	int cycle=baseCycle;
	double loss=0;

	for(int cframe=0;cframe<frameRange;cframe++) {
		// IMPORTANT: push sid, THEN call resid!
		pushSid(currentSid, testSid, cframe);
		cycle+=palFrame;
		genSmpl+=testSid->clock(cycle, workBuffer[t]+genSmpl, DUMBARRAYSIZE, 1);
	}
	// some extra, just in case - more window is good window :)
	cycle+=(palFrame/4);
	genSmpl+=testSid->clock(cycle, workBuffer[t]+genSmpl, DUMBARRAYSIZE, 1);

	// compare samples because fuck everything
	double dstMax=0-32768.0;
	double dstMin= 32768.0;
	for(int c=0;c<FFTSAMPLES;c++) {
		// if this cast isn't here, the algorythm breaks.  why??!
		// seriously, this should work...
		dstIn[t][c][0]=(double)workBuffer[t][((c*FFTSCALERATE)>>16)];
		dstIn[t][c][1]=0;
		if(dstMax<dstIn[t][c][0])
			dstMax=dstIn[t][c][0];
		if(dstMin>dstIn[t][c][0])
			dstMin=dstIn[t][c][0];
			
	}
	double dstCentre=(dstMax+dstMin)/2.0;
	for(int c=0;c<FFTSAMPLES;c++) {
		dstIn[t][c][0]-=dstCentre;
	}


	fftw_execute(dstPlan[t]);

	// skip the DC offset - it'll be balls out wrong anyhow...
	for(int c=1;c<((FFTSAMPLES/2)-1);c++) {
#if 0
		double ia=sqrt(srcOut[c][0]*srcOut[c][0]+srcOut[c][1]*srcOut[c][1]);
		double ib=sqrt(dstOut[t][c][0]*dstOut[t][c][0]+dstOut[t][c][1]*dstOut[t][c][1]);
		double smplLoss=(ia-ib);
		smplLoss*=smplLoss;
		smplLoss*=loudMult[c];
#else
		// This is wrong, but for some reason i think the
		// other is too, so it can stay here as a comment
		double ia=srcOut[c][0]-dstOut[t][c][0];
		double ib=srcOut[c][1]-dstOut[t][c][1];
		ia*=loudMult[c];
		ib*=loudMult[c];
		// Normalise output, to help the humans a little
		ia/=FFTSAMPLES;
		ib/=FFTSAMPLES;

		// break this out, because reasons...
		double smplLoss=((ia*ia)+(ib*ib));
#endif
		loss+=smplLoss;
	}
	return loss;
}

inline bool sfunc(struct gpop a, struct gpop b) {
	return a.fitness < b.fitness;
}

static int handler(void *user, const char *section, const char *name, const char *value) {
	#define MATCH(s, n) (strcmp(section, s)==0 && strcmp(name, n)==0)
	if(MATCH("encoder", "maxIter"))
		maxIter=atoi(value);
	else if(MATCH("encoder", "inputAmp"))
		inputAmp=atof(value);
	else if(MATCH("encoder", "frameRange"))
		frameRange=atoi(value);
	else if(MATCH("encoder", "windowSize"))
		preWindowFrames=atoi(value);
	else if(MATCH("encoder", "numChannels"))
		numChannels=atoi(value);
	else if(MATCH("encoder", "numOperators"))
		numOperators=atoi(value);
	else if(MATCH("encoder", "updateTime"))
		updateTime=atoi(value);
	else if(MATCH("search", "popSize"))
		popSize=atoi(value);
	else if(MATCH("search", "eliteSize"))
		eliteSize=atoi(value);
	else if(MATCH("speedhacks", "speedScale"))
		speedScale=atoi(value);
	else
		return 0;
	return 1;
}

void calculateLoudness() {
	loudMult=new double[FFTSAMPLES];
	for(int c=0;c<FFTSAMPLES;c++)
	{
#if 0
		double freq=44100.0/FFTSAMPLES;
		freq*=c;
		// speedscale drops freq buckets
		freq/=speedScale;

		int bm=0;
		double d=999;
		// sizeof is the wrong tool here - the doubles should be
		// a vector, not an array :)
		// also, this is terrible :)
		bm=(sizeof(freqTable)/sizeof(double))-1;
		for(unsigned int e=0;e<(sizeof(freqTable)/sizeof(double));e++)
		{
			double fdif=fabs(freq-freqTable[e]);
			if(fdif<d) {
				d=fdif;
				bm=e;
			}
		}
		//printf("%d, %d, %d\n", sizeof(freqTable), sizeof(loudnessDB), bm);////

		//printf("%f - %f\n", freq, freqTable[bm]);
		
		// for now, set our multiplier such that a 40db tone would be appropriate
		// later we'll import something real!
		//double dbDif=loudnessDB[bm]-40.0;
		//double multiplier=pow(2, dbDif/10.0);
		//double multiplier=loudnessDB[bm]/40.0;
		double multiplier=1.0;
		//multiplier=sqrt(multiplier);
		//printf("%f: %f, %f\n", freq, dbDif, multiplier);
		loudMult[c]=1.0/multiplier;
#else
		loudMult[c]=1.0;
#endif
	}
	//exit(1);
}

int main(int argc, char **argv) {
	int cycle=0;
	if(argc<4)
	{
		printf("Usage: %s [infile.raw] [outfile.asm] [encoderconfig.ini]\n", argv[0]);
		exit(1);
	}
	
	printf("Config: %s\n", argv[3]);
	if(ini_parse(argv[3], handler, NULL) < 0) {
		printf("Could not parse %s\n", argv[3]);
		exit(1);
	}
	
	printf("Super algorythm codec v2 - %d threads\n", omp_get_max_threads());
	printf("FFT: %d. Frame: %d (LastSmpl: %d)\n", FFTSAMPLES, FRAMESAMPLES, (FFTSAMPLES*FFTSCALERATE)>>16);
	printf("FFTSCALERATE: %d\n", FFTSCALERATE);
	printf("Possible commands: %d\n", NUMCOMMANDS);

	srcIn=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES * 2);
	srcOut=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES * 2);
	srcPlan=fftw_plan_dft_1d(FFTSAMPLES, srcIn, srcOut, FFTW_FORWARD, FFTW_MEASURE);
	for(int t=0;t<omp_get_max_threads();t++)
	{
		dstIn[t]=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES * 2);
		dstOut[t]=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES * 2);
		dstPlan[t]=fftw_plan_dft_1d(FFTSAMPLES, dstIn[t], dstOut[t], FFTW_FORWARD, FFTW_MEASURE);
		workBuffer[t]=(short *)calloc(DUMBARRAYSIZE,sizeof(short));
	}

	calculateLoudness();

	palFrame=PALFRAME;

	outputBuffer=(short *)calloc(DUMBARRAYSIZE,sizeof(short));
	throwawayBuffer=(short *)calloc(DUMBARRAYSIZE,sizeof(short));
	inputBuffer=(short *)calloc(DUMBARRAYSIZE,sizeof(short));

	// time for a goddamned sid!
	// we have two sid types we use - a test sid, configured to
	// be fast, and a normal sid that is configued to be good
	SID *testSid[MAXTHREADS];
	for(int t=0;t<omp_get_max_threads();t++)
	{
	testSid[t]=new SID();
#if 0
	testSid[t]->set_chip_model(MOS8580);
	testSid[t]->set_voice_mask(0x0f);
#endif
#ifndef USE_FASTSID
	testSid[t]->input(0);
	testSid[t]->enable_filter(false);
#endif
	//testSid[t]->enable_external_filter(false);
	testSid[t]->set_sampling_parameters(985248, SAMPLETYPE, 44100/speedScale);
	testSid[t]->reset();
	for(int c=0;c<32;c++)
		testSid[t]->write(c, 0);
	testSid[t]->write(0x18, 15);
	cycle=palFrame*50;
	testSid[t]->clock(cycle, workBuffer[t], DUMBARRAYSIZE, 1);
	}

	SID *outSid=new SID();
	SID *lowresSid=new SID();
#if 0
	outSid->set_chip_model(MOS8580);
	outSid->set_voice_mask(0x0f);
	lowresSid->set_chip_model(MOS8580);
	lowresSid->set_voice_mask(0x0f);
#endif

#ifndef USE_FASTSID
	outSid->input(0);
	outSid->enable_filter(false);
	lowresSid->input(0);
	lowresSid->enable_filter(false);
#endif
	//outSid->enable_external_filter(false);
	outSid->set_sampling_parameters(985248, SAMPLETYPE, 44100);
	outSid->reset();
	lowresSid->set_sampling_parameters(985248, SAMPLETYPE, 44100/speedScale);
	lowresSid->reset();
	for(int c=0;c<32;c++)
	{
		lowresSid->write(c, 0);
		outSid->write(c, 0);
	}
	outSid->write(0x18, 15);
	lowresSid->write(0x18, 15);

	// do this twice to nuke instabilities in resid startup
	cycle=palFrame*100;
	outSid->clock(cycle, outputBuffer, DUMBARRAYSIZE, 1);
	cycle=palFrame*100;
	lowresSid->clock(cycle, outputBuffer, DUMBARRAYSIZE, 1);
	cycle=palFrame*100;
	outSid->clock(cycle, outputBuffer, DUMBARRAYSIZE, 1);
	cycle=palFrame*100;
	lowresSid->clock(cycle, outputBuffer, DUMBARRAYSIZE, 1);


	//FILE *inputFile=fopen("harry.raw", "rb");
	FILE *inputFile=fopen(argv[1], "rb");
	FILE *fml=fopen("output.raw", "wb");
	FILE *ref=fopen("referece.raw", "wb");
	mmmm=fopen("mmmm.raw", "wb");
	FILE *outasm=fopen(argv[2], "wt");
	//fseek(inputFile, 88200*10, 0);

	int readLen;
	cycle=0;
	int frameCount=0;
	// these are all 44100 buffers
	int bufferOffset=(preWindowFrames-1)*FRAMESAMPLES;
	outBufferOffset=(preWindowFrames-1)*FRAMESAMPLES;
	int readWindow=FRAMESAMPLES*(preWindowFrames+1);
#ifdef BENCHMARK
	#error broken benchmarking!
	while(readLen=fread(inputBuffer, sizeof(short), FRAMESAMPLES, inputFile) && frameCount<BENCHMARK) {
#else
	while((readLen=fread(inputBuffer+bufferOffset, sizeof(short), readWindow-bufferOffset, inputFile))) {
#endif
		bufferOffset+=readLen;
		printf("Read: %d\n", readLen);
		printf("BufferTail: %d / %d\n", bufferOffset, outBufferOffset);
		struct gpop mypop[popSize];
		int baseCycle=cycle;
		// we get our base state from a sid in the same
		// resolution as testsid, since some impls (like fastsid)
		// have state dependent on freq
		SID::State baseState=lowresSid->read_state();

		// do the FFT and processing on the source
		setupSrcGlobal();

		// do initial pop
		for(int c=0;c<popSize;c++)
		{
			// this causes memory errors, so we're just
			// not doing it this time - it makes vibrato harder
			// to track though!
			//if(frameCount==0 || c>THROWAWAY) {
			fullyRandomPop(mypop[c].population);
			//}
		}
		
		/* parallel test fitness */
#pragma omp parallel for shared(baseState,mypop,testSid,baseCycle)
		for(int c=0;c<popSize;c++)
		{
			int t=omp_get_thread_num();
			testSid[t]->write_state(baseState);
			mypop[c].fitness=testFitness(testSid[t], mypop[c].population, baseCycle);
		}
		


		printf("go - %d\n", maxIter);
		std::sort(mypop, mypop+popSize, sfunc);
		double bf=mypop[0].fitness;
		for(int i=0;i<maxIter;i++)
		{
			//printf("Best: %f\n", mypop[0].fitness);
			// sort the population
			std::sort(mypop, mypop+popSize, sfunc);
			if(mypop[0].fitness < bf)
			{
				bf=mypop[0].fitness;
				printf("BF: %f (%d) (%ld)\n", bf, i, GSIZE);
				//printf("%x\n", possibleCtrl[mypop[0].population.s.ctrl2%sizeof(possibleCtrl)]);
				i=0;
			}

			// replace the bottom half - seperate this so it's non-parallel
			for(int c=eliteSize;c<popSize;c++)
			{
				unsigned int splitPoint=rand()%(sizeof(sidOutput)*8);
				//int p1=rand()%popSize;
				//int p2=rand()%popSize;
				int p1=rand()%eliteSize;
				int p2=rand()%eliteSize;
				for(unsigned int bit=0;bit<GSIZE*8;bit++)
				{
					int cbyte=bit>>3;
					int cbit=bit&7;
					int f;
					if(bit<splitPoint)
						f=p1;
					else
						f=p2;
					//mypop[c].population.gval[cbyte]=mypop[f].population.gval[cbyte];
					if(cbit==0)
						mypop[c].population.gval[cbyte]=0;
					mypop[c].population.gval[cbyte]|=(mypop[f].population.gval[cbyte])&(1<<(7-cbit));
				}
				filterValues(mypop[c].population);

				// mutate one command of each of these
				int randBit=rand()%numOperators;
				singleRandom(mypop[c].population, randBit);
			}

#pragma omp parallel for shared(baseState,mypop,testSid)
			for(int c=eliteSize;c<popSize;c++)
			{
				int t=omp_get_thread_num();
				testSid[t]->write_state(baseState);
				mypop[c].fitness=testFitness(testSid[t], mypop[c].population, baseCycle);
			}
			/* parallel ends here! */
		}
		printf("PCOMPLETE\n");
		// final sort the population
		std::sort(mypop, mypop+popSize, sfunc);
		//std::sort(mypop[0].population.s, mypop[0].population.s+numOperators, opSort);
		printf("Fitness:%f\n", mypop[0].fitness);
		for(int c=0;c<numOperators;c++)
		{
			int cval=mypop[0].population.s[c].value;
			if(mypop[0].population.s[c].command == set_ctrl)
			{
				if(mypop[0].population.s[c].channel<simpleChannels)
					cval=simpleCtrl[cval];
				else
					cval=possibleCtrl[cval];
			}
			printf("%d/%d/%s/%x\n", mypop[0].population.s[c].frame, mypop[0].population.s[c].channel, commandToString(mypop[0].population.s[c].command), cval);
		}

		frameCount++;
		int genSmpl=0;
		cycle=baseCycle;
		int lcycle=baseCycle;
		printf("basecyc:%d\n", baseCycle);
		for(int cframe=0;cframe<frameRange;cframe++) {
			pushSid(mypop[0].population, outSid, cframe);
			pushSid(mypop[0].population, lowresSid, cframe);
			// keep cycle/lcycle seperate, even though we're throwing it out, 
			// since it's updated by reference
			lcycle+=palFrame;
			lowresSid->clock(lcycle, throwawayBuffer+(genSmpl+outBufferOffset), DUMBARRAYSIZE, 1);
			// now the real one!
			cycle+=palFrame;
			genSmpl+=outSid->clock(cycle, outputBuffer+(genSmpl+outBufferOffset), DUMBARRAYSIZE, 1);
		}
		outBufferOffset+=genSmpl;
#ifndef USE_FASTSID
		printf("VOL REG: %x\n", outSid->filter.vol);
#endif

		// write our buffers from the current position of read, rather
		// than the head
		printf("Generated %d samples (%d)\n", genSmpl, frameCount*2*frameRange);
		fwrite(outputBuffer+(outBufferOffset-genSmpl), sizeof(short), genSmpl, fml);
		//fwrite(outputBuffer, sizeof(short), genSmpl, fml);
		fflush(fml);
		fwrite(inputBuffer+(preWindowFrames-1)*FRAMESAMPLES, sizeof(short), genSmpl, ref);
		fflush(ref);
		writeAsm(outasm, mypop[0].population);
		fflush(outasm);
		
		// move the input window back
		memmove(inputBuffer, inputBuffer+genSmpl, (readWindow-genSmpl)*sizeof(short));
		bufferOffset-=genSmpl;
		
		// move the output window back
		memmove(outputBuffer, outputBuffer+genSmpl, (outBufferOffset-genSmpl)*sizeof(short));
		outBufferOffset-=genSmpl;
		//if(frameCount==4)
		//	exit(1);
	}


	return 0;
}
