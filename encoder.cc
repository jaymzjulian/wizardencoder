#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fftw3.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>

#include <omp.h>
#include "sid.h"

#define DEDUPEOP

#define MAXTHREADS 32

#define NUMCHANNELS 3

#define MEASUREINDIVIDUAL

using namespace reSID;
short *outputBuffer;
short *inputBuffer;
short *workBuffer[MAXTHREADS];

// benctime:
// 13 seconds for '2'
// 25 seconds for '4'
int maxIter=8;

#define INPUTAMP 1.0

// 13 seconds for '1'
// 24 seconds for '2' (doubles OPERATORS, POPSIZE!)
#define FRAMERANGE 1

// Benchtime:
// 13 seconds for '4'
// 20 seconds for '2'
// 30 seconds for '1'
#define SPEEDSCALE 4

#undef FUNKYCTRL


//#define BENCHMARK (100/FRAMERANGE)

#define FREQSCALE 3

#define POPSIZE (1024*FRAMERANGE)
#define KEEPSIZE (64*FRAMERANGE)
#define THROWAWAY (POPSIZE-KEEPSIZE)


#define PALFRAME 19656

#ifndef MEASUREINDIVIDUAL
	#define FRAMESAMPLES (879*FRAMERANGE)
	#define FFTSAMPLES (((1024*FRAMERANGE)/SPEEDSCALE))
#else
	#define FRAMESAMPLES (879)
	#define FFTSAMPLES (((1024)/SPEEDSCALE))
#endif
#define FFTSCALERATE ((256*FRAMESAMPLES)/FFTSAMPLES)

// 4 operators per frame by default
#define OPERATORS 4*NUMCHANNELS*FRAMERANGE



// WRONG WRONG WRONG
int palFrame;
fftw_complex *srcIn, *srcOut, *dstIn[MAXTHREADS], *dstOut[MAXTHREADS];
fftw_plan srcPlan, dstPlan[MAXTHREADS];

#define NUMCOMMANDS 5

enum commands {
	set_ad,
	set_sr,
	set_freq,
	set_ctrl,
	set_pw
};

const char *commandToString(int command) {
	switch(command) {
		case set_ad: return "set_ad";
		case set_sr: return "set_sr";
		case set_freq: return "set_freq";
		case set_ctrl: return "set_ctrl";
		case set_pw: return "set_pw";
	}

}	


typedef unsigned char byte;
	struct so {
		byte frame ;
		byte channel ;
		byte command ;
		byte value ;
	};

// we're going to pack the shit out of this instead
union sidOutput {
	struct so s[OPERATORS];
	byte gval[OPERATORS*sizeof(struct so)];
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
	for(int cframe=0;cframe<FRAMERANGE;cframe++) {
		for(int c=0;c<OPERATORS;c++)
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

void pushSid(sidOutput o, SID &testSid, int f)
{
	for(int c=0;c<OPERATORS;c++)
	{
		if(o.s[c].frame==f) {
		int baseChannel=o.s[c].channel;
		int command=o.s[c].command;
		int val=o.s[c].value;
		int igate;
		switch(command) {
			case set_ad:
				testSid.write(5+baseChannel*7, val);
				break;
			case set_sr:
				testSid.write(6+baseChannel*7, val);
				break;
			case set_pw:
				testSid.write(3+baseChannel*7, val);
				break;
			case set_freq:
				igate=testSid.voice[baseChannel].wave.waveform;
				if(igate==1){
					testSid.write(0+baseChannel*7, ((val)<<(8-(FREQSCALE-1)))&0xff);
					testSid.write(1+baseChannel*7, ((val)>>((FREQSCALE-1)))&0xff);
				} else {
					testSid.write(0+baseChannel*7, ((val)<<(8-(FREQSCALE)))&0xff);
					testSid.write(1+baseChannel*7, ((val)>>((FREQSCALE)))&0xff);
				}
				break;
			case set_ctrl:
				// only allow modulation in channel 2!
				int realGate;
				if(baseChannel==0)
					realGate=simpleCtrl[val%sizeof(simpleCtrl)];
				else
					realGate=possibleCtrl[val%sizeof(possibleCtrl)];
				testSid.write(4+baseChannel*7, realGate);
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
			if(currentSid.s[c].channel==0)
				currentSid.s[c].value%=sizeof(simpleCtrl);
			else
				currentSid.s[c].value%=sizeof(possibleCtrl);
		}
		if(currentSid.s[c].command==set_pw)
			currentSid.s[c].value%=8;
		currentSid.s[c].frame%=FRAMERANGE;
		currentSid.s[c].command%=NUMCOMMANDS;
#ifdef DEDUPEOP
		for(int d=0;d<OPERATORS;d++)
		{
			if(c!=d) {
				if(currentSid.s[c].frame == currentSid.s[d].frame && currentSid.s[c].command == currentSid.s[d].command && currentSid.s[c].channel==currentSid.s[d].channel) {
					currentSid.s[c].command=rand()%NUMCOMMANDS;
					currentSid.s[c].channel=rand()%NUMCHANNELS;
					currentSid.s[c].frame=rand()%FRAMERANGE;
					filterSingle(currentSid, c);
				}
			}
		}
#endif

}

void filterValues(sidOutput &currentSid) {
	for(int c=0;c<OPERATORS;c++) {
		filterSingle(currentSid, c);
	}	
}

void singleRandom(sidOutput &currentSid, int c) {
	currentSid.s[c].frame=rand()%FRAMERANGE;
	currentSid.s[c].command=rand()%NUMCOMMANDS;
	currentSid.s[c].channel=rand()%NUMCHANNELS;
	currentSid.s[c].value=rand()%256;
	filterSingle(currentSid, c);
}

void fullyRandomPop(sidOutput &currentSid) {
	for(int c=0;c<OPERATORS;c++) {
		currentSid.s[c].frame=rand()%FRAMERANGE;
		currentSid.s[c].command=rand()%NUMCOMMANDS;
		currentSid.s[c].channel=rand()%NUMCHANNELS;
		currentSid.s[c].value=rand()%256;
	}
	filterValues(currentSid);
	// help it along a little!
	std::sort(currentSid.s, currentSid.s+OPERATORS, opSort);

}

// Do the src FFT once
void setupSrcGlobal() {
#ifdef MEASUREINDIVIDUAL
	for(int cframe=0;cframe<FRAMERANGE;cframe++)
#endif
	{
		for(int c=0;c<FFTSAMPLES;c++) {
#ifdef MEASUREINDIVIDUAL
			srcIn[c][0]=((double)inputBuffer[(cframe*FRAMESAMPLES)+((c*FFTSCALERATE)>>8)])*INPUTAMP;
#else
			srcIn[c][0]=((double)inputBuffer[(c*FFTSCALERATE)>>8])*INPUTAMP;
#endif
			srcIn[c][1]=0;
			// noramalise
			srcIn[c][0]/=32768.0;
		}
	}
	fftw_execute(srcPlan);
}

double testFitness(SID &testSid, sidOutput currentSid, int baseCycle) {
	// IMPORTANT: push sid, THEN call resid!
	int t=omp_get_thread_num();
	int genSmpl=0;
	int cycle=baseCycle;
	double loss=0;
	for(int cframe=0;cframe<FRAMERANGE;cframe++) {
		pushSid(currentSid, testSid, cframe);
		cycle+=palFrame;
		genSmpl+=testSid.clock(cycle, workBuffer[t]+genSmpl, 131072, 1);
	}

	// compare samples because fuck everything
#ifdef MEASUREINDIVIDUAL
	for(int cframe=0;cframe<FRAMERANGE;cframe++)
#endif
	{
		for(int c=0;c<FFTSAMPLES;c++) {
#ifdef MEASUREINDIVIDUAL
			dstIn[t][c][0]=workBuffer[t][(cframe*FRAMESAMPLES)+(((c*FFTSCALERATE)/SPEEDSCALE)>>8)];
#else
			dstIn[t][c][0]=workBuffer[t][((c*FFTSCALERATE)/SPEEDSCALE)>>8];
#endif

			dstIn[t][c][1]=0;

			// normalise
			dstIn[t][c][0]/=32768.0;
		}

		fftw_execute(dstPlan[t]);

		for(int c=1;c<FFTSAMPLES/2;c++) {
			//double ia=sqrt(srcOut[t][c][0]*srcOut[t][c][0]+srcOut[t][c][1]*srcOut[t][c][1]);
			//double ib=sqrt(dstOut[t][c][0]*dstOut[t][c][0]+dstOut[t][c][1]*dstOut[t][c][1]);
			//double myLoss=(ia-ib);
			//loss+=myLoss*myLoss;

			// This is wrong, but for some reason i think the
			// other is too, so it can stay here as a comment
			double ia=srcOut[c][0]-dstOut[t][c][0];
			double ib=srcOut[c][1]-dstOut[t][c][1];
			loss+=(ia*ia)+(ib*ib);
		}
	}
	return loss;
}

bool sfunc(struct gpop a, struct gpop b) {
	return a.fitness < b.fitness;
}

int main(int argc, char **argv) {
	printf("Super algorythm codec v2 - %d threads\n", omp_get_max_threads());
	printf("FFT: %d. Frame: %d\n", FFTSAMPLES, FRAMESAMPLES);
	printf("FFTSCALERATE: %d\n", FFTSCALERATE);
	if(argc<3)
	{
		printf("Usage: %s [infile.raw] [outfile.asm]\n", argv[0], argv[1]);
		exit(1);
	}
	srcIn=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES);
	srcOut=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES);
	srcPlan=fftw_plan_dft_1d(FFTSAMPLES, srcIn, srcOut, FFTW_FORWARD, FFTW_MEASURE);
	for(int t=0;t<omp_get_max_threads();t++)
	{
	dstIn[t]=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES);
	dstOut[t]=(fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FFTSAMPLES);
	dstPlan[t]=fftw_plan_dft_1d(FFTSAMPLES, dstIn[t], dstOut[t], FFTW_FORWARD, FFTW_MEASURE);
	workBuffer[t]=(short *)malloc(131072*sizeof(short));
	}

	palFrame=PALFRAME;

	outputBuffer=(short *)malloc(131072*sizeof(short));
	inputBuffer=(short *)malloc(131072*sizeof(short));

	// time for a goddamned sid!
	// we have two sid types we use - a test sid, configured to
	// be fast, and a normal sid that is configued to be good
	SID testSid[MAXTHREADS];
	for(int t=0;t<omp_get_max_threads();t++)
	{
	testSid[t].enable_filter(false);
	testSid[t].reset();
	testSid[t].set_sampling_parameters(985248, SAMPLE_FAST, 44100/SPEEDSCALE);
	testSid[t].write(3, 0x08);
	testSid[t].write(3+7, 0x08);
	testSid[t].write(0x18, 15);
	}

	SID outSid;
	outSid.enable_filter(false);
	outSid.reset();
	outSid.set_sampling_parameters(985248, SAMPLE_FAST, 44100);
	outSid.write(3, 0x08);
	outSid.write(3+7, 0x08);
	outSid.write(0x18, 15);




	//FILE *inputFile=fopen("harry.raw", "rb");
	FILE *inputFile=fopen(argv[1], "rb");
	FILE *fml=fopen("output.raw", "wb");
	FILE *ref=fopen("referece.raw", "wb");
	FILE *outasm=fopen(argv[2], "wt");
	fseek(inputFile, 88200*10, 0);

	int readLen;
	int cycle=0;
	int frameCount=0;
#ifdef BENCHMARK
	while(readLen=fread(inputBuffer, sizeof(short), FRAMESAMPLES, inputFile) && frameCount<BENCHMARK) {
#else
	while(readLen=fread(inputBuffer, sizeof(short), FRAMESAMPLES, inputFile)) {
#endif
		printf("Read: %d\n", readLen);
		double bestLoss=9999999999;
		struct gpop mypop[POPSIZE];
		sidOutput bestSid;
		sidOutput currentSid;
		SID::State baseState=outSid.read_state();
		int baseCycle=cycle;
		int allCtrls[4]={0x0,0x2,0x4,0x6};

		setupSrcGlobal();

		// do initial pop
		for(int c=0;c<POPSIZE;c++)
		{
			if(frameCount==0 || c>THROWAWAY) {
			fullyRandomPop(mypop[c].population);
			}
		}
		
		/* parallel test fitness */
#pragma omp parallel for shared(baseState,mypop,testSid)
		for(int c=0;c<POPSIZE;c++)
		{
			int t=omp_get_thread_num();
#ifdef FUNKYCTRL
			double bp=999999;
			int cnum=0;
			for(int ctype=0;ctype<3;ctype++) {
				mypop[c].population.s.ctrl2&0xf9;
				mypop[c].population.s.ctrl2|allCtrls[ctype];

				testSid[0].write_state(baseState);
				mypop[c].fitness=testFitness(testSid[t], mypop[c].population, baseCycle);
				if(mypop[c].fitness < bp) {
					cnum=ctype;
					bp=mypop[c].fitness;
				}
			}
			mypop[c].fitness=bp;
			mypop[c].population.s.ctrl2&0xf9;
			mypop[c].population.s.ctrl2|allCtrls[cnum];
#else
			testSid[t].write_state(baseState);
			mypop[c].fitness=testFitness(testSid[t], mypop[c].population, baseCycle);
#endif
		}
		


		memset(&bestSid, 0, sizeof(bestSid));
		printf("go - %d\n", maxIter);
		std::sort(mypop, mypop+POPSIZE, sfunc);
		double bf=mypop[0].fitness;
		for(int i=0;i<maxIter;i++)
		{
			//printf("Best: %f\n", mypop[0].fitness);
			// sort the population
			std::sort(mypop, mypop+POPSIZE, sfunc);
			if(mypop[0].fitness < bf)
			{
				bf=mypop[0].fitness;
				printf("BF: %f (%d) (%d)\n", bf, i, sizeof(mypop[0].population.gval));
				//printf("%x\n", possibleCtrl[mypop[0].population.s.ctrl2%sizeof(possibleCtrl)]);
				i=0;
			}

			// replace the bottom half - seperate this so it's non-parallel
			for(int c=KEEPSIZE;c<POPSIZE;c++)
			{
				int splitPoint=rand()%(sizeof(sidOutput)*8);
				//int p1=rand()%POPSIZE;
				//int p2=rand()%POPSIZE;
				int p1=rand()%KEEPSIZE;
				int p2=rand()%KEEPSIZE;
				for(int bit=0;bit<sizeof(mypop[c].population.gval)*8;bit++)
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
				int randBit=rand()%OPERATORS;
				singleRandom(mypop[c].population, randBit);
			}

#pragma omp parallel for shared(baseState,mypop,testSid)
			for(int c=KEEPSIZE;c<POPSIZE;c++)
			{
				int t=omp_get_thread_num();
#ifdef FUNKYCTRL
				double bp=999999;
				int cnum=0;
				for(int ctype=0;ctype<3;ctype++) {
					mypop[c].population.s.ctrl2&0xf9;
					mypop[c].population.s.ctrl2|allCtrls[ctype];

					testSid[t].write_state(baseState);
					mypop[c].fitness=testFitness(testSid[t], mypop[c].population, baseCycle);
					if(mypop[c].fitness < bp) {
						cnum=ctype;
						bp=mypop[c].fitness;
					}
				}
				mypop[c].fitness=bp;
				mypop[c].population.s.ctrl2&0xf9;
				mypop[c].population.s.ctrl2|allCtrls[cnum];
#else
				testSid[t].write_state(baseState);
				mypop[c].fitness=testFitness(testSid[t], mypop[c].population, baseCycle);
#endif
			}
			/* parallel ends here! */
		}
		printf("PCOMPLETE\n");
		// final sort the population
		std::sort(mypop, mypop+POPSIZE, sfunc);
		std::sort(mypop[0].population.s, mypop[0].population.s+OPERATORS, opSort);
		bestSid=mypop[0].population;
		printf("Fitness:%f\n", mypop[0].fitness);
		for(int c=0;c<OPERATORS;c++)
		{
			int cval=mypop[0].population.s[c].value;
			if(mypop[0].population.s[c].command == set_ctrl)
			{
				if(mypop[0].population.s[c].channel==0)
					cval=simpleCtrl[cval];
				else
					cval=possibleCtrl[cval];
			}
			printf("%d/%d/%s/%x\n", mypop[0].population.s[c].frame, mypop[0].population.s[c].channel, commandToString(mypop[0].population.s[c].command), cval);
		}

		frameCount++;
		outSid.write_state(baseState);
		int genSmpl=0;
		cycle=baseCycle;
		for(int cframe=0;cframe<FRAMERANGE;cframe++) {
			pushSid(bestSid, outSid, cframe);
			cycle+=palFrame;
			genSmpl+=outSid.clock(cycle, outputBuffer+genSmpl, 131072, 1);
		}

		printf("Generated %d samples (%d)\n", genSmpl, frameCount*2*FRAMERANGE);
		fwrite(outputBuffer, sizeof(short), genSmpl, fml);
		fflush(fml);
		fwrite(inputBuffer, sizeof(short), genSmpl, ref);
		fflush(ref);
		writeAsm(outasm, bestSid);
		fflush(outasm);
	}


	return 0;
}
