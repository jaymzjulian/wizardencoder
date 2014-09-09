#include "fastsid.h"
#include "wrapper.h"
#include <string.h>
#include <stdio.h>

fastsid::State fastsid::read_state() {
	// if this does not work, we're fucked
	return me;
}


fastsid::fastsid() {
	memset(&me.psid, 0, sizeof(me.psid));
}

void fastsid::write_state(State &state) {
	// copy the data how the compiler things ke should...
	// me.psid=state.psid;
	memcpy(&me.psid, &state.psid, sizeof(sound_t));

	// fixup the pointers in the psid structure
	for(int v=0;v<3;v++)
	{
		me.psid.v[v].s=(&me.psid);
		me.psid.v[v].vnext=(&me.psid.v[(v+1)%3]);
		me.psid.v[v].vprev=(&me.psid.v[(v+2)%3]);
		me.psid.v[v].d=me.psid.d+v*7;
	}
}

void fastsid::reset() {
	fastsid_reset(&me.psid, 0);
}

void fastsid::write(int reg, int val) {
	fastsid_store(&me.psid, reg, val);
}

int fastsid::clock(int &cycle, short *buffer, int maxsamples, int interleave) {
	// FIXME: of course this needs a better thing....
	int samplecount=(cycle*samplefreq)/clockspeed;
	if(samplecount>maxsamples)
		samplecount=maxsamples;
	me.psid.maincpu_clk+=cycle;
	cycle-=((samplecount*clockspeed)/samplefreq);
	me.psid.maincpu_clk-=cycle;
	return fastsid_calculate_samples(&me.psid, buffer, samplecount, interleave, 0);
}

void fastsid::set_sampling_parameters(int cpufreq, int type, int sampleFreq) {
	printf("FASTSID: clock=%d, sample=%d\n", cpufreq, sampleFreq);
	clockspeed=cpufreq;
	samplefreq=sampleFreq;
	fastsid_init(&me.psid, sampleFreq, cpufreq);
}
