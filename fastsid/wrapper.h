#include "fastsid.h"

class fastsid {
public:
	class State {
	public:
		sound_t psid;
	};
	State me;
	State read_state();
	fastsid();
	void write_state(State &state);
	void reset();
	void write(int reg, int val);
	int clock(int &cycle, short *buffer, int maxsamples, int interleave);
	void set_sampling_parameters(int cpufreq, int type, int sampleFreq);
	int clockspeed;
	int samplefreq;
};


