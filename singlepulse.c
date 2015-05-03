#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pulse_dm.h"

#define maxloopsize 1048576

void singlepulse(struct data_alldm *da, int spthresh, int ncandsmax, int nsmax, int nr, double tsamp)
// spthresh - single-pulse search threshold
// ncandsmax - maximum number of single-pulse candidates per DM channel
// nsmax - number of times to smooth time series for single-pulse search
// nr - the time series number of the data 
{
	int done_pulse, pulse;
	int loopsize, lsd1, llog, npulses;
	float realdata[nr];
	char best_pulses[132], hist_pulses[132], scrdisk[80];
	int dmidx;
	int ntim2;
	
	loopsize= nr/8;
	strncpy(scrdisk, "./", 2);
	lsd1 = 2;
	ntim2 = nr;
	dmidx = 1;
	/*
	 * Search for pulses in this DM file and write out to file best.tmp 
	 */
	 printf("%d\n", ntim2);
	 printf("%d\n", loopsize);
	 ntim2 = ntim2 - 10./tsamp;
	 printf("%d\n", ntim2);
	 done_pulse = pulse_(realdata, &spthresh, &nsmax, &dmidx, &ntim2,
						da->data_adm, &loopsize, scrdisk, &lsd1);
}
