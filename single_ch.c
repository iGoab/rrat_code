/*
 * This routine applies a threshold test on a time series. Matched 
 * boxcar filters are applied to the data, with different boxcar lengths,
 * and the threshold test is reapplied, in an effeort to find broadened
 * pulses.
 */
#define MAXPULSE 81920 /* Maximum number of pulses allowed in each DM channel */

#include <math.h>
#include <stdio.h>
#include "pulse.h"

void decimation(ndata, data) 
int *ndata;
float data[];
{
	int i;
	
	for (i = 0;i < *ndata; i+= 2) {
		data[i/2] = data[i];
	}
	
	*ndata = *ndata/2 + *ndata%2;
	return;
}

void single_ch(ndata, realdata, nsmax, thresh, ndm, scrdsk)
int ndata, nsmax;
char scrdsk[];
float realdata[], thresh;
{
	int i, ns, npulse, npoints;
	int ierr;
	float sample_int, width_1, width_2;
	FILE *fp;
	Pulsus pulse[MAXPULSE];
	
	thresh_1d(ndata, realdata, 0, &npulse, pulse, thresh, scrdsk);
	
	/* fprintf(stderr, "%d pulse found in dm channel %d smoothing %d\n", npulse, ndm, ns);*/
	
	if (npulse) {
		ierr = write_pulses(ns, npulse, pulse, ndm, scrdsk);
	}
	/* Note: all errors are written out in write_pulses. */
	npoints -= 1;
	decimation(&npoints, &realdata[0]);
	return;
}
