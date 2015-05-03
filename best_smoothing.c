#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define maxcands pow(2, 21) 
void best_smoothing(int ncandsmax, int dmidx, length, char *best_pulses,
					char *hist_pulses, thresh, nsmax, 
					scrdsk1, npulses, double tsamp)
{
	int ndm, ns, timebin, ncandsmax, dmidx, npulses;
	float snr, mean, rms, dm;
	double tsamp, time, width;
	int ndm1, ns1, timebin1, npulse, n_to_print;
	float snr1, mean1, rms1, thresh;
	
	int nsmax;
	
	int length, power, lun;
	
	float rlength;
	
	int best_dm[maxcands], best_ns[maxcands], best_time[maxcands], best_snr[maxcands];
	
	int itolerance, io_error, isamp;
	
	char inline1[80], inline2[80];
	char scrdsk1[132];
	
	int ind;
	int latest = 0;
	
	int lnblnk;
	int lunin, lunout;
	
	int nit, nitmax = 4;
	
	int l0, l1, l2, l3, l4, l5, i;
	npulses = 0;
}
