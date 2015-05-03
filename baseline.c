#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pulse_dm.h"

#define nrav 32
int min(int a, int b)
{
	if (a > b) {return b;}
	else {return a;}
}
/*****************************************************************
 * subtracts a sloping baseline from the time series contained in*
 * the array series(ntim). Writes out the mean which is then   *
 * subracted 
 * The time series is then normalized so that the rms is unity *
 * *************************************************************/
void baseline(struct data_alldm *da, double dt, int nr)
{
	int i, j, nrun;
	
	float mean, sums, sumd,var, mea, msq, rms;
	float rssq[nrav], rsum[nrav];
	float x[nrav], y[nrav], z[nrav];
	float e[nrav], slope, inter, eslo, eint;
	mean = 0.0f;
	sums = 0.0f;
	sumd = 0.0f;
	var = 0.0f;
	mea = 0.0f;
	msq = 0.0f;
	rms = 0.0f;
	slope = 0.0f;
	inter = 0.0f;
	eslo = 0.0f;
	eint = 0.0f;
	for (i = 0; i < nrav; i++) {
		rsum[i] = 0.0f;
		rssq[i] = 0.0f;
		z[i] = 0.0f;
	}
	nrun = nr/nrav;
	for (i = 0; i < nr; i++) {
		j = min(nrav, (i + 1)/nrun + 1);
		rsum[j - 1] = rsum[j - 1] + da->data_adm[i];
		rssq[j - 1] = rssq[j - 1] + da->data_adm[i] * da->data_adm[i];
		sumd = sumd + da->data_adm[i];
	}

	mean = sumd/(float)nr;
	for (i = 0; i < nrav; i++) {
		mea = rsum[i]/(float)(nrun);
		msq = rssq[i]/(float)(nrun);
		var = sqrt(msq - mea*mea);
		j = nrun/2 + (i) * nrun;
		x[i] = i + 1;
		y[i] = mea;
		e[i] = var/sqrt((float)(nrun));
	}
	slfit(x, y, z, nrav, e, 0, &inter, &slope, &eint, &eslo);
	printf("Subtracting mean from function: %f %f\n", inter, slope); 
	for (i = 0; i < nr; i++) {
		mean = inter + slope * ((float)(i + 1)) * ((float)nrav) / ((float)nr);
		da->data_adm[i] = da->data_adm[i] - mean;
		
	}
	sums = 0.0f;
	for (i = 0; i < nr; i++) {
		sums = sums + da->data_adm[i] * da->data_adm[i];
	}
	
	rms = sqrt(sums/((float)nr));
	for (i = 0; i < nr; i++) {
		da->data_adm[i] = da->data_adm[i]/rms;
		if (da->dm_data_max < *(da->data_adm + i)) {
			da->dm_data_max = *(da->data_adm + i);
		}
		if (da->dm_data_min > *(da->data_adm + i)) {
			da->dm_data_min = *(da->data_adm + i);
		}
		
		*(da->dm_x_coor + i) = (float)i * dt;
	}
	printf("Time series now has zero mean and unit rms\n");
}
