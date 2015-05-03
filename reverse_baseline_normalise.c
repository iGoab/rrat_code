#include <stdio.h>

void reverseBaselineNormalise(float *data, unsigned long long sample0, int nsamps, int nchans) {
	FILE *file;
	unsigned long long sample;
	int isamp, ichan;
	float empty;
	unsigned long long sampstart, sampend;
	char string[256];
	double inmean;
	
	switch(nbits) {
		case 1:
			inmean = 0.5;
			break;
		case 2:
			inmean = 1.5;
			break;
		case 4:
			inmean = 7.5;
			break;
		case 8:
			inmean = 127.5;
			break;
		default:
			fprintf(stderr, "Can't use reverse baseline for nbits != 1, 2, 4, 8");
			return;
	}
	
	sample = sample0;
	if (baselineValidStart == )
}
