#include <stdio.h>
#include <stdlib.h>

int *dmshift(float f1, double df, int nchans, int nbands, double dm, double refrf, double tsamp, float frequency[])
{
	int i, cpb, *shift;
	float fi;
printf("%f\n", *(frequency));
	printf("1\n");
	shift = (int *)malloc(nchans * sizeof(int));
	fi = f1;
	printf("%f\n", fi);
	cpb = nchans/nbands;
	printf("%d %d %d\n", cpb, nchans, nbands);

	if (frequency[0] != 0.0) f1 = frequency[0];
	printf("%f\n", frequency[0]);
	for (i = 0; i < nchans; i++) {
		if (refrf > 0.0) f1 = refrf;
		printf("%d\n", i);
		if (frequency[0] != 0.0) fi = frequency[0];
 		shift[i]=(int)(dmdelay(fi,f1,dm)/tsamp);
		fi += df;
		if (!((i + 1) % cpb)) f1 += (double)cpb * df;
	}
	return shift;
}
