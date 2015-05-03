#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "psrfits.h"

int main(int argc, char *argv[])
{
	struct psrfits pf;
	//printf("%s\n", *(argv + optind));
	sprintf(pf.basefilename, *(argv + optind));
	printf("%s\n", pf.basefilename);
	psrfits_set_files(&pf, argc - optind, argv + optind);
	pf.tot_rows = pf.N = pf.T = pf.status = 0;
	//pf.hdr.chan_dm = 0.0;
	int rv = psrfits_open(&pf);
	if (rv) { fits_report_error(stderr, rv);  exit(1);}
	printf("1\n");
	printf("%lf\n", pf.hdr.MJD_epoch);
	FILE *fp;
	if ((fp = fopen("./data.txt", "w")) == NULL) {
		printf("Cannot create the file.\n");
		exit(0);
	}
	fprintf(fp, "pulsar_parameter:\n");
	fprintf(fp, "mjd: %ld \n", (int)pf.hdr.MJD_epoch);
	fprintf(fp, "fmjd: %f\n", fmod(pf.hdr.MJD_epoch, 1.0));
	fprintf(fp, "bw: %f\n", pf.hdr.df);
	fprintf(fp, "polyco_struct:\n");
	fprintf(fp, "pc.mjd: %d\n", (int)pf.hdr.MJD_epoch);
	fprintf(fp, "pc.fmjd: %f\n", fmod(pf.hdr.MJD_epoch, 1.0));
	fprintf(fp, "pc.nmin = 24 * 60\n");
	fprintf(fp, "pc.rphase = 0.0\n");
	fprintf(fp, "pc.nc = 1\n");
	fprintf(fp, "pc.rf = %f\n", pf.hdr.fctr);
	fprintf(fp, "pc.c[0] = 0.0\n");
	fprintf(fp, "pc.used = 0\n");
	
	fclose(fp);
}
