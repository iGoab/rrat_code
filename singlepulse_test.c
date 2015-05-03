#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pulse_dm.h"

int main(int argc, char *argv[])
{
	struct data_alldm *data_adm;
	int i = 0, j = 0;
	int nr;
	char source[24];
	double sample_dt = 0.0f;
	double real_dm = 0.0f;
	int dm_flag = 0;
	
	/* Malloc the output data buf */
	struct output_dmdata *opd;
	opd = (struct output_dmdata *)malloc(sizeof(struct output_dmdata));
	
	/* Malloc the dedispersion data buf to process the dedispersion data */
	struct data_alldm *da;
	da = (struct data_alldm *)malloc(sizeof(struct data_alldm));
	struct output_dmdata *opd1;
	opd1 = (struct output_dmdata *)malloc(sizeof(struct output_dmdata));
	
	sprintf(opd->filename, argv[1]);
	printf("%s\n", opd->filename);
	opd->fp = open_file(opd->filename, "r");
	fread(source, sizeof(char), 24, opd->fp);
	
	fread(&sample_dt, sizeof(double), 1, opd->fp);
	fread(&(da->numpoint), sizeof(int), 1, opd->fp);
	fread(&(da->ndm), sizeof(int), 1, opd->fp);
	fread(&dm_flag, sizeof(int), 1, opd->fp);
	if (dm_flag) {
		fread(&(da->real_dm), sizeof(double), 1, opd->fp);
	}
	da->dm = (double *)malloc(sizeof(double) * da->ndm);
	fread(da->dm, sizeof(double), da->ndm, opd->fp);
	
	fclose(opd->fp);
	
	da->data_adm = (float *)malloc(sizeof(float) * da->numpoint);
	da->dm_data_max = 0.0f;
	da->dm_data_min = 0.0f;
	da->dm_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
	
	sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, *(da->dm + 0));
	opd->fp = open_file(opd->filename, "r");
	fread(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
	fread(&(da->dm_data_min), sizeof(float), 1, opd->fp);
	fread(&(da->dm_data_max), sizeof(float), 1, opd->fp);
	fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
	fclose(opd->fp);
	int spthresh = 4;
	int ncandsmax = 100000;
	int nsmax = 10;
	nr = da->numpoint;
	singlepulse(da,  spthresh, ncandsmax, nsmax, nr, sample_dt);
	
}
