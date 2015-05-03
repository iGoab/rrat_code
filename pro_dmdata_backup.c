#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pulse_dm.h"

int main(int argc, char *argv[])
{
	struct data_alldm *data_adm;
	int i = 0, j = 0;
	int float_flag = 0;
	int gauss_flag = 0;
	char message[256];
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
	
	/* Malloc the float smooth data buf to store smooth data */
	struct smooth_data *fltsd;
	fltsd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	fltsd->dura_point = 15;
	fltsd->smooth_data_min = 0.0f;
	fltsd->smooth_data_max = 0.0f;
	fltsd->smooth_sum = 0.0f;
	fltsd->smooth_count = 0;
	
	struct smooth_data *gaussd;
	gaussd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	gaussd->sigma = 10;
	gaussd->smooth_data_min = 0.0f;
	gaussd->smooth_data_max = 0.0f;
	gaussd->smooth_sum = 0.0f;
	gaussd->smooth_count = 0;
	
	
	if (argc < 2) {
		pro_dmdata_help();
		exit(0);
	}
	
	if (help_required(argv[1])) {
		pro_dmdata_help();
		exit(0);
	}
	
	i = 1;


	if (strings_equal(argv[i], "-f")) {
		float_flag = 1;
		fltsd->dura_point = atoi(argv[++i]);
		gauss_flag = 0;
		i++;
	} else if (strings_equal(argv[i], "-g")) {
		gauss_flag = 1;
		float_flag = 0;
		gaussd->sigma = atoi(argv[++i]);
		i++;
	} 

	gaussd->dura_point = 2 * (int)rint(gaussd->sigma * 3.5) + 1;
	gaussd->gauss_coe = (float *)malloc(sizeof(float) * gaussd->dura_point);
	sprintf(opd->filename, argv[i]);
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
	fltsd->data_smooth = (float *)malloc(sizeof(float) * da->numpoint);
	fltsd->smooth_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
	gaussd->data_smooth = (float *)malloc(sizeof(float) * da->numpoint);
	gaussd->smooth_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
	printf("%d\n", gaussd->sigma);
	printf("%d\n", gaussd->dura_point);
	printf("%d\n", fltsd->dura_point);
	

	printf("float: %d, gauss: %d\n", float_flag, gauss_flag);
	for (j = 0; j < da->ndm; j++) {
		sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", source, *(da->dm + j));
		printf("\rProcess %lf dm file", *(da->dm + j));
		fflush(stdout);
		opd->fp = open_file(opd->filename, "r");
		fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
		fclose(opd->fp);
		calc_data_sub_avg(da, sample_dt, da->numpoint);	
		sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, *(da->dm + j));
		//for (j = 0; j < da->numpoint; j++) *(da->dm_x_coor + j) = j;
		//printf("1\n");
		
		opd->fp = open_file(opd->filename, "wb");
		fwrite(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
		fwrite(&(da->dm_data_min), sizeof(float), 1, opd->fp);
		fwrite(&(da->dm_data_max), sizeof(float), 1, opd->fp);
		fwrite(da->data_adm, sizeof(float), da->numpoint, opd->fp);
		fclose(opd->fp);
		
		if (float_flag == 1) {
			calc_float_avg_smooth(fltsd, da, sample_dt, da->numpoint);
			sprintf(opd->filename, "./DmFluxFile/%s/%f.fltsd", source, *(da->dm + j));
			opd->fp = open_file(opd->filename, "wb");
			fwrite(fltsd->smooth_x_coor, sizeof(float), da->numpoint, opd->fp);
			fwrite(&(fltsd->smooth_data_min), sizeof(float), 1, opd->fp);
			fwrite(&(fltsd->smooth_data_max), sizeof(float), 1, opd->fp);
			fwrite(fltsd->data_smooth, sizeof(float), da->numpoint, opd->fp);
			fclose(opd->fp);
		}
		
		if (gauss_flag == 1) {
			calc_gauss_smooth(gaussd, da, sample_dt, da->numpoint);
			sprintf(opd->filename, "./DmFluxFile/%s/%f.gaussd", source, *(da->dm + j));
			opd->fp = open_file(opd->filename, "wb");
			fwrite(gaussd->smooth_x_coor, sizeof(float), da->numpoint, opd->fp);
			fwrite(&(gaussd->smooth_data_min), sizeof(float), 1, opd->fp);
			fwrite(&(gaussd->smooth_data_max), sizeof(float), 1, opd->fp);
			fwrite(gaussd->data_smooth, sizeof(float), da->numpoint, opd->fp);
			fclose(opd->fp);
		}
	}
	//draw_flux_pic(da->data_adm, da->dm_x_coor, 0.0, da->numpoint, da->dm_data_min, da->dm_data_max, sample_dt);
	printf("\n");
}
