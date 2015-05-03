#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pulse_dm.h"
#include "/home/igoab/pulsar_software/pgplot/cpgplot.h"

float get_max_point(float var, float fix, int *point, int *int_point) {
	if ( var < fix)  {
		return(fix); 
	} else { 
		*point = *int_point; 
		return(var);
	}
}

float get_min_point(float var, float fix, int *point, int *int_point) {
	if ( var > fix)  {
		return(fix); 
	} else { 
		*point = *int_point; 
		return(var);
	}
}

int main(int argc, char *argv[])
{
	
	char source[24];
	double sample_dt = 0.0f;
	int float_flag = 0;
	int raw_flag = 0;
	float dm = 0.0f;
	int gauss_flag = 0;
	int avg_flag = 0;
	int i = 0, j = 0;
	int gray_flag = 0;
	int dm_flag = 0;
	int dm_map = 0;
	int draw_num_point = 5000;
	struct output_dmdata *opd = (struct output_dmdata *)malloc(sizeof(struct output_dmdata));
	struct data_alldm *da = (struct data_alldm *)malloc(sizeof(struct data_alldm ));
	struct smooth_data *fltsd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	struct smooth_data *gaussd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	
	
	if (argc < 2) {
		draw_pic_help();
		exit(0);
	}
	
	if (help_required(argv[1])) {
		draw_pic_help();
		exit(0);
	}
	i = 1;
	if (strings_equal(argv[i], "-f")) {
		float_flag = 1;
		dm = atof(argv[++i]);
		i++;
	} else if (strings_equal(argv[i], "-g")) {
		gauss_flag = 1;
		dm = atof(argv[++i]);
		i++;
	} else if (strings_equal(argv[i], "-a")) {
		avg_flag = 1;
		dm = atof(argv[++i]);
		i++;
	} else if (strings_equal(argv[i], "-r")) {
		raw_flag = 1;
		dm = atof(argv[++i]);
		i++;
	} else if (strings_equal(argv[i], "-y")) {
		gray_flag = 1;
		i++;
	} else if (strings_equal(argv[i], "-d")) {
		dm_map = 1;
		i++;
	}
	
	if (strings_equal(argv[i], "-n")) {
		draw_num_point = atoi(argv[++i]);
		i++;
	}
	//open the dm header file and read data from it to get ready for 
	//malloc data_buf
	sprintf(opd->filename, argv[i]);

	opd->fp = open_file(opd->filename, "r");
	
	fread(source, sizeof(char), 24, opd->fp);
	fread(&sample_dt, sizeof(double), 1, opd->fp);
	fread(&(da->numpoint), sizeof(int), 1, opd->fp);
	fread(&(da->ndm), sizeof(int), 1, opd->fp);
	fread(&dm_flag, sizeof(int), 1, opd->fp);
	if (dm_flag)
		fread(&(da->real_dm), sizeof(double), 1, opd->fp);
	da->dm = (double *)malloc(sizeof(double) * da->ndm);
	fread(da->dm, sizeof(double), da->ndm, opd->fp);
	
	fclose(opd->fp);
	


	//malloc data to get ready to plot gray map
	
	if (float_flag == 1) {
		fltsd->data_smooth = (float *)malloc(sizeof(float) * da->numpoint);
		fltsd->smooth_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
		sprintf(opd->filename, "./DmFluxFile/%s/%f.fltsd", source, dm);
		opd->fp = open_file(opd->filename, "r");
		fread(fltsd->smooth_x_coor, sizeof(float), da->numpoint, opd->fp);
		fread(&(fltsd->smooth_data_min), sizeof(float), 1, opd->fp);
		fread(&(fltsd->smooth_data_max), sizeof(float), 1, opd->fp);
		fread(fltsd->data_smooth, sizeof(float), da->numpoint, opd->fp);
		fclose(opd->fp);
		draw_flux_pic(fltsd->data_smooth, fltsd->smooth_x_coor, 0.0, da->numpoint, fltsd->smooth_data_min, fltsd->smooth_data_max, sample_dt, draw_num_point);
		free(fltsd->data_smooth);
		free(fltsd->smooth_x_coor);
	}
	
	if (gauss_flag == 1) {
		gaussd->data_smooth = (float *)malloc(sizeof(float) * da->numpoint);
		gaussd->smooth_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
		sprintf(opd->filename, "./DmFluxFile/%s/%f.gaussd", source, dm);
		opd->fp = open_file(opd->filename, "r");
		fread(gaussd->smooth_x_coor, sizeof(float), da->numpoint, opd->fp);
		fread(&(gaussd->smooth_data_min), sizeof(float), 1, opd->fp);
		fread(&(gaussd->smooth_data_max), sizeof(float), 1, opd->fp);
		fread(gaussd->data_smooth, sizeof(float), da->numpoint, opd->fp);
		fclose(opd->fp);
		draw_flux_pic(gaussd->data_smooth, gaussd->smooth_x_coor, 0.0, da->numpoint, gaussd->smooth_data_min, gaussd->smooth_data_max, sample_dt, draw_num_point);
		free(gaussd->data_smooth);
		free(gaussd->smooth_x_coor);
	}
	
	if (avg_flag == 1) {
		da->data_adm = (float *)malloc(sizeof(float) * da->numpoint);
	    da->dm_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
		sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, dm);
		opd->fp = open_file(opd->filename, "r");
		fread(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
		fread(&(da->dm_data_min), sizeof(float), 1, opd->fp);
		fread(&(da->dm_data_max), sizeof(float), 1, opd->fp);
		fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
		fclose(opd->fp);
		draw_flux_pic(da->data_adm, da->dm_x_coor, 0.0, da->numpoint, da->dm_data_min, da->dm_data_max, sample_dt, draw_num_point);
		free(da->data_adm);
		free(da->dm_x_coor);
	}
	
	if (raw_flag == 1) {
		da->data_adm = (float *)malloc(sizeof(float) * da->numpoint);
		da->dm_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
		sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", source, dm);
		opd->fp = open_file(opd->filename, "r");
		fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
		set_rawdata(da, da->numpoint, sample_dt);
		fclose(opd->fp);
		draw_flux_pic(da->data_adm, da->dm_x_coor, 0.0, da->numpoint, da->dm_data_min, da->dm_data_max, sample_dt, draw_num_point);
		free(da->data_adm);
		free(da->dm_x_coor);
	}
	
	if (gray_flag == 1) {
		da->data_adm = (float *)malloc(sizeof(float) * da->numpoint);
		da->dm_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
		float *gray_data = (float *)malloc(sizeof(float) * da->numpoint * da->ndm);
		for (i = 0; i < da->ndm; i++) {
			printf("%f\n", *(da->dm + i));
			sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, *(da->dm + i));
			opd->fp = open_file(opd->filename, "r");
			fread(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
			fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
			memcpy(gray_data + i * da->numpoint, da->data_adm, da->numpoint);
			for (j = 0; j < da->numpoint; j++) *(gray_data + j + i * da->ndm) = *(da->data_adm + j);
			fclose(opd->fp);
		}
		draw_gray_map(gray_data, da->numpoint, da->ndm, sample_dt);
		free(da->data_adm);
		free(gray_data);
	}
	if (dm_flag) {
		if (dm_map == 1) {
			float *diff_dm_flux;
			float realdm_max = 0.0f, realdm_min = 0.0f;
			float otherdm_max = 0.0f;
			int realdm_max_samp = 0, realdm_min_samp = 0;
			int otherdm_max_samp = 0;
			da->data_adm = (float *)malloc(sizeof(float) * da->numpoint);
			da->dm_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
			diff_dm_flux = (float *)malloc(sizeof(float) * da->ndm);
			printf("%f\n", da->real_dm);
			sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, da->real_dm);
			opd->fp = open_file(opd->filename, "r");
			fread(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
			fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
			realdm_max = *(da->data_adm + 0);
			fclose(opd->fp);
			*(da->data_adm + 1) = 0.0f;
			*(da->data_adm + da->numpoint - 1) = 0.0f;
			*(da->data_adm + da->numpoint - 2) = 0.0f;
			printf("%d\n", da->numpoint);
			for (i = 1; i < da->numpoint; i++) {
				realdm_max = get_max_point(da->data_adm[i], realdm_max, &realdm_max_samp, &i);
			}
			for (i = 1; i < da->numpoint; i++) {
				realdm_min = get_min_point(da->data_adm[i], realdm_min, &realdm_min_samp, &i);
			}
			printf("%d\n", realdm_max_samp);
			printf("%f\n", realdm_max);
				printf("%d\n", realdm_min_samp);
			printf("%f\n", realdm_min);
			for (i = 0; i < da->ndm; i++) {
				printf("\rProces %f.dav", *(da->dm + i));
				fflush(stdout);
				if(*(da->dm + i) == da->real_dm) { *(diff_dm_flux + i) = realdm_max; continue; }
				sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, *(da->dm + i));
				opd->fp = open_file(opd->filename, "r");
				fread(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
				fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
				otherdm_max = *(da->data_adm + realdm_max_samp);
				otherdm_max_samp = realdm_max_samp;
				for (j = realdm_max_samp - 1000; j <= realdm_max_samp + 1000; j++)
					otherdm_max = get_max_point(*(da->data_adm + j), otherdm_max, &otherdm_max_samp, &j);
				*(diff_dm_flux + i) = otherdm_max;
				fclose(opd->fp);
			}
			for (i = 0; i < da->ndm; i++) {
				printf("%f\n", *(diff_dm_flux + i));
			}
			if (cpgbeg(0, "/xs", 1, 1) != 1) {
				return EXIT_FAILURE;
			}
			float x[da->ndm];
			for (i = 0; i < da->ndm; i++) x[i] = i;
			cpgenv(0.0, da->ndm, 0.0, 1.0, 0, 0);
			cpgline(da->ndm, x, diff_dm_flux);
			cpgend();
			
		}
	} else {
		printf("no real dm input.\n");
	}
}
