#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include "/home/igoab/pulsar_software/pgplot/cpgplot.h"
#include "./pulse_dm.h"
#include "psrfits.h"

#define PI M_PI

int main(int argc, char *argv[])
{
	static struct option long_opts[] = {
		{"output", 1, NULL, 'o'},
		{"dedispersion", 1, NULL, 'd'},	
		{"read", 0, NULL, 'r'},
		{"zerodm", 0, NULL, 'z'},
		{"help", 0, NULL, 'h'},
		{0, 0, 0, 0}
	};
	
	int opt, opti;
	int i = 0, j = 0, k = 0;
	double dm = 0.0;
	char output_base[256] = "\0";
	int read_file = 0;
	int apply_zerodm = 0;
	while((opt = getopt_long(argc, argv, "o:d:rzh", long_opts, &opti)) != -1) {
		switch(opt) {
			case 'o':
				strncpy(output_base, optarg, 255);
				output_base[255] = '\0';
				break;
			case 'd':
				dm = atof(optarg);
				break;
			case 'r':
				read_file = 1;
				break;	
			case 'z':
				apply_zerodm = 1;
				break;
			case 'h':
			default:
				usage();
				exit(0);
				break;
		}
	}
	if (optind == argc) {
		usage(0);
		exit(1);
	}
	
	/**************************************************
	 * Print the header of the file and the header of *
	 * the subintegration							  *
	 * ************************************************/
	if (read_file) {
		fitsfile *fptr;
		char card[FLEN_CARD];
		int status = 0, nkeys, ii;
		fits_open_file(&fptr, *(argv + optind), READONLY, &status);
		fits_get_hdrspace(fptr, &nkeys, NULL, &status);
		printf("*************************************************\n");
		printf("FILE_HEADER");
		for (ii = 0; ii < nkeys; ii++) {
			fits_read_record(fptr, ii, card, &status);
			printf("%s\n", card);
		} 
		printf("END\n");
		fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
		printf("***************************************************\n");
		printf("SUBINT:");
		fits_get_hdrspace(fptr, &nkeys, NULL, &status);
		for (ii = 0; ii < nkeys; ii++) {
			fits_read_record(fptr, ii, card, &status);
			printf("%s\n", card);
		}
		printf("END\n");
		printf("*****************************************************\n");
		fits_close_file(fptr, &status);
		
		if (status) {
			fits_report_error(stderr, status);
		}
		exit(1);
	}
	
	struct psrfits pf;
	psrfits_set_files(&pf, argc - optind, argv + optind);
	pf.tot_rows = pf.N = pf.T = pf.status = 0;
	pf.hdr.chan_dm = 0;
	
	int rv = psrfits_open(&pf);
	
	if ( rv ) { fits_report_error(stderr, rv); exit(1);}
	int nr = pf.hdr.nsblk * pf.rows_per_file;

	//float data_max = 0.0f;
	float data_nx = (float) nr;
	float dm_step = 1.0f;
	float dm_max = 100.0f;
	float *data;
	float *data_zerodm;
	double *chan_delays;
	//dm_step = calc_dm_step(&pf);
	//printf("tsamp: %f fctr: %f BW: %f\n", pf.hdr.dt, pf.hdr.fctr, pf.hdr.BW);
	//printf("%f\n", dm_step);
	/* Alloc data buffers */
	pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_offsets = (float *)malloc(sizeof(float) 
			* pf.hdr.nchan * pf.hdr.npol);
	pf.sub.dat_scales = (float *)malloc(sizeof(float) 
			* pf.hdr.nchan * pf.hdr.npol);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);	
	data = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf.hdr.nchan);
	data_zerodm = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf.hdr.nchan);
	chan_delays = (double *)malloc(sizeof(double) * pf.hdr.nchan);		
	
	int ndm = (int)(dm_max/dm_step) + 1;
	/* Malloc dedispersion data */
	struct pulse_draw *pd;
	pd = (struct pulse_draw *)malloc(sizeof(struct pulse_draw));
	struct data_alldm *da;
	pd->da = (struct data_alldm *)malloc(sizeof(struct data_alldm ) * ndm);
	for (i = 0; i < ndm; i++) {
		if (i == ndm - 1) {
			pd->da[i].dm = dm;
		} else {
			pd->da[i].dm = i * dm_step;
		}
		pd->da[i].numpoint = nr;
		pd->da[i].avg_flux = 0.0f;
		pd->da[i].dm_data_max = 0.0f;
		pd->da[i].dm_data_min = 0.0f;
		pd->da[i].idelays = (int *)malloc(sizeof(int) * pf.hdr.nchan);
		pd->da[i].data_adm  = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf.rows_per_file);	
		pd->da[i].dm_x_coor = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf. rows_per_file);
	}
	printf("%d\n", pf.hdr.nsblk);
	/* use straight insertion to sort dm */
	float a;
	for (i = 1; i < ndm; i++) {
		a = pd->da[i].dm;
		j = i;
		while(i > 0 && pd->da[i - 1].dm > a) {
			pd->da[i].dm = pd->da[i - 1].dm;
			i--;
		}
		pd->da[i].dm = a;
	}
	/* Malloc the fold_buf data to calculate the fold */
	/************************************************/
	struct fold_buf *fb;
	fb = (struct fold_buf *)malloc(sizeof(struct fold_buf));
	fb->pulsar_freq = 1.399541538720;
	fb->freq_point = (int)((1/fb->pulsar_freq) / pf.hdr.dt);
	//fb->fold_data = (float *)malloc(sizeof(float) * fb->freq_point);
	//fb->fold_x_coor = (float *)malloc(sizeof(float) * fb->freq_point);
	fb->fold_data_max = 0.0f;
	fb->fold_data_min = 0.0f;
	fb->nbin = 256;
	
	/*************************************************/
	
	/* Malloc the smooth_data to smooth the data */
	struct smooth_data *sd;
	sd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	sd->dura_point = 200;
	sd->sum_point = (int)(nr / sd->dura_point);
	sd->smooth_count = 0;
	sd->smooth_data_max = 0.0f;
	sd->smooth_data_min = 0.0f;
	sd->smooth_sum = 0.0f;
	
	/* Malloc the float smooth data to store after float smooth */
	struct smooth_data *fltsd;
	fltsd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	fltsd->dura_point = 151;   // it must be odd
	fltsd->smooth_data_min = 0.0f;
	fltsd->smooth_data_max = 0.0f;
	fltsd->smooth_sum = 0.0f;
	fltsd->smooth_count = 0;
	fltsd->data_smooth = (float *)malloc(sizeof(float) * nr);
	fltsd->smooth_x_coor = (float *)malloc(sizeof(float) * nr);
	
	/* Malloc the gauss smooth data to store after gauss smooth */
	struct smooth_data *gaussd;
	gaussd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	gaussd->sigma = 100;
	gaussd->dura_point = 2 * (int)rint(gaussd->sigma * 3.5) + 1;
	gaussd->smooth_data_min = 0.0f;
	gaussd->smooth_data_max = 0.0f;
	gaussd->smooth_sum = 0.0f;
	gaussd->smooth_count = 0;
	gaussd->data_smooth = (float *)malloc(sizeof(float) * nr);
	gaussd->smooth_x_coor = (float *)malloc(sizeof(float) * nr);
	gaussd->gauss_coe = (float *)malloc(sizeof(float) * gaussd->dura_point);
	
	/* Malloc the exponential smooth data to store after exponential smooth */
	/* The bellow is just double exponential smoothing */
	struct smooth_data *double_expsd;
	double_expsd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	double_expsd->exp_alpha = 0.26;
	double_expsd->exp_beta = 0.24;
	double_expsd->smooth_data_max = 0.0f;
	double_expsd->smooth_data_min = 0.0f;
	double_expsd->smooth_count = 0;
	double_expsd->data_smooth = (float *)malloc(sizeof(float) * nr);
	double_expsd->data_trend = (float *)malloc(sizeof(float) * nr);
	double_expsd->smooth_x_coor = (float *)malloc(sizeof(float) * nr);
	
	/* Malloc the single_pulse_map to draw the picture */
	//struct pulse_num *pn;
	pd->pn = (struct pulse_num *)malloc(sizeof(struct pulse_num) * ndm);
	for (i = 0; i < ndm; i++) {
		pd->pn[i].pulse_count = 0;
		pd->pn[i].rms = 0.0f;
		pd->pn[i].count_threshold = 0;
		pd->pn[i].threshold_value = 0.0f;
		pd->pn[i].pulse_peak_max = 0.0f;
		pd->pn[i].pulse_peak_min = 0.0f;
	}

	/* Malloc the output data file pointer */
	struct out_put_data *opd;
	opd = (struct out_put_data *)malloc(sizeof(struct out_put_data));
	sprintf(opd->filename, "data.dat");
	
	int run = 1, row_count = 0, chan_count = 0;
	
	int calc_data_dm = 1;
	int flag_check_freq = 0;
	float sum = 0.0f, avg = 0.0f;
	if ((strncmp(pf.hdr.poln_order, "AABBCRCI", 8) == 0) || 
			  (strncmp(pf.hdr.poln_order, "AABB", 4) == 4)) {
		while(run) {
			rv = psrfits_read_subint(&pf);
			if (rv) {
				if (rv == FILE_NOT_OPENED) rv = 0;
				run = 0;
				break;
			}
			
			if (calc_data_dm) {
				calc_data_dm = 0;
				for (i = 0; i < ndm; i++) {
					flag_check_freq = calc_dm(&pf, chan_delays, pd->da[i].idelays, pd->da[i].dm);	
				}	
			}
			
			get_AABB_I(&pf, data);	
			for (i = 0; i < ndm; i++) {
				delay_data(&pf, data, &(pd->da[i]), row_count, flag_check_freq, apply_zerodm, i);
			}
			row_count++;
			printf("\rProcess subint(file %d row %d/%d)",
							pf.filenum, pf.rownum - 1, pf.rows_per_file);
		} 	
	} else if ((strncmp(pf.hdr.poln_order, "IQUV", 4) == 0) ||
				  (strncmp(pf.hdr.poln_order, "AA", 2) == 0)) {
		while(run) {
			rv = psrfits_read_subint(&pf);
			if (rv) {
				if (rv == FILE_NOT_OPENED) rv = 0;
				run = 0;
				break;
			}
			
			if (calc_data_dm) {
				calc_data_dm = 0;
				for (i = 0; i < ndm; i++) {
					flag_check_freq = calc_dm(&pf, chan_delays, pd->da[i].idelays, pd->da[i].dm);			
				}	
			}
			//for(j = 0; j < ndm; j++)
			get_only_I(&pf, data);
			for (i = 0; i < ndm; i++) {
				delay_data(&pf, data, &(pd->da[i]), row_count, flag_check_freq, apply_zerodm, i);
			}
			//}
			//for (i = 0; i < pf.hdr.nsblk; i++)
					//	printf("%f\n", *(da[9].data_adm + i + row_count * pf.rows_per_file));
			row_count++;
			printf("\rProcess subint(file %d row %d/%d)",
							pf.filenum, pf.rownum - 1, pf.rows_per_file);

		}
	}
	/*for (i = 0; i < pf.hdr.nsblk; i++)
					for (j = 0; j < pf.hdr.nchan; j++)
						printf("%f\n", *(da[9].data_adm + j + i * pf.hdr.nchan));*/
	//output_datatotxt(&pf, opd, &da[9]);
	FILE *dm_flux_fp;
	for (i = 0; i < ndm; i++) 
		//set_rawdata(&pf, &(pd->da[i]), nr);
		calc_data_sub_avg(&pf, &(pd->da[i]), nr);
	char dm_flux_file[256];
	struct stat buf;
	/*for (i = 0; i < ndm; i++) {

		if (stat("./DmFluxFile", &buf) == -1) {
			system("mkdir DmFluxFile");
		}
	
		sprintf(dm_flux_file, "./DmFluxFile/dm%.3f_fluxfile", pd->da[i].dm);
		if ((dm_flux_fp = fopen(dm_flux_file, "wb")) == NULL) {
			fprintf(stderr, "Cannot create the file\n");
			exit(1);
		}
		fprintf(dm_flux_fp, "number points: %d\n", nr);
		fprintf(dm_flux_fp, "the value of dm: %f\n", pd->da[i].dm);
		fprintf(dm_flux_fp, "the max flux of the data: %f\n", pd->da[i].dm_data_max);
		fprintf(dm_flux_fp, "the min flux of the data: %f\n", pd->da[i].dm_data_min);
		fprintf(dm_flux_fp, "the average flux of the data: %f\n", pd->da[i].avg_flux);
		fprintf(dm_flux_fp, "dm_x_coor\tdata_value\n");
		for (j = 0; j < nr; j++) {
			fprintf(dm_flux_fp, "%f\t%s%f\n", *(pd->da[i].dm_x_coor + j), *(pd->da[i].data_adm + j) <0?"":" ", 
												*(pd->da[i].data_adm + j));
		}

		fclose(dm_flux_fp);
		sprintf(dm_flux_file, "");
		
	}*/
	//set_rawdata(&da[9], nr);
	//calc_fold(fb, &da[9], nr);
	
	//output_datatotxt(&pf, opd, &da[9]);
	//The code to smooth the data
	//calc_smooth(&pf, sd, &da[9], nr);
	//calc_fold_smooth(fb, sd);
	//calc_smooth_fold(sd, fb);
	for (i = 0; i < ndm; i++) {
		calc_float_avg_smooth(&pf, fltsd, &(pd->da[i]), nr);
		//printf("%d\n", i);
		calc_pulse_num(&pf, &(pd->pn[i]), fltsd, nr);
		free(pd->da[i].dm_x_coor);
		free(pd->da[i].idelays);
		free(pd->da[i].data_adm);
		//printf("%d\n", i);
	}
	//calc_pulse_num(&pf, pn, fltsd, nr);
	//calc_gauss_smooth(&pf, gaussd, &da[9], nr);
	
	//calc_dbexp_smooth(double_expsd, &da[9], nr);
	//calc_smooth_fold(&pf, fltsd, fb);
	//calc_smooth_fold(&pf, gaussd, fb);
	
	//fold_phase(&pf, &da[9], fb, nr);
	//printf("threshold: %f value_count: %d\n", da[9].rms, count_threshold);
	//printf("%f \n", da[9].dm_data_min);
	//printf("%f %f\n", fltsd->smooth_data_min, fltsd->smooth_data_max);
	//draw_flux_pic(&pf, pd->da[0].data_adm, pd->da[0].dm_x_coor, 0.0, nr, pd->da[0].dm_data_min, pd->da[0].dm_data_max);
	
	//draw_flux_pic(&pf, sd->data_smooth, sd->smooth_x_coor, 0.0, sd->smooth_count, sd->smooth_data_min, sd->smooth_data_max);

	//draw_flux_pic(&pf, fltsd->data_smooth, fltsd->smooth_x_coor, 0.0, fltsd->smooth_count, fltsd->smooth_data_min, fltsd->smooth_data_max);

	//draw_flux_pic(double_expsd->data_smooth, double_expsd->smooth_x_coor, 0.0, double_expsd->smooth_count, double_expsd->smooth_data_min, double_expsd->smooth_data_max);
	//draw_flux_pic(&pf, gaussd->data_smooth, gaussd->smooth_x_coor, 0.0, gaussd->smooth_count, gaussd->smooth_data_min, gaussd->smooth_data_max);
	
	//draw_flux_pic(&pf, fb->fold_data, fb->fold_x_coor, 0.0, fb->freq_point, fb->fold_data_min, fb->fold_data_max);
	//draw_histo_pic(pn);
	
	//draw_dmpic(&pf, pd, dm_max, ndm, nr);
	draw_flux_dm_sample(&pf, pd, dm_max, ndm, nr);
	//draw_fb_pulse_dis(&pf, fb, pn);
	psrfits_close(&pf);
	if (rv > 100) {fits_report_error(stderr, rv); }
	exit(0);
}

