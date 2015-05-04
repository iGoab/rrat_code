#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>
#include <getopt.h>
#include "/home/pulsar_software/pgplot/cpgplot.h"
#include "psrfits.h"
#include "pulse_dm.h"

#define LINE 1024

struct rrat {
	char name[24];
	char p[8];
	char pdot[8];
	char dm[8];
	char ra[16];
	char dec[16];
	char l[8];
	char b[8];
	char rate[8];
	char logB[8];
	char logts[8];
	char dhat[8];
	char fluxd[8];
	char pulse_width[8];
	char survey[24];
};

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

void get_axis_mm(float *x, float *y, float *x_min, float *x_max, float *y_min, float *y_max, int num)
{
	int i = 0;
	*x_min = x[0];
	*x_max = x[num - 1];
	*y_min = y[0];
	*y_max = y[0];
	for (i = 1; i < num; i++) {
		if (*y_min > *(y + i)) *y_min = *(y + i);
		if (*y_max < *(y + i)) *y_max = *(y + i);
	}
}

void get_raw_chan(struct psrfits *pf, float *raw_chan, int num)
{
	struct hdrinfo	*hdr = &(pf->hdr);
	struct subint *sub = &(pf->sub);
	unsigned char *data, *data1;
	
	int i = 0, j = 0;
	if ((strncmp (hdr->poln_order, "AABBCRCI", 8) == 0) ||
					(strncmp(hdr->poln_order, "AABB", 4) == 0)) {
		for (i = 0; i < hdr->nsblk; i++) {
			data = sub->rawdata + i * hdr->nchan * hdr->npol;
			data1 = data + hdr->nchan;
			for (j = 0; j < hdr->nchan; j++, data++, data1++)
				*(raw_chan + j + i * hdr->nchan + num * hdr->nsblk * hdr->nchan) = (float)(*data + *data1);
		}
	} else if ((strncmp(hdr->poln_order, "IQUV", 4) == 0) ||
					(strncmp(hdr->poln_order, "AA", 2) == 0) ||
					(strncmp(hdr->poln_order, "AA+BB", 5) == 0)) {
		for (i = 0; i < hdr->nsblk; i++) {
			data = sub->rawdata + i * hdr->nchan * hdr->npol;
			for (j = 0; j < hdr->nchan; j++, data++) {
				*(raw_chan + j + i * hdr->nchan + num * hdr->nsblk * hdr->nchan) = (float)*data;
			}	
		}	
	}
}
void get_plot_chan(struct psrfits *pf, float *raw_chan, float *plot_chan, 
						int left_row, int get_row_point, int num)
{
	struct hdrinfo *hdr = &(pf->hdr);

	int i = 0, j = 0, k = 0;
	k = left_row * hdr->nsblk + get_row_point;
	printf("%d %d %d %d\n", (k - num/2) * hdr->nchan, (k + num/2) * hdr->nchan, k, hdr->nsblk);
	for (i = 0; i < num; i++) {
		for (j = 0; j < hdr->nchan; j++) {
			*(plot_chan + i + j * num) = *(raw_chan + j + (k - num/2 + i) * hdr->nchan);
		}
	}
	
}

int main(int argc, char *argv[])
{
	static struct option long_opts[] = {
		{"read", 0, NULL, 'r'},
		{"dm_header", 1, NULL, 'd'},
		{"thresh", 1, NULL, 't'},
		{"sample", 1, NULL, 's'},
		{"help", 0, NULL, 'h'},
		{0, 0, 0, 0}
	};
	int opt, opti;
	struct rrat rr;
	char source[24];
	char dmfile_name[256];
	char fitsfile_name[256];
	double sample_dt = 0.0f;
	int float_flag = 0;
	int raw_flag = 0;
	float dm = 0.0f;
	float P = 0.0f;
	int downsample = 0;
	float downsample_time = 0.0f;
	float pw = 0.0f;
	int gauss_flag = 0;
	int avg_flag = 0;
	int i = 0, j = 0, k = 0, l = 0;
	int gray_flag = 0;
	int dm_flag = 0;
	int read_dm = 0;
	int read_head = 0;
	int set_thresh = 9;
	char rrat_source_dir[64] = "/home/igoab/RRAT_code/rrat_source.dat";
	struct psrfits pf;
	fitsfile *fptr;
	int status;
	int downsample_flag = 1;
	struct output_dmdata *opd = (struct output_dmdata *)malloc(sizeof(struct output_dmdata));
	struct data_alldm *da = (struct data_alldm *)malloc(sizeof(struct data_alldm ));
	struct smooth_data *fltsd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	struct smooth_data *gaussd = (struct smooth_data *)malloc(sizeof(struct smooth_data));
	
	while(((opt) = getopt_long(argc, argv, "d:t:s:rh", long_opts, &opti)) != -1) {
		switch(opt) {
			case 'r':
				read_head = 1;
				break;
			case 'd':
				read_dm = 1;
				strncpy(opd->filename, optarg, 255);
				opd->filename[255] = '\0';
				break;
			case 't':
				set_thresh = atoi(optarg);
				break;
			case 's':
				downsample = atoi(optarg);
				downsample_flag = 0;
				break;
			case 'h':
			default:
				find_peak_help();
				exit(0);
				break;
		}
	}
	/*if (argc < 2) {
		find_peak_help();
		exit(0);
	}
	opt = 1;
	if (help_required(argv[opt])) {
		find_peak_help();
		exit(0);
	}
	
	if (strings_equal(argv[opt], "-r")) {
		read_head = 1;
		strncpy(fitsfile_name, argv[++opt], 255);
	} else {		
		strncpy(fitsfile_name, argv[opt++], 255);
	}
	if (strings_equal(argv[opt], "-d")) {
		read_dm = 1;
		strncpy(opd->filename, argv[++opt], 255);
	}
	
	if (strings_equal(argv[++opt], "-t")) {
		set_thresh = atoi(argv[++opt]);
	}
	
	if (opt + 1 <  argc) {
		
		if (strings_equal(argv[++opt], "-s")) {
			downsample= atoi(argv[++opt]);
			downsample_flag = 0;
		}
	}*/
	if (optind == argc) {
		find_peak_help();
		exit(1);
	}
	printf("set_thresh: %d\n", set_thresh);
	
	/***************/
	/* the bottome code is to read the file header */
	/***************/
	if (read_head == 1) {
		fits_open_file(&fptr, *(argv+optind), READONLY, &status);
		char card[FLEN_CARD];
		int nkeys, ii;
		fits_get_hdrspace(fptr, &nkeys, NULL, &status);
		puts("**************************************************");
		puts("FILE_HEADER");
		for (ii = 0; ii < nkeys; ii++) {
			fits_read_record(fptr, ii, card, &status);
			printf("%s\n", card);
		}
		printf("END\n");
		fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
		puts("**************************************************");
		puts("SUBINT_HEADER:");
		fits_get_hdrspace(fptr, &nkeys, NULL, &status);
		for (ii = 0; ii < nkeys; ii++) {
			fits_read_record(fptr, ii, card, &status);
			printf("%s\n", card);
		}
		puts("END");
		puts("**************************************************");
		fits_close_file(fptr, &status);
		print_fits_error(status);
		exit(0);
	}
	
	/*************/
	/* the bottom code is to find the peak of the data */
	/*************/
	if (read_dm == 1) {
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
		
		//printf("%s\n", source);
		
		/**********************************/
		/* the bottom code is to read the rrat_source file and get the parameter */
		FILE *fp;
		if ((fp = fopen(rrat_source_dir, "r")) == NULL) {
			fprintf(stderr, "Can not read the file\n");
			exit(1);
		}
		
		char *buf;
		int num;
		int source_exist = 0;
		buf = (char *)malloc(sizeof(char) * LINE);
	//	printf("%s\n", source);
		while(fgets(buf, LINE, fp)) {
			  fscanf(fp, "%s", rr.name);
			  
			  if (!strcmp(source, rr.name)) {
				//  printf("yes\n"); 
				  source_exist = 1;
				  fscanf(fp, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
							 rr.p, rr.pdot, rr.dm, rr.ra, rr.dec, rr.l,
							 rr.b, rr.rate, rr.logB, rr.logts, rr.dhat, rr.fluxd, rr.pulse_width, rr.survey); 
				  P = atof(rr.p);
				  pw = atof(rr.pulse_width);
				  num = 5 * (pw/sample_dt/1000);
				 // printf("num: %d\n", num);
				  if ( 10 <= num && num <= 100) {
					  num = num * 10;
				  } else if ( 0 < num && num < 10){
					  num = num * 100;  
				  } else if (num == 0){
					   num = 300;
				  }
				  break;
			  } else if (!strcmp(source, "B0531+21")) {
				  source_exist = 1;
				  num = 200;
			  } 
		}
		
		//printf("num: %d\n", num);
		if (!source_exist) {
			printf("Source is not exist in the source file, Please change it!\n");
			exit(1);
		}
		fclose(fp);
		if (downsample_flag == 1) {
			downsample = (int)(pw/sample_dt/1000/8);
			if (downsample == 0) {
				printf("The downsample due to the pulse width and the sample_time is zero, Please set is by manual:\n");
				scanf("%d", &downsample);
			}
		}
		
		downsample_time = downsample * sample_dt;
		//printf("downsample: %d\n", downsample);
		/************************************************************/
		
		float *diff_dm_flux;
		float realdm_max = 0.0f, realdm_min = 0.0f;
		float otherdm_max = 0.0f;
		int realdm_max_samp = 0, realdm_min_samp = 0;
		int otherdm_max_samp = 0;
		int val_over_rms = 0;
		float rms;
		da->data_adm = (float *)malloc(sizeof(float) * da->numpoint);
		da->dm_x_coor = (float *)malloc(sizeof(float) * da->numpoint);
		diff_dm_flux = (float *)malloc(sizeof(float) * da->ndm);
		printf("the real dm:%f\n", da->real_dm);
		sprintf(opd->filename, "./DmFluxFile/%s/%f.dav", source, da->real_dm);
		opd->fp = open_file(opd->filename, "r");
		fread(da->dm_x_coor, sizeof(float), da->numpoint, opd->fp);
		fread(da->data_adm, sizeof(float), da->numpoint, opd->fp);
		realdm_max = *(da->data_adm + 0);
		fclose(opd->fp);
		*(da->data_adm + 1) = 0.0f;
		*(da->data_adm + da->numpoint - 1) = 0.0f;
		*(da->data_adm + da->numpoint - 2) = 0.0f;
		
		int num_after_downsample = da->numpoint/downsample + 1;
		int num1 = da->numpoint/downsample;
		int num2 = da->numpoint%downsample;
		//printf("%d\n", num_after_downsample);
		float *downsample_data;
		downsample_data = (float *)malloc(sizeof(float) * num_after_downsample);
		float downsample_temp = 0.0f;
		for (i = 0; i < num1; i++) {
			for(j = 0; j < downsample; j++) {
				downsample_temp += *(da->data_adm + j + i * downsample);
			}
			downsample_temp /= downsample;
			downsample_data[i] = downsample_temp;
			downsample_temp = 0;	
		}
		for (i = 0; i < num2; i++) {
			downsample_temp +=*(da->data_adm + i + num1 * downsample); 
		}
		downsample_temp /=downsample;
		downsample_data[num1] = downsample_temp;
		da->numpoint = num_after_downsample; 
		calc_rms(downsample_data, &rms, da->numpoint);
		
		for (i = 0; i < da->numpoint; i++) {
			if (*(downsample_data + i) > (set_thresh * rms)) {
				val_over_rms++;
			}
		}

		float dat_over_rms[val_over_rms];
		int point_over_rms[val_over_rms];
		int point_count = 0;
		j = 0;
		for (i = 0; i < da->numpoint; i++) {
			if (*(downsample_data + i) > (set_thresh * rms)) {
				dat_over_rms[j] = *(downsample_data + i);
				point_over_rms[j] = i;
				j++;
			}
		}
		//printf("j: %d\n", j);
		if (j == 0) {
			printf("Please reset the threshold or the downsample parameter, cannot find pulse\n");
			exit(1);
		}
		for (i = 0; i < j - 1; i++) {
			if ( (point_over_rms[i] + 1) == point_over_rms[i + 1]) {
				continue;
			} 
			point_count++;			
		}
		point_count++;
		float dat_over_rms1[point_count];
		int point_over_rms1[point_count];
		float temp = dat_over_rms[0];
		int samp_num = point_over_rms[0];
		k = 0;
		l = 0;
		for (i = 0; i < j - 1; i++) {
			if ((point_over_rms[i] + 1) == point_over_rms[i + 1]) {
				if (temp < dat_over_rms[i]) {
					temp = dat_over_rms[i];
					samp_num = point_over_rms[i];
				}
				l++;
				continue;
			} 
			dat_over_rms1[k] = l > 0? temp:dat_over_rms[i];
			point_over_rms1[k] = l > 0? samp_num:point_over_rms[i];
			temp = 0;
			l = 0;
			k++;
		}
		if ((point_over_rms[j - 2] + 1) == point_over_rms[j - 1]) {
			if (temp < dat_over_rms[j - 1]) {
				dat_over_rms1[k] = dat_over_rms[j - 1];
				point_over_rms1[k] = point_over_rms[j - 1];
			} else {
				dat_over_rms1[k] = temp;
				point_over_rms1[k] = samp_num;
			}
		} else {
			dat_over_rms1[k] = dat_over_rms[j - 1];
			point_over_rms1[k] = point_over_rms[j - 1];
		}
		//for(i = 0;i < val_over_rms; i++) printf("1: %d\n", point_over_rms[i]);
		//for(i = 0; i < point_count; i++) printf("2: %d\n", point_over_rms1[i]);
		printf("total point to plot:%d\n", point_count);
		/*for (i = 0; i < point_count; i++) { 
			point_over_rms1[i] = point_over_rms1[i] * downsample;
			printf("%f\n", dat_over_rms1[i]);
		}*/
		/******************/
		/* the bottom code to use the point from above to get data from the fits file */
		/* and use the data to plot flux and gray map*/
		/******************/		
		
		if (num % 100) {
			num = num - num%100;
		}
		
		int downsample_num = num * downsample;
		//printf("%d\n", downsample_num);
		int pw_num = (int)(pw/sample_dt/1000);
		
		//printf("%d %f\n", downsample, downsample_time);
		float plot_dat[num];
		float dat_x[num];
		int sum_point = point_count - 1;
	
		float x_min = 0.0f, x_max = 0.0f, y_min = 0.0f, y_max = 0.0f;
		
		psrfits_set_files(&pf, argc - optind, argv + optind);

		pf.tot_rows = pf.N = pf.T = pf.status = 0;
		pf.hdr.chan_dm = 0;
		
		/* open the fits file */
		int rv = psrfits_open(&pf);
		if (rv) { fits_report_error(stderr, rv); exit(1);}
		
		float *plot_chan = (float *)malloc(sizeof(float) * pf.hdr.nchan * num);
		int get_row_num;
		int get_row_point;
		int left_row = 0;
		int right_row = 0;
		int total_row = 0;
		int colnum = 0;
		float plot_max;
		float plot_min;
		/** the below is pgplot parameter **/
		float hi, lo, scale, x1, x2, xleft, xright, xscale;
		float y1, y2, ybottom, yscale, ytop;
		float xmin, xmax; 
		float ymin, ymax;
		float zmin, zmax;
		int nx, ny;
		float tr[6]; //tr[2] & tr[4] usually 0
		char output_ps[24];
		float *raw_chan;
		/***************/
		/** set the color index **/
		float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
		float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
		float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
		float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
		float contrast = 1.0;
		float brightness = 0.5;
		/*********/
		k = 0;
		system("rm ./*.ps");
		pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
		pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
		pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan
							* pf.hdr.npol);
		pf.sub.dat_scales = (float *)malloc(sizeof(float)
							*pf.hdr.nchan * pf.hdr.npol);
		pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	
		fits_get_colnum(pf.fptr, 0, "DAT_FREQ", &colnum, &pf.status);
		fits_read_col(pf.fptr, TFLOAT, colnum, 1, 1, pf.hdr.nchan, NULL, pf.sub.dat_freqs,
            NULL, &pf.status);
				
		while (k <= sum_point) {

			get_row_num = point_over_rms1[k] * downsample/pf.hdr.nsblk + 1;
			//printf("%d\n", point_over_rms1[k]);
			get_row_point = (point_over_rms1[k] * downsample)%pf.hdr.nsblk;
			if ((downsample_num/2 - get_row_point) > 0) {
				left_row = (downsample_num/2 - get_row_point)/pf.hdr.nsblk + 1;
			} else {
				left_row = 0;
			}
			
			if ((downsample_num/2 - (pf.hdr.nsblk - get_row_point)) > 0) {
				right_row = (downsample_num/2 - (pf.hdr.nsblk - get_row_point))/pf.hdr.nsblk + 1;
			} else {
				right_row = 0;
			}
			
			for (i = 0; i < pf.hdr.nchan * num; i++)
							*(plot_chan + i) == 0.0f;
			i = 0 - left_row;
			j = 0;
			total_row = left_row + right_row + 1;
			raw_chan = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol * pf.hdr.nsblk * total_row);
			if ((get_row_num + i) < 1) {
				k++;
				continue;
			}
			int count_0 = 0;
			int m = 0;
			
			while(i <= right_row) {
				fits_get_colnum(pf.fptr, 0, "DATA", &colnum, &pf.status);
			
				/*if (get_row_num + i < 1) {
					printf("Please reset the threshold, the row num of the fits file is over the first row!\n");
					exit(1);
				}*/
				fits_read_col(pf.fptr, TBYTE, colnum, get_row_num + i, 1, pf.sub.bytes_per_subint,
						  NULL, pf.sub.rawdata, NULL, &pf.status);
				if (pf.hdr.nbits==4) pf_4bit_to_8bit(pf);
				get_raw_chan(&pf, raw_chan, j);
				print_fits_error(pf.status);
				j++;
				i++;
			}
			//printf("left_row:%d, right_row:%d, sizeof_raw_chan:%d\n", left_row, right_row,  pf.hdr.nchan * pf.hdr.npol * pf.hdr.nsblk * total_row);
			//printf("total_row:%d, left_row:%d, get_row_point: %d, num: %d\n", total_row, left_row, get_row_point, num);
			
			//get_plot_chan(&pf, raw_chan, plot_chan, left_row, get_row_point, num);
			l = left_row * pf.hdr.nsblk + get_row_point;
			//printf("%d %d %d %d\n", (l - num/2) * pf.hdr.nchan, (l + num/2) * pf.hdr.nchan, l, pf.hdr.nsblk);
			//printf("%d\n", k);
			m = 0;
			float chan_temp = 0.0f;
			//printf("%d %d %d %s\n", pf.hdr.nchan, pf.hdr.npol, pf.hdr.nsblk, pf.hdr.poln_order);
			//if (k == 0) {
				/*for (i = 0; i < pf.hdr.nchan * pf.hdr.npol * pf.hdr.nsblk * total_row; i++) {
						//if (*(raw_chan + i) == 0.0)
						printf("raw_chan: %f\n", *(raw_chan + i));
					}
			//	}*/
			//for (i = 0; i < pf.hdr.nchan; i++) chan_temp[i] = 0.0f;
			//printf("num:%d nchan:%d\n", num, pf.hdr.nchan);
			//printf("%d %d\n", num, downsample);
			for (i = 0; i < num; i++) {
				
				for (j = 0; j < pf.hdr.nchan; j++) {
					/*if (k == 13) {
				
					for (i = 0; i < pf.hdr.nchan * pf.hdr.npol * pf.hdr.nsblk * total_row; i++)
						printf("raw_chan:%f\n", *(raw_chan + i));
					}*/
					for (m = 0; m < downsample; m++) chan_temp += *(raw_chan + j + (l - downsample_num/2 + m + i * downsample) * pf.hdr.nchan);
					chan_temp /= 4;
					*(plot_chan + i + j * num) = chan_temp;
					//printf("%f\n", *(plot_chan + i + j * num));
					chan_temp = 0.0f;
					/*if (k == 0) {
							printf("plot_chan:%f\n", *(plot_chan + i + j * num));
					}*/
				}
			}
			l = 0;
			for (i = 0; i < num; i++) {
				plot_dat[i] = *(downsample_data + point_over_rms1[k] - num/2  + i);
				dat_x[i] = (point_over_rms1[k] - num/2 + i) * sample_dt;
				downsample_temp = 0;
			}
			
			get_axis_mm(dat_x, plot_dat, &x_min, &x_max, &y_min, &y_max, num); 
			
			plot_max = *(plot_chan + 0);
			plot_min = *(plot_chan + 0);
			for (i = 0; i < pf.hdr.nchan * num; i++) {
				if (plot_max < *(plot_chan + i)) plot_max = *(plot_chan + i);
				if (plot_min > *(plot_chan + i)) plot_min = *(plot_chan + i);
			}
			sprintf(output_ps, "%d.ps/cps", k);
			if (cpgopen(output_ps) <= 0)
				return EXIT_FAILURE;
			cpgsubp(1, 2);
			cpgctab (heat_l, heat_r, heat_g, heat_b, 5, contrast, brightness);
			//cpgscf(1);
			//cpgsch(2);
			cpgenv(x_min, x_max, y_min, y_max, 0, 1);
			//cpgpage();
			//cpgvstd();
			//cpgswin(x_min, x_max, y_min, y_max);
			
			//cpgmtxt("B", 2.3, 0.5,1.0, "Time(s)");
			//cpgbox("1", 0.0, 0, "1", 0.0, 1);
			
			cpglab("Time(s)", "Flux(abitrary unit)", "Pulse Plot");
			cpgline(num, dat_x, plot_dat);
			
			
		
			nx = num;
			//ny = pf.sub.dat_freqs[pf.hdr.nchan];
			ny = pf.hdr.nchan;
			hi = plot_max;
			lo = plot_min;
			xmin = x_min;
			xmax = x_max;
			//ymin = pf.sub.dat_freqs[0];
			//ymin = 2256;
			//ymax = 2195;
			ymin = pf.sub.dat_freqs[0];
			ymax = pf.sub.dat_freqs[pf.hdr.nchan - 1];
			zmin = 0; zmax = 3;
			
			cpgenv(xmin, xmax, 1, ny + 5, 0, 1);
			cpglab("Time(s)", "Channel(MHz)","");
	
			tr[0] = x_min;
			tr[1] = (x_max - x_min)/nx;
			tr[2] = 0.0;
			tr[3] = ymin;
			tr[4] = 0.0;
			tr[5] = (ymax - ymin)/ny;
			float ylims[num];
			for (i = 0; i < num; i++) ylims[i] = i;
			//baseline_data(plot_chan, num * pf.hdr.nchan);
			cpghi2d(plot_chan, nx, ny, 1, nx, 1, ny, dat_x, 0, 1, 1, ylims); 
			//cpgimag(plot_chan, nx, ny, 1, nx, 1, ny, lo, hi, tr);
			//cpgaxis("1", xmin, ymin, xmax, ymin, x_min, x_max, 0.0001, 2, 0.3, 0.2, 0.1, 0.0, 0.0);
			cpgbbuf();
			cpgend();
			k++;
			free(raw_chan);
		}
		psrfits_close(&pf);
	}
}
