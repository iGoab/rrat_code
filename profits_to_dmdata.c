#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include "./pulse_dm.h"
#include "psrfits.h"

#define PI M_PI


int main(int argc, char *argv[])
{
	static struct option long_opts[] = {
		{"output", 1, NULL, 'o'},
		{"dedispersion", 1, NULL, 'd'},
		{"dmmin", 1, NULL, 'i'},
		{"nthread", 1, NULL, "n"},
		{"dmmax", 1, NULL, 'a'},
		{"dmstep", 1, NULL, 's'},
		{"read", 0, NULL, 'r'},
		{"zerodm", 0, NULL, 'z'},
		{"help", 0, NULL, 'h'},
		{0, 0, 0, 0}
	};
	
	int opt, opti, dm_flag = 0, status = 0;
	double dm_min = 0.0f, dm_max = 1000.0f, dm_step = 5.0f; 
	int ndm = 0; // set how many dm will dedisperse data
	// set dm_flag equals 1 if use the real dm, or set the dm_flag equals 0
	int i = 0, j = 0, k = 0;
	float temp = 0.0f;
	float *data, *data1;
	char line[256], line1[256];
	double dm = 0.0;
	char output_base[256] = "\0";
	int read_file = 0;
	int nthread = 4;
	int apply_zerodm = 0;
	while((opt = getopt_long(argc, argv, "o:d:i:a:s:n:rzh", long_opts, &opti)) != -1) {
		switch(opt) {
			case 'o':
				strncpy(output_base, optarg, 255);
				output_base[255] = '\0';
				break;
			case 'd':
				dm = atof(optarg);
				dm_flag = 1;
				break;
			case 'i':
				dm_min = atof(optarg);
				break;
			case 'a':
				dm_max = atof(optarg);
				break;
			case 'n':
				nthread = atoi(optarg);
				break;
			case 's':
				dm_step = atof(optarg);
				break;
			case 'r':
				read_file = 1;
				break;
			case 'z':
				apply_zerodm = 1;
				break;
			case 'h':
			default:
				profits_to_dmdata_help();
				exit(0);
				break;
		}
	}
	if (optind == argc) {
		profits_to_dmdata_help(0);
		exit(1);
	}
	
	/******************************************************
	 * print the header of the file and the header of     *
	 * the subintegration header 						  *
	 * ****************************************************/
	if (read_file) {
		fitsfile *fptr;
		fits_open_file(&fptr, *(argv + optind), READONLY, &status);
		status = print_fits_header(fptr);
		print_fits_error(status);
		exit(0);
	}
	
	struct psrfits pf;
	psrfits_set_files(&pf, argc - optind, argv + optind);
	pf.tot_rows = pf.N = pf.T = pf.status = 0;
	pf.hdr.chan_dm = 0;
	
	/* open the fits file */
	int rv = psrfits_open(&pf);
	if (rv) {fits_report_error(stderr, rv); exit(1);}
	
	int nr = pf.hdr.nsblk * pf.rows_per_file;

	/* malloc data buf */
	pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan
						* pf.hdr.npol);
	pf.sub.dat_scales = (float *)malloc(sizeof(float)
						*pf.hdr.nchan * pf.hdr.npol);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	data = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf.hdr.nchan);
	data1 = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf.hdr.nchan);

	if (dm_flag)
		ndm = (int)((dm_max - dm_min)/dm_step) + 1;
	else 
		ndm = (int)((dm_max - dm_min)/dm_step);
	//printf("%d\n", ndm);
	/* Malloc data_alldm to store the data after dedispersion*/
	struct data_alldm da;
	da.ndm = ndm;
	da.dm = (double *)malloc(sizeof(double) * ndm);
	da.numpoint = 0;
	da.avg_flux = 0.0f;
	da.dm_data_max = 0.0f;
	da.dm_data_min = 0.0f;
	da.real_dm = dm;
	da.idelays = (int **)malloc(sizeof(int *) * ndm);
	for (i = 0; i < ndm; i++) {
		da.idelays[i] = (int *)malloc(sizeof(int *) * pf.hdr.nchan);
	}
	da.data_adm = (float *)malloc(sizeof(float) * pf.hdr.nsblk);
	for (i = 0; i < ndm; i++) {
		if (dm_flag) {
			if (i == ndm - 1)
				*(da.dm + i) = dm;
			else 
				*(da.dm + i) = dm_min + i * dm_step;
		} else {
			*(da.dm + i) = dm_min + i * dm_step;
		}
	}
	
	if (dm_flag) {
		/* use straight insertion to sort dm*/
		for (i = 1; i < ndm; i++) {
			j = i;
			temp = *(da.dm + i);
			while(j > 0 && da.dm[j - 1] > temp) {
				da.dm[j] = da.dm[j - 1];
				j--;
			}
			*(da.dm + j) = temp;
		}
	}

	int run = 1, row_count = 0, chan_count = 0;
	
	int calc_data_dm = 1;
	int flag_check_freq = 0;
	float sum = 0.0f, avg = 0.0f;
	struct output_dmdata *opd;
	opd = (struct output_dmdata *)malloc(sizeof(struct output_dmdata));
	struct stat buf;
	printf("Process data now and using dm to dedisperse\n");

	sprintf(line, "rm -rf DmFluxFile/%s",
            pf.hdr.source);

    system(line);
    sprintf(line, "mkdir DmFluxFile");
    system(line);
    sprintf(line, "mkdir ./DmFluxFile/%s", pf.hdr.source);
    system(line);
	sprintf(line1, "./DmFluxFile/%s", pf.hdr.source);
	sprintf(line, "mkdir ./DmFluxFile/%s", pf.hdr.source);
	//system("rm ./RawData/*.ch");
	if (stat(line1, &buf) == -1) {system(line);
	}
	/* output the dm information to get ready for the following process */
	sprintf(opd->filename, "./DmFluxFile/%sdm.header", pf.hdr.source);
	opd->fp = open_file(opd->filename, "wb");
	fwrite(pf.hdr.source, sizeof(char), 24, opd->fp);
	fwrite(&(pf.hdr.dt), sizeof(double), 1, opd->fp);
	fwrite(&nr, sizeof(int), 1, opd->fp);
	fwrite(&ndm, sizeof(int), 1, opd->fp);
	fwrite(&dm_flag, sizeof(int), 1, opd->fp);
	float dm_temp = 0.0f;
	dm_temp = dm - *(da.dm + 0);
	if (dm_temp < 0) dm_temp = 0 - dm_temp;
	//printf("%f\n", dm_temp);
	for (i = 1; i < ndm; i++) {
		dm_temp = dm - *(da.dm + i);
		
		if (dm_temp < 0) dm_temp = 0 - dm_temp;
		//printf("%f\n", dm_temp);
		//printf("%f\n", abs(dm - *(da.dm + i)));
		if (dm_temp< 1e-3) {
			k = i;
		}
	}
	//printf("%d\n", k);
	if (dm_flag == 1) {
		fwrite((da.dm + k), sizeof(double), 1, opd->fp);
	}
	//printf("%f\n", *(da.dm + k));
	//fprintf(opd->fp, "%s\n", pf.hdr.source);
	//fprintf(opd->fp, "%d\n", nr);
	//fprintf(opd->fp, "%d\n", ndm);

	fwrite(da.dm, sizeof(double), ndm, opd->fp);
		//fprintf(opd->fp, "%f\n", da.dm[i]);
	fclose(opd->fp);
	int getdata = 1;
	int getdata1 = 0;
	int cur_thread = 0;
	pthread_t *thread_id;
	struct output_buf *opb = (struct output_buf *)malloc(sizeof(struct output_buf) * nthread);
	for (i = 0; i < nthread; i++) {
		opb[i].odata = (float *)malloc(sizeof(float) * pf.hdr.nsblk);
		opb[i].size = pf.hdr.nsblk;
	}
	
	if ((strncmp (pf.hdr.poln_order, "AABBCRCI", 8) == 0) ||
					(strncmp(pf.hdr.poln_order, "AABB", 4) == 0)) {
		while(run) {
			/* read two row then process the data */
					
					rv = psrfits_read_subint(&pf);
					if (rv) {
						if (rv == FILE_NOT_OPENED) rv = 0;
						run = 0;
						break;
					}
					
					if (getdata == 1) {
						getdata = 0;
						get_AABB_I(&pf, data);
					} else {
						getdata = 1;
						get_AABB_I(&pf, data1);
					}
					/* calculate the delays, just one time */
					if (calc_data_dm) {
						calc_data_dm = 0;
						for (i = 0; i < ndm; i++) {
							flag_check_freq = calc_dm (&pf, da.idelays, da.dm, ndm);
						}
					}
					/*for (j = 0; j < ndm; j++)
					for (i = 0; i < pf.hdr.nchan; i++) {
						printf("%d: %d\n", j, da.idelays[j][i]);
					}*/
					/* if want to clip the rfi, use the zerodm algorithm*/	
					if (apply_zerodm == 1) {
						if (getdata == 0) {
							zerodm(&pf, data);
						} else {
							zerodm(&pf, data1);
						}
					}
					row_count++;
					/* use dm to process the data and output the data to .dm file*/
					if (row_count >= 2) {
						
						for (i = 0; i < ndm; i++) {
							if (row_count & 1) {
								delay_data(&pf, data1, data, &da, flag_check_freq, i);
							} else {

								delay_data(&pf, data, data1, &da, flag_check_freq, i);
							}
							
							sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opd->fp = open_file(opd->filename, "a+");						
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							fclose(opd->fp);
						}
					}
					
					if (row_count == pf.rows_per_file) {

						for (i = 0; i < pf.hdr.nsblk; i++) {
							*(data + i) = 0.0f;
						}
						for (i = 0; i < ndm; i++) {
							delay_data(&pf, data1, data, &da, flag_check_freq, i);
							sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opd->fp = open_file(opd->filename, "a+");						
								//fprintf(opd->fp, "%f\n", *(da.data_adm + j));
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							fclose(opd->fp);
						}
					}
					printf("\rProcess subint(file %d row %d/%d)",
											pf.filenum, pf.rownum - 1, pf.rows_per_file);
					fflush(stdout);
		}				
	} else if ((strncmp(pf.hdr.poln_order, "IQUV", 4) == 0) ||
					(strncmp(pf.hdr.poln_order, "AA", 2) == 0) ||
					(strncmp(pf.hdr.poln_order, "AA+BB", 5) == 0)) {
				while(run) {
					
					/* read two row then process the data */
					
					rv = psrfits_read_subint(&pf);
					if (rv) {
						if (rv == FILE_NOT_OPENED) rv = 0;
						run = 0;
						break;
					}
					
					if (getdata == 1) {
						getdata = 0;
						get_only_I(&pf, data);
					} else {
						getdata = 1;
						get_only_I(&pf, data1);
					}
					/* calculate the delays, just one time */
					if (calc_data_dm) {
						calc_data_dm = 0;
						for (i = 0; i < ndm; i++) {
							flag_check_freq = calc_dm (&pf, da.idelays, da.dm, ndm);
						}
					} 
					/*for (j = 0; j < ndm; j++)
					for (i = 0; i < pf.hdr.nchan; i++) {
						printf("%d: %d\n", j, da.idelays[j][i]);
					}*/
					/* if want to clip the rfi, use the zerodm algorithm*/	
					if (apply_zerodm == 1) {
						if (getdata == 0) {
							zerodm(&pf, data);
						} else {
							zerodm(&pf, data1);
						}
					}
					row_count++;
					
					/*if (row_count >= 2) {
						
						for (i = 0; i < ndm; i++) {
							if (row_count & 1) {
								delay_data(&pf, data1, data, &da, flag_check_freq, i);
							} else {

								delay_data(&pf, data, data1, &da, flag_check_freq, i);
							}

							sprintf(opb[cur_thread].filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opb[cur_thread].fp = open_file(opb[cur_thread].filename, "a+");	
								printf("%s\n", opb[cur_thread].filename);
							memcpy(opb[cur_thread].odata, da.data_adm, sizeof(float) * pf.hdr.nsblk);
							rv = pthread_create(&thread_id[cur_thread], NULL, (void *)thread_function, &opb[cur_thread]);
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							if (rv) fprintf(stderr, "Create thread error.\n");
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							cur_thread++;
							if (cur_thread == nthread || i + 1 == ndm) {
								for (j = 0; j < cur_thread; j++) {
									rv = pthread_join(thread_id[j], NULL);
									if (rv) fprintf(stderr, "Join thread error.\n");
									thread_id[j] = 0;
									memset(opb[j].odata, 0, sizeof(float) * pf.hdr.nsblk);
								} 
								cur_thread = 0;
							}
						}
					}
					
					if (row_count == pf.rows_per_file) {

						for (i = 0; i < pf.hdr.nsblk; i++) {
							*(data + i) = 0.0f;
						}
						for (i = 0; i < ndm; i++) {
							delay_data(&pf, data1, data, &da, flag_check_freq, i);
							/*sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opd->fp = open_file(opd->filename, "a+");						
								//fprintf(opd->fp, "%f\n", *(da.data_adm + j));
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							fclose(opd->fp);
							sprintf(opb[cur_thread].filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opb[cur_thread].fp = open_file(opb[cur_thread].filename, "a+");	
							printf("%s\n", opb[cur_thread].filename);
							memcpy(opb[cur_thread].odata, da.data_adm, sizeof(float) * pf.hdr.nsblk);
							rv = pthread_create(&thread_id[cur_thread], NULL, (void *)thread_function, &opb[cur_thread]);
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							if (rv) fprintf(stderr, "Create thread error.\n");
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							cur_thread++;
							if (cur_thread == nthread || i + 1 == ndm) {
								for (j = 0; j < cur_thread; j++) {
									rv = pthread_join(thread_id[j], NULL);
									if (rv) fprintf(stderr, "Join thread error.\n");
									thread_id[j] = 0;
									memset(opb[j].odata, 0, sizeof(float) * pf.hdr.nsblk);
								} 
								cur_thread = 0;
							}
						}
					}*/
					/* use dm to process the data and output the data to .dm file*/
					if (row_count >= 2) {
						for (i = 0; i < ndm; i++) {
							if (row_count & 1) {
								delay_data(&pf, data1, data, &da, flag_check_freq, i);
							} else {
								delay_data(&pf, data, data1, &da, flag_check_freq, i);
							}
							sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opd->fp = open_file(opd->filename, "a+");						
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							fclose(opd->fp);
						}
					}
					
					if (row_count == pf.rows_per_file) {

						for (i = 0; i < pf.hdr.nsblk; i++) {
							*(data + i) = 0.0f;
						}
						for (i = 0; i < ndm; i++) {
							delay_data(&pf, data1, data, &da, flag_check_freq, i);
							sprintf(opd->filename, "./DmFluxFile/%s/%f.dm", pf.hdr.source, da.dm[i]);
							opd->fp = open_file(opd->filename, "a+");						
								//fprintf(opd->fp, "%f\n", *(da.data_adm + j));
							fwrite(da.data_adm, sizeof(float), pf.hdr.nsblk, opd->fp);
							
							for (j = 0; j < pf.hdr.nsblk; j++) {
								*(da.data_adm + j) = 0.0f;
							}
							fclose(opd->fp);
						}
					}
					printf("\rProcess subint(file %d row %d/%d)",
											pf.filenum, pf.rownum - 1, pf.rows_per_file);
					fflush(stdout);
				}		
		} 
		free(data);
		free(data1);
		free(da.data_adm);
}

