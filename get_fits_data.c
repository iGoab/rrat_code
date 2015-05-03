#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "psrfits.h"
#include "./pulse_dm.h"

int main(int argc, char *argv[])
{
	static struct option long_opts[] = {
		{"read", 0, NULL, 'r'},
		{"help", 0, NULL, 'h'},
		{0, 0, 0, 0}
	};
	
	int opt, opti, dm_flag = 0, status = 0;
	
	int i = 0, j = 0, k = 0;
	char line[256], line1[256];
	int read_file = 0;
	while((opt = getopt_long(argc, argv, "rh", long_opts, &opti)) != -1) {
		switch(opt) {
				case 'r':
					read_file = 1;
					break;
				case 'h':
				default:
					get_fits_data_help();
					exit(0);
					break;
			}
	}
	if (optind == argc) {
		get_fits_data_help();
		exit(1);
	}
	
	/**************************************************
	 * print the header of the file and the header of *
	 * the subintegration header					  *
	 * ************************************************/
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
	pf.hdr.chan_dm  = 0;
	
	/* open the fits file */
	int rv = psrfits_open(&pf);
	if (rv) {fits_report_error(stderr, rv); exit(1);}

	int nr = pf.hdr.nsblk * pf.rows_per_file;
	float *data;
	/* malloc data buf */
	pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
	pf.sub.dat_scales = (float *)malloc(sizeof(float)
								*pf.hdr.nchan * pf.hdr.npol);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	data = (float *)malloc(sizeof(float) * pf.hdr.nsblk * pf.hdr.nchan);
	int run = 1, row_count = 0, chan_count = 0;
	
	struct output_dmdata *opd;
	
	opd = (struct output_dmdata *)malloc(sizeof(struct output_dmdata));
	struct stat buf;
	printf("Get psrfits raw data and output as *.ch \n");
	sprintf(line, "rm ./RawData/%s/*.ch",
            pf.hdr.source);

    system(line);
	sprintf(line1, "./RawData/%s", pf.hdr.source);
	sprintf(line, "mkdir ./RawData/%s", pf.hdr.source);
	//system("rm ./RawData/*.ch");
	if (stat(line1, &buf) == -1) system(line);
	
	
	if ((strncmp(pf.hdr.poln_order, "AABBCRCI", 8) == 0) ||
			(strncmp(pf.hdr.poln_order, "AABB", 4) == 0)) {
		while(run) {
			rv = psrfits_read_subint(&pf);
			if (rv) {
				if (rv == FILE_NOT_OPENED) rv = 0;
				run = 0;
				break;
			}
			get_AABB_I(&pf, data);
			
			for (j = 0; j < pf.hdr.nsblk; j++)
				for(i = 0; i < pf.hdr.nchan; i++) 
					fprintf(opd[i].fp, "%d\t%d\t%f\n", pf.rownum - 1, j, *(data + i + j * pf.hdr.nchan));
			
			printf("\rProcess subint(file %d row %d/%d)",
												pf.filenum, pf.rownum - 1, pf.rows_per_file);
		}		
	} else if ((strncmp(pf.hdr.poln_order, "IQUV", 4) == 0) ||
						(strncmp(pf.hdr.poln_order, "AA", 2) == 0) ||
						(strncmp(pf.hdr.poln_order, "AA+BB", 5) == 0)) {
		while(run) {
			rv = psrfits_read_subint(&pf);
			if (rv) {
				if (rv == FILE_NOT_OPENED) rv = 0;
				run = 0;
				break;
			}
			
			get_only_I(&pf, data);
	
			for(i = 0; i < pf.hdr.nsblk; i++) {
				for (j = 0; j < pf.hdr.nchan; j++) {
					sprintf(opd->filename, "./RawData/%s/%d.ch", pf.hdr.source, j + 1);
					opd->fp = open_file(opd->filename, "a+");
					fprintf(opd->fp, "%f\n", *(data + j + i * pf.hdr.nchan));
					fclose(opd->fp);
				}
			}
			printf("\rProcess subint(file %d row %d/%d)",
												pf.filenum, pf.rownum - 1, pf.rows_per_file);
		}
	}
	for (i = 0; i < pf.hdr.nchan; i++) fclose(opd[i].fp);
	free(data);
	psrfits_close(&pf);
	if (rv > 100) {fits_report_error(stderr, rv); }
	exit(0);
}
