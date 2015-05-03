#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <getopt.h>
#include "psrfits.h"


void usage()
{
	printf(
		  "-h --help print this \n"
		  "-o --output output a file\n"
		  "-d --delchan delete channel (it is using \" \" like (\"1 2:13\"))\n"
		  "-r --readheader read header file from the fits file \n"
		  );
}

/* get the stokes I of data, dat_scales, dat_offsets  */
int get_AABB_I(struct psrfits *pf)
{
	struct subint *sub = &(pf->sub);
	struct hdrinfo *hdr = &(pf->hdr);
	unsigned char *data, *data1;
	float *dat_scales, *dat_scales1;
	float *dat_offsets, *dat_offsets1;
	int i = 0, j = 0;
	

	for(i = 0; i < hdr->nsblk; i++) {
		data = sub->rawdata + i * hdr->nchan * hdr->npol;
		data1 = data + hdr->nchan;
		for (j = 0; j < hdr->nchan; j++, data++, data1++) {
			*(sub->data + j + i * hdr->nchan) = *data + *data1;
		}
	}
	dat_scales = sub->dat_scales;
	dat_scales1 = dat_scales + hdr->nchan; 
	dat_offsets = sub->dat_offsets;
	dat_offsets1 = dat_offsets + hdr->nchan;
	for (j = 0; j < hdr->nchan; j++, dat_scales++, dat_scales1++, dat_offsets++, 
							dat_offsets1++) {
		*(sub->datI_scales + j) = 0.5 * ( *dat_scales + *dat_scales1);
		*(sub->datI_offsets + j) = 0.5 * (*dat_offsets + *dat_offsets1);	
	}
	return(0);
}

int get_only_I(struct psrfits *pf)
{
	struct subint *sub = &(pf->sub);
	struct hdrinfo *hdr = &(pf->hdr);
	unsigned char *data;
	float *dat_scales, *dat_offsets;
	int i = 0, j = 0;
	
	for (i = 0; i < hdr->nsblk; i++) {
		data = sub->rawdata + i * hdr->nchan * hdr->npol;
		for (j = 0; j < hdr->nchan; j++, data++) {
			*(sub->data + j + i * hdr->nchan) = *data;
		}
	}
	
	dat_scales = sub->dat_scales;
	dat_offsets = sub->dat_offsets;

	for (j = 0; j < hdr->nchan; j++, dat_scales++, dat_offsets++) {
		*(sub->datI_scales + j ) = *dat_scales;
		*(sub->datI_offsets + j ) = *dat_offsets;
	}
	
	return(0);
}
void del_chan(struct psrfits *pf_out, struct psrfits *pf, int *c[], int c_k, int sum_dchan, int *del)
{
	struct subint *sub = &(pf->sub);
	struct subint *sub1 = &(pf_out->sub);
	struct hdrinfo *hdr = &(pf->hdr);
	struct hdrinfo *hdr1 = &(pf_out->hdr);
	sub1->tsubint = sub->tsubint;
	sub1->offs = sub->offs;
	sub1->lst = sub->lst;
	sub1->ra = sub->ra;
	sub1->dec = sub->dec;
	sub1->glon = sub->glon;
	sub1->glat = sub->glat;
	sub1->feed_ang = sub->feed_ang;
	sub1->pos_ang = sub->pos_ang;
	sub1->par_ang = sub->par_ang;
	sub1->tel_az = sub->tel_az;
	sub1->tel_zen = sub->tel_zen;
	int i = 0, j = 0, k = 0, l = 0, m = 0;
	
	if (sub->dat_freqs[0] > sub->dat_freqs[1]) {
		for (i = 0; i < (hdr1->nchan); i++) {
			for (j = k; j < (hdr->nchan); j++) {
				if(del[j] == 1) {
					*(sub1->dat_freqs + i) = *(sub->dat_freqs + j);
					*(sub1->dat_weights + i) = *(sub->dat_weights + j);
					k = j + 1;
					break;
				}

			}
		}
		k = 0;
		for (i = 0; i < (hdr1->nchan); i++) {
			for(j = k; j < (hdr->nchan); j++) {
				if(del[j] == 1) {
					*(sub1->datI_offsets + i) = 
							*(sub->datI_offsets + j);
					*(sub1->datI_scales + i + l) = 
							*(sub->datI_scales + j);
					k = j + 1;
					break;
				}
			}
		}
		for (m = 0; m < hdr->nsblk; m++) {
			k = 0;
			for (i = 0; i < hdr1->nchan; i++) {
				for(j = k; j < hdr->nchan; j++) {
					if(del[j] == 1) {
						*(sub1->data + i + m * hdr1->nchan) = 
						 *(sub->data + j + m * hdr->nchan);
						 k = j + 1;
						 break;
					}
				}
			}
		}
	} else {
		for (i = hdr1->nchan - 1; i >= 0; i--) {
			for (j = k; j < (hdr->nchan); j++) {
				if(del[j] == 1) {
					*(sub1->dat_freqs + i) = *(sub->dat_freqs + j);
					*(sub1->dat_weights + i) = *(sub->dat_weights + j);
					k = j + 1;
					break;
				}

			}
		}
		k = 0;
		for (i = hdr1->nchan -1; i >= 0; i--) {
			for(j = k; j < (hdr->nchan); j++) {
				if(del[j] == 1) {
					*(sub1->datI_offsets + i) = 
							*(sub->datI_offsets + j);
					*(sub1->datI_scales + i) = 
						 *(sub->datI_scales + j);
					k = j + 1;
					break;
				}
			}
		}
		for (m = 0; m < hdr->nsblk; m++) {
			k = 0;
			for (i = hdr1->nchan -1 ; i >= 0; i--) {
				for(j = k; j < hdr->nchan; j++) {
					if(del[j] == 1) {
						*(sub1->data + i + m * hdr1->nchan) = 
						 *(sub->data + j + m * hdr->nchan);
						 k = j + 1;
						 break;
					}
				}
			}
		}
	}
}

int main(int argc, char *argv[])
{
	static struct option long_opts[] = {
		{"output", 1, NULL, 'o'},	
		{"delchan", 1, NULL, 'd'},
		{"help", 0, NULL, 'h'},
		{"readheader", 0, NULL, 'r'},
		{0, 0, 0, 0}
	};

	int opt, opti;
	int cal = 0;
	char output_base[256] = "\0";
	char del_str[256] = "";
	int readfile = 0;
	while ((opt = getopt_long(argc, argv, "o:d:hr", long_opts, &opti)) != -1)
	{
		switch(opt)
		{
			case 'o':
				strncpy(output_base, optarg, 255);
				output_base[255] = '\0';
				break;
			case 'd':
				strncpy(del_str, optarg, 255);
				del_str[255] = '\0';
				break;
			case 'r':
				readfile = 1;
				break;
			case 'h':
			default:
				usage();
				exit(0);
				break;
		}
	}
	if (optind == argc)
	{
		usage();
		exit(1);
	}

	/* Open file*/
	struct psrfits pf;
	psrfits_set_files(&pf, argc - optind, argv + optind);	
	pf.tot_rows = pf.N = pf.T = pf.status = 0;
	pf.hdr.chan_dm = 0.0;
	int rv = psrfits_open(&pf);
	if (rv) { fits_report_error(stderr, rv);  exit(1);}
	
	/**************************************************
	 * Print the header of the file and the header of *
	 * the subintegration.							  *
	 * ************************************************/
	 
	if (readfile) {
		fitsfile *fptr;
		char card[FLEN_CARD];
		int status = 0, nkeys, ii;
		fits_open_file(&fptr, *(argv + optind), READONLY, &status);
		fits_get_hdrspace(fptr, &nkeys, NULL, &status);
		printf("***************************************************\n");
		printf("FILE_HEADER:");
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
		printf("***************************************************\n");
		fits_close_file(fptr, &status);
		
		if (status) {
			fits_report_error(stderr, status);
		}
		exit(1);
	}
	
	char *del_str1, *saveptr, *saveptr1;
	int i = 0, j = 0, c_k= 0, sum_dchan = 0, k = 0;
	int *c[100];
	for (i = 0; i < 100; i++) { 
		c[i] = (int *)malloc(sizeof(int) * 2); 
		c[i][0] = c[i][1] = 0;
	}
	char *del_str2[50];
	char *buf = del_str;
	
	/* set the channel will be deleted */
	// use strtok_r to split the array
	while((del_str2[i] = strtok_r(buf, " ", &saveptr)) != NULL)
	{
		buf = del_str2[i];
		while((del_str2[i] = strtok_r(buf, ":", &saveptr1)) != NULL)
		{
			c[c_k][j] = atoi(del_str2[i]);
			i++;
			buf = NULL;
			j++;
		}
		j = 0;
		c_k++;
		buf = NULL;
	}

	for (j = 0; j < c_k; j++) {		
		sum_dchan += (((c[j][1] - c[j][0]) > 0)? (c[j][1] - c[j][0] + 1) : 1);
	}
		
	/* Set up output file*/	
	struct psrfits pf_out;
	struct hdrinfo *hdr = &(pf.hdr);
	memcpy(&pf_out, &pf, sizeof(struct psrfits));

	/* Set up default output filename */
	if (output_base[0] == '\0') {
		sprintf(output_base, "%s_%s_%5.5d_%5.5d%s", pf.hdr.backend,
				pf_out.hdr.source, pf_out.hdr.start_day,
				(int)pf_out.hdr.start_sec, cal ? "_cal" : "");
	}
	sprintf(pf_out.basefilename, output_base);
	pf_out.fptr = NULL;
	pf_out.filenum = 0;
	pf_out.status = 0;
	pf_out.quiet = 0;
	pf_out.hdr.onlyI = 1;
	pf_out.hdr.nchan = pf.hdr.nchan - sum_dchan;
	pf_out.hdr.npol = 1;
	pf_out.sub.FITS_typecode = TFLOAT;
	if (pf_out.hdr.df > 0) pf_out.hdr.df = 0 - pf_out.hdr.df; // set the channel bandwidth to negative to invert the band
	long long lltmp = pf_out.hdr.nsblk;
	lltmp = (lltmp * pf_out.hdr.nbits * pf_out.hdr.nchan * pf_out.hdr.npol) / 8L;
	pf_out.sub.bytes_per_subint = (int) lltmp;
	
	//rv = psrfits_create(&pf_out);
	if (rv) { fits_report_error(stderr, rv); exit(1);}
	
	/* set one array to store the channel will be deleted */
	int *del = (int *)malloc((hdr->nchan)*sizeof(int));
	for (i = 0; i < (hdr->nchan); i++) { del[i] = 1;}
	for(j = 0; j < c_k; j++) {
		if (c[j][0] > c[j][1]) {
		   del[c[j][0] - 1] = 0;
		} else {
			for(i = c[j][0] - 1; i <=c[j][1] - 1; i++) {
				del[i] = 0;
			}
		}
	}
	/* Alloc data buffers */

	pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf_out.sub.dat_freqs = (float *)malloc(sizeof(float) * pf_out.hdr.nchan);
	pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf_out.sub.dat_weights = (float *)malloc(sizeof(float) * pf_out.hdr.nchan);
	pf.sub.dat_offsets = (float *)malloc(sizeof(float)
			* pf.hdr.nchan * pf.hdr.npol);
	pf.sub.datI_offsets = (float *)malloc(sizeof(float)
			* pf.hdr.nchan);
	pf_out.sub.dat_offsets = (float *)malloc(sizeof(float)
			* pf_out.hdr.nchan * pf_out.hdr.npol);
	pf_out.sub.datI_offsets = (float *)malloc(sizeof(float)
			* pf_out.hdr.nchan);
	pf.sub.dat_scales = (float *)malloc(sizeof(float)
			* pf.hdr.nchan * pf.hdr.npol);
	pf.sub.datI_scales = (float *)malloc(sizeof(float)
			* pf.hdr.nchan);
	pf_out.sub.dat_scales = (float *)malloc(sizeof(float)
			* pf_out.hdr.nchan * pf_out.hdr.npol);
	pf_out.sub.datI_scales = (float *)malloc(sizeof(float)
			* pf_out.hdr.nchan);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	pf_out.sub.rawdata = (unsigned char *)malloc(pf_out.sub.bytes_per_subint);
	pf.sub.data = (unsigned char *)malloc(pf.hdr.nsblk * pf.hdr.nchan);
	pf_out.sub.data = (unsigned char *)malloc(pf_out.hdr.nsblk * pf_out.hdr.nchan);
	int run = 1;
	rv = 0;
	int first_cal_cf = 1; // set to caculate the center frequency
	while(run) {
		rv = psrfits_read_subint(&pf);
		if (rv) {
			if (rv == FILE_NOT_OPENED) rv = 0;
			run = 0;
			break;
		}

		if ((strncmp(pf.hdr.poln_order, "AABBCRCI", 8) == 0) || 
			 (strncmp(pf.hdr.poln_order, "AABB", 4) == 0)) {
			get_AABB_I(&pf);

		} else if ((strncmp(pf.hdr.poln_order, "IQUV", 4) == 0) ||
				(strncmp(pf.hdr.poln_order, "AA", 2) == 0) ||
				(strncmp(pf.hdr.poln_order, "AA+BB") == 0)) {
			get_only_I(&pf);
		}	
		
		del_chan(&pf_out, &pf, c, c_k, sum_dchan, del);	

		if (first_cal_cf) {
			if (pf_out.hdr.nchan == 1) {
				pf_out.hdr.fctr = (double)(pf_out.sub.dat_freqs[pf_out.hdr.nchan - 1]);
			}else if ((pf_out.hdr.nchan % 2) == 0) {
				pf_out.hdr.fctr = (double)(0.5 * (pf_out.sub.dat_freqs[pf_out.hdr.nchan/2] +
										pf_out.sub.dat_freqs[pf_out.hdr.nchan/2 -1 ]));
			} else {
				pf_out.hdr.fctr = (double)(pf_out.sub.dat_freqs[(pf_out.hdr.nchan - 1)/2]);
			}
			first_cal_cf = 0;
		}

		psrfits_write_subint(&pf_out);
		printf("\rWrite subint(file %d, row %d/%d)",
				pf_out.filenum, pf_out.rownum - 1, pf_out.rows_per_file);
		fflush(stdout);
	}
	psrfits_close(&pf);
	psrfits_close(&pf_out);
	
	if (rv > 100) {fits_report_error(stderr, rv);}
	exit(0);
}
