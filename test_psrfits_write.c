#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <getopt.h>
#include <signal.h>
#include <pthread.h>

#include "psrfits.h"

struct dc_args {
	int **c;
	int c_k;
	struct psrfits *pf_out;
	struct psrfits *pf;
};

void usage()
{
	printf(
		  "-h, --help print this \n"
		  "-o name, --output=name output base filename (auto-generate)\n"
		  "-d --delchan delete channel (it is using \" \" like (\"1 2:13\"))\n"
		  "-j --nthread=nn Max number of threads (4)\n"
		  );
}
/* Singnal handler*/
int run = 1;
void cc(int sig) { run = 0; }

int copy_data(struct subint *sub1, struct subint *sub)
{
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
	return(0);
}
int del_chan(struct psrfits *pf_out, struct psrfits *pf, int **c, int c_k)
{
	int i = 0, j = 0, k = 0, l = 0, m = 0;
	struct hdrinfo *hdr = &(pf->hdr);
	struct hdrinfo *hdr1 = &(pf_out->hdr);
	struct subint *sub = &(pf->sub);
	struct subint *sub1 = &(pf_out->sub);
	copy_data(sub1, sub);
	int del[hdr->nchan];
	//printf("2\n");
	for (i = 0; i < hdr->nchan; i++) { del[i] = 1;}
	//printf("3\n");
	for(j = 0; j < c_k; j++) {
		if (c[j][0] > c[j][1]) {
		   del[c[j][0] - 1] = 0;
		} else{
			for(i = c[j][0] - 1; i <=c[j][1] - 1; i++) {
				del[i] = 0;
			}
		}
	}
	//printf("4\n");
	for (i = 0; i < hdr1->nchan; i++) {
		for (j = k; j < hdr->nchan; j++) {
			if(del[j] == 1) {
				sub1->dat_freqs[i] = sub->dat_freqs[j];
				sub1->dat_weights[i] = sub->dat_weights[j];
				k = j + 1;
				break;
			}

		}
	}
	k = 0;
	for (l = 0; l < hdr->npol; l++) {
		k = 0;
		for (i = 0; i < hdr1->nchan; i++) {
			for(j = k; j < hdr->nchan; j++) {
				if(del[j] == 1) {
						sub1->dat_offsets[i +  l * hdr1->nchan] = 
								sub->dat_offsets[j + l * hdr->nchan];
						sub1->dat_scales[i + l * hdr1->nchan] = 
								sub->dat_scales[j + l * hdr->nchan];
						k = j + 1;
						break;
					}
				}
			}
		}
	for (m = 0; m < hdr->nsblk; m++) {
		for (l = 0; l < hdr->npol; l++) {
			k = 0;
			for (i = 0; i < hdr1->nchan; i++) {
				for(j = k; j < hdr->nchan; j++) {
					if(del[j] == 1) {
						sub1->rawdata[i + l * hdr1->nchan + 
							m * hdr1->nchan * hdr1->npol] = 
						 sub->rawdata[j + l * hdr->nchan + 
							m * hdr->nchan * hdr->npol];
						 k = j + 1;
						 break;
					}
				}
				
			}
		}
	}
	
	return(0);
}

void *del_chan_thread(void *_args) {
	struct dc_args *args = (struct dc_args *)_args;
	int rv = psrfits_read_subint(args->pf);
		printf("Read subint(file %d, row %d/%d)\n",
				args->pf->filenum, args->pf->rownum - 1, args->pf->rows_per_file);
	if (rv) {
			if (rv == FILE_NOT_OPENED) rv = 0;
			run = 0;
		}
	rv = del_chan(args->pf_out, args->pf, args->c, args->c_k);
	psrfits_write_subint(args->pf_out);
		printf("Write subint(file %d, row %d/%d)\n",
				args->pf_out->filenum, args->pf_out->rownum - 1, args->pf_out->rows_per_file);
	if (rv) {
		fits_report_error(stderr, rv);
		exit(1);
	}
	pthread_exit(&rv);
}
int main(int argc, char *argv[])
{
	static struct option long_opts[] = {
		{"output", 1, NULL, 'o'},	
		{"delchan", 1, NULL, 'd'},
		{"help", 0, NULL, 'h'},
		{"ntread", 1, NULL, 'j'},
		{0, 0, 0, 0}
	};

	int opt, opti;
	int cal = 0;
	char output_base[256] = " ";
	char del_str[256] = " ";
	int nthread = 4;

	while ((opt = getopt_long(argc, argv, "o:d:f", long_opts, &opti)) != -1)
	{
		switch(opt)
		{
			case 'o':
				strncpy(output_base, optarg, 255);
				printf("%s\n", output_base);
				output_base[255] = '\0';
				break;
			case 'd':
				strncpy(del_str, optarg, 255);
				del_str[256] = '\0';
				break;
			case 'j':
				nthread = atoi(optarg);
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
	printf("optind:%d,argc - optind:%d, argv + optind:%s\n", optind, argc-optind, *(optind + argv));
	//scanf("%d%d", del[0], del[1]);
	/* Open file*/
	struct psrfits pf;
	psrfits_set_files(&pf, argc - optind, argv + optind);	
	pf.tot_rows = pf.N = pf.T = pf.status = 0;
	pf.hdr.chan_dm = 0.0;
	int rv = psrfits_open(&pf);
	if (rv) { fits_report_error(stderr, rv);  exit(1);}
	

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
	//printf("%s\n", output_base);
	sprintf(pf_out.basefilename, output_base);
	pf_out.fptr = NULL;
	pf_out.filenum = 0;
	pf_out.status = 0;
	pf_out.quiet = 0;
	pf_out.hdr.nchan = pf.hdr.nchan - sum_dchan;
	pf_out.sub.FITS_typecode = TFLOAT;
	long long lltmp = pf_out.hdr.nsblk;
	lltmp = (lltmp * pf_out.hdr.nbits * pf_out.hdr.nchan * pf_out.hdr.npol) / 8L;
	pf_out.sub.bytes_per_subint = (int) lltmp;
	
	rv = psrfits_create(&pf_out);
	if (rv) { fits_report_error(stderr, rv); exit(1);}
	
	/* Alloc data buffers */

	pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf_out.sub.dat_freqs = (float *)malloc(sizeof(float) * pf_out.hdr.nchan);
	pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf_out.sub.dat_weights = (float *)malloc(sizeof(float) * pf_out.hdr.nchan);
	pf.sub.dat_offsets = (float *)malloc(sizeof(float)
			* pf.hdr.nchan * pf.hdr.npol);
	pf_out.sub.dat_offsets = (float *)malloc(sizeof(float)
			* pf_out.hdr.nchan * pf.hdr.npol);
	pf.sub.dat_scales = (float *)malloc(sizeof(float)
			* pf.hdr.nchan * pf.hdr.npol);
	pf_out.sub.dat_scales = (float *)malloc(sizeof(float)
			* pf_out.hdr.nchan* pf.hdr.npol);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	pf_out.sub.rawdata = (unsigned char *)malloc(pf_out.sub.bytes_per_subint);
	//printf("%d\n", pf.sub.bytes_per_subint);
	//printf("%d\n", pf_out.sub.bytes_per_subint);
	printf("1\n");
	pthread_t *thread_id;
	struct dc_args *dargs;
	thread_id = (pthread_t *)malloc(sizeof(pthread_t) * nthread);
	dargs = (struct dc_args *)malloc(sizeof(struct dc_args) * nthread);
	printf("2\n");
	for (i = 0; i < nthread; i++) {
		thread_id[i] = 0;
		
		memcpy(dargs[i].pf, &pf, sizeof(struct psrfits));
		
		memcpy(dargs[i].pf_out, &pf_out, sizeof(struct psrfits));
		dargs[i].c = c;
		dargs[i].c_k = c_k;		
	}
	printf("3\n");
	rv = 0;
	signal(SIGINT, cc);
	int cur_thread = 0;
	while(run ) {
		rv = pthread_create(&thread_id[cur_thread], NULL,
							del_chan_thread, &dargs[cur_thread]);
		if(rv) {
			fprintf(stderr, "Thread creation error.\n");
			exit(1);
		}
		cur_thread++;
		if (cur_thread == nthread) {
			for(i = 0; i < cur_thread; i++) {
				rv = pthread_join(thread_id[i], NULL);
				if (rv) {
					fprintf(stderr, "Thread join error.\n");
					exit(1);
				}
				thread_id[i] = 0;
			}
			
			cur_thread = 0;
		}
		
		/* If we'v passed final file,exit */
		//if(fnum_)
		//printf("\n");
		
		//if((pf.rownum - 1) == pf.rows_per_file) { break; }
	}
	psrfits_close(&pf);

}
