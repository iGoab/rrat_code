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

double dmdelay(float f1, float f2, double dm) /* includefile */
{
  return(4148.741601*((1.0/f1/f1)-(1.0/f2/f2))*dm);
}

int *dmshift(float f1, double df, int nchans, int nbands, double dm, double refrf, double tsamp, float frequency[])
{
	int i, cpb, *shift;
	double a;
	float fi;
	shift = (int *)malloc(nchans * sizeof(int));
	fi = f1;
	cpb = nchans/nbands;
	if (frequency[0] != 0.0) f1 = frequency[0];
	for (i = 0; i < nchans; i++) {
		if (refrf > 0.0) f1 = refrf;
		
		if (frequency[0] != 0.0) fi = frequency[i];
		shift[i]=(int)(dmdelay(f1,fi,dm)/tsamp);
		fi += df;
		if (!((i + 1) % cpb)) f1 += (double)cpb * df;
	}
	return shift;
}


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
	dm_temp = abs(dm - *(da.dm + 0));
	for (i = 1; i < ndm; i++) {
		if (dm_temp > abs(dm - *(da.dm + i))) {
			k = i;
			dm_temp = abs(dm - *(da.dm + i));
		}
	}
	if (dm_flag == 1) {
		fwrite((da.dm + k), sizeof(double), 1, opd->fp);
	}
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
	
	/* read two row then process the data */
				
	rv = psrfits_read_subint(&pf);
	if (rv) {
		if (rv == FILE_NOT_OPENED) rv = 0;
		run = 0;
		exit(0);
	}
	int *ishift;
	float fch1 = 0.0f;
	double foff = pf.hdr.df;
	int nchans = pf.hdr.nchan;
	int nbands = 1;
	double userdm  = 10.0;
	double refrf = 0.0;
	double tsamp = pf.hdr.dt;
	int maxshift;
	int nifs = pf.hdr.npol;
	float *buff[2];
	int readnext = 0, isamp, bnum, nsamp, s, c, b, indx, ns[2], soffset, ddidx;
	int ic, ixnb, ixnc, nsblk, nsmax, cpb, d, spb, x, nsout, nxb;
	float nextbaseline, realtime; 
	//for (i = 0; i < pf.hdr.nchan; i++) printf("%f\n", pf.sub.dat_freqs[i]);
	
	ishift = dmshift(fch1, foff, nchans, nbands, userdm, refrf, tsamp, pf.sub.dat_freqs);		
	
	maxshift = ishift[nchans - 1];
	
	nsblk = 256 * 2048; nsout = 32 * nchans;
	nsmax = maxshift * nifs * nchans;
	printf("%d %d %d\n", nsblk, nsout, nsmax);
	if (nsmax > nsblk) nsblk = nsmax;
	nxb = nifs * nbands;
	printf("%d\n", nxb);
	char message[80];
	float *dedisp = (float *)malloc(nxb * nsout * sizeof(float));
	float *offset = (float *)malloc(nxb * sizeof(float));
	float *tmpblk = (float *)malloc(nsout * sizeof(float));
	buff[0] = (float *)malloc(nsblk * sizeof(float));
	buff[1] = (float *)malloc(nsblk * sizeof(float));
	for (i = 0; i < nxb; i++) offset[i] = 0.0;
	d = bnum = isamp = 0;
	ic = nchans * nifs;
	nextbaseline = realtime = 0.0;
	cpb = nchans/nbands;
	if ((cpb*nbands != nchans)) error_message("silly sub-band selection!\n");
	while(1) {
		if (!readnext) {
			if (ns[bnum] = read_block(input, nbits, buff[bnum], nsblk))
		}
	}
	free(data);
	free(data1);
	free(da.data_adm);
}

