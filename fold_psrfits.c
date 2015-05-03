/* fold_psrfits.c
 * 
 * Fold PSRFITS search data into PSRFITS folded format.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <getopt.h>
#include <signal.h>
#include <pthread.h>
#include "polyco.h"
#include "fold.h"
#include "psrfits.h"

#define nthread 4
#define num_bin 256
/* Signal handler */
int run = 1;
void cc(int sig) { run = 0;}
int main(int argc, char *argv[])
{
	psrfits pf;
	sprintf(pf.filename, "/home/igoab/igoab/guppi_56590_B1133+16_0003_0001.fits");
	int rv, i = 0;
	rv = psrfits_open(&pf);
	if (rv) { fits_report_error(stderr, rv); exit(1);}
	
	/* Alloc sub data */
	pf.sub.dat_freq = (float *)malloc(sizeof(float) * pf.hdr.nchan );
	pf.sub.dat_wts = (float *)malloc(sizeof(float) * pf.hdr.nchan);
	pf.sub.dat_offs = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
	pf.sub.dat_scl = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint); 
	pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	
	psrfits_read_subint(&pf);
	/* Alloc total fold buf */
	foldbuf fb;
	fb.nchan = pf.hdr.nchan;
	fb.npol = pf.hdr.npol;
	fb.nbin = num_bin;
	malloc_foldbuf(&fb);
	clear_foldbuf(&fb);
	
	/* Set up thread management */
	pthread_t *thread_id;
	fold_args *fargs;
	thread_id = (pthread_t *)malloc(sizeof(pthread_t) * nthread);
	fargs = (fold_args *)malloc(sizeof(fold_args) * nthread);
	for (i = 0; i < nthread; i++) {
		thread_id[i] = 0;
		fargs[i].data = (char *)malloc(
				sizeof(char)*pf.sub.bytes_per_subint);
		fargs[i].fb = (foldbuf *)malloc(sizeof(foldbuf));
		fargs[i].fb->nbin = num_bin;
		fargs[i].fb->nchan = pf.hdr.nchan;
		fargs[i].fb->npol = pf.hdr.npol;
		fargs[i].nsamp = pf.hdr.nsblk;
		fargs[i].tsamp = pf.hdr.tbin;
		malloc_foldbuf(fargs[i].fb);
		clear_foldbuf(fargs[i].fb);
		fargs[i].scale = (float *)malloc(sizeof(float) 
				* pf.hdr.nchan * pf.hdr.npol);
		fargs[i].offset = (float *)malloc(sizeof(float)
			    * pf.hdr.nchan * pf.hdr.npol);
	}
	
	/* Main loop */
	rv = 0;
	int imjd;
	double fmjd, fmjd0 = 0, fmjd_next = 0, fmjd_epoch;
	double offs0 = 0, offs1 = 0;
	int first = 1, subcount = 0;
	int cur_thread = 0;
	signal(SIGINT, cc);
}

