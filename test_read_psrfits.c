#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <malloc.h>
#include "fitsio.h"
#include "psrfits.h"

int main(int argc, char *argv[])
{
	struct psrfits pf;
	psrfits_set_files(&pf, argc - optind, argv + optind);
	
	int rv = psrfits_open(&pf);
	
	int ii = 0;
    //printf("%s\n", pf.hdr.backend);
    /*pf.sub.dat_freq = (float *)malloc(sizeof(float) * pf.hdr.nchan);

	pf.sub.dat_wts = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    for (ii = 0; ii < pf.hdr.nchan; ii++) {
	     *(pf.sub.dat_freq + ii) = ii;
	}
    for (ii = 0; ii < pf.hdr.nchan; ii++) {
	    printf("%f\n", pf.sub.dat_freq[ii]);
	}
	pf.sub.dat_offs = (float *)malloc(sizeof(float) 
					  * pf.hdr.nchan * pf.hdr.npol);
	pf.sub.dat_scl = (float *)malloc(sizeof(float)
	                  * pf.hdr.nchan * pf.hdr.npol);
	pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	
	/*(while((rv = psrfits_read_subint(&pf)) == 0) {
		 printf("Read subint (file %d, row %d/%d)\n",
		        pf.filenum, pf.rownum - 1, pf.rows_per_file);
		 if ((pf.rownum -1) == pf.rows_per_file) { break; }
	}*/
	/*for (ii = 0; ii < 200; ii++) {
	    printf("%c\n", *(pf.sub.rawdata + ii));	
	}*/
	/*free(pf.sub.dat_wts);
	free(pf.sub.dat_freq);
	free(pf.sub.dat_scl);
	free(pf.sub.dat_offs);
	free(pf.sub.rawdata);
	free(pf.sub.data);*/
	if(rv) { fits_report_error(stderr, rv); }
	exit(0);
}
