#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"

int print_fits_header(fitsfile *fptr)
{
	char card[FLEN_CARD];
	int status = 0, nkeys, ii;
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
	
	return status;
}
