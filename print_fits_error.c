#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
void print_fits_error(int status)
{
	if (status) {
		fits_report_error(stderr, status);
	}
}
