#include <stdio.h>
#include <stdlib.h>

int file_exist(char *filename)
{
	FILE *fp;
	if ((fp = fopen(filename, "wb")) == NULL) {
		return(0);
	} else {
		return(1);
	}
}
