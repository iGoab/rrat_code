#include <stdio.h>
#include <stdlib.h>

FILE *open_file(char *filename, char *descriptor)
{
	FILE *fopen(), *fptr;
	if ((fptr = fopen(filename, descriptor)) == NULL) {
		fprintf(stderr, "Error in opening file:%s\n", filename);
		exit(1);
	}
	
	return fptr;
}
