#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[])
{
	char source[24];
	int numpoint;
	FILE *fp;
	if ((fp = fopen(argv[1], "r")) == NULL) {
		exit(1);
	}
	fread(source, sizeof(char), 24, fp);
	printf("%s\n", source);
	fread(&numpoint, sizeof(int), 1, fp);
	printf("%d\n", numpoint);
}
