#include <stdio.h>
#include <stdlib.h>

void error_message(char *message)
{
	fprintf(stderr, "ERROR: %s", message);
	exit(1);
}
