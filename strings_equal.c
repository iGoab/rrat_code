#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int strings_equal(char *str1, char *str2)
{
	if (strcmp(str1, str2) == 0) {
		return 1;
	} else {
		return 0;
	}
}
