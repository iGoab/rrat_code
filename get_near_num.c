#include <stdlib.h>
#include <stdio.h>

float get_near_num(float var, float *var_fix, float *var_fix_before) {
	if (var_fix < var_fix_before) {
		*var_fix_before = *var_fix;
		return (var);
	}
}
