#include <stdio.h>
/*
   return the delay in seconds between two sky frequencies f1 and f2 (MHz)
   N.B. constants of proprotionality derived from e^2 pcm2 / (2pi c me) where
   the elementary charge is assumed to be in Gaussian units: 4.8032068e-20
   and pc2m = 3.0856776e16 is the conversion between metres and parsecs
*/
double dmdelay(float f1, float f2, double dm)
{
	return(4148.741601 * ((1.0/f1/f1) - (1.0/f2/f2)) * dm);
}
