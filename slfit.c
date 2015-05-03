#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*********************************************************************
 * Fits the data passed down in arrays x(ndata) y(ndata) to the      *
 * straight line y = a + bx by minimising chi**2. Standard deviations*
 * in y are passed down by sig(ndata) and can be weighted into       *
 * the fit if the logical switch weight is on. a and b together      *
 * with their uncertainties siga and sigb are returned.               *
 * adapted from the fit routine given in numerial recipes to be used *
 * in the psrflux software.                                          *
 * *******************************************************************/
void slfit(float *x, float *y, float *z, int ndata,
           float *sig, int weight, float *a, float *b, 
           float *siga, float *sigb)
{
	float sx = .0f;
	float sy = .0f;
	float st2 = .0f;
	float ss = .0f;
	*b = .0f;
	int q = 0;
	float sigdat = 0.0f;
	float t = 0.0f;
	float wt = 0.0f;
	float chi2 = 0.0f;
	float sxoss = 0.0f;
	int useweight = weight;
	int i = 0;
/**
 * for upper limits opt for an unweighted fit 
 * and use only three quarters the flux values
 * */
	for (i = 0; i < ndata; i++) {
		*(z + i) = *(y + i);
		if (*(sig + i) == -999.0) {
		  useweight = 0;
		  *(z + i) = 0.75 * (*(y + i));
	    }
	}
	
	if (useweight) {
		for (i = 0; i < ndata; i++) {
			wt = 1.0/(*(sig + i) * (*(sig + i)));
			ss = ss + wt;
			sx = sx + *(x + i) * wt;
			sy = sy + *(z + i) * wt;
		}
	} else {
		for (i = 0; i < ndata; i++) {
			sx = sx + *(x + i);
			sy = sy + *(z + i);
		}
		ss = (float)ndata;
	}
	sxoss = sx/ss;
	if (useweight) {
		for (i = 0; i < ndata; i++) {
			t = (*(x + i) - sxoss) / *(sig + i);
			*b = *b + t * (*(z + i)/ * (sig + i));
			st2 = st2 + t * t;
		}
	} else {
		for (i = 0 ; i < ndata; i++) {
			t = (*(x + i) - sxoss);
			st2 = st2 + t * t;
			*b = *b + t * (*(z + i));
		}
	}
	*b = *b /st2;
	*a = (sy -sx * (*b)) / ss;
	*siga = sqrt((1.0 + sx * sx/(ss * st2))/ss);
	*sigb = sqrt(1.0/st2);
	chi2 = 0.0;
	if (!useweight) {
		for (i = 0; i < ndata; i++) {
			chi2 = pow((chi2 + ((*(z + i)) - *a - *b * (*(x + i)))), 2);
		}
		q = 1;
		sigdat = 0.0f;
		if (ndata > 2) sigdat = sqrt(chi2 / (ndata - 2));
		*siga = *siga * sigdat;
		*sigb = *sigb * sigdat;
	}
}
