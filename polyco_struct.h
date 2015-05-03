/* polyco struct */
struct polyco {
	 char psr[15];          // The pulsar source name
	 double fmjd;   		// The modified julian date
	 double rphase;			// Reference phase
	 double f0;				// Reference rotation frequency
	 double c[15];			// Coefficients
	 float rf;				//
	 int mjd;				// Modified julian date
	 int nsite;
	 int nmin;
	 int nc;				// Number of coefficients
	 int used;
	 long long rphase_int;
};
