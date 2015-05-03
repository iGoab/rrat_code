#include "psrfits.h"
/* struct to store the dedispersion data */
struct data_alldm
{
	int numpoint;
	double *dm;                    // to store the dm
	int ndm;
	double real_dm;
	float dm_data_max;           // biggiest value in data
	float dm_data_min;
	float avg_flux;
	float *dm_x_coor;
	int **idelays;                //to store the delay
	float *data_adm;	         // to store the dm data
};

/* Struct to store the fold_data*/
struct fold_buf
{
	int freq_point;
	int *ibin_count;
	int nbin;
	float fold_data_max;
	float fold_data_min;
	double pulsar_freq;
	float *fold_data;
	float *fold_x_coor;
	float *phase_ibin;
};

/* struct to store the smooth data */
struct smooth_data
{
	int sum_point;
	int dura_point;               
	int smooth_count;            
	int sigma;                  // the gaussian function parameter
	float exp_alpha;				// the exponential alpha coefficient
	float exp_beta;					// the exponential beta coefficient
	float exp_gamma;				// the exponential gamma coefficient
	float smooth_data_max;
	float smooth_data_min;
	float smooth_sum;
	float *data_trend;				// the exponential trend data
	float *data_seasonality;		// the exponential seasonlity data
	float *data_smooth;
	float *smooth_x_coor;
	float *gauss_coe;           // the gaussian coefficients
};

struct output_buf
{
	FILE *fp;
	char filename[256];
	float *odata;
	int size;
};
/* struct to define the output file pointer */
struct output_dmdata
{
	FILE *fp;
	char filename[256];
};

/* struct to store the data greater rms */
/* and the number of pulses             */
struct pulse_num
{
	int pulse_count;
	int count_threshold;
	float threshold_value;
	float rms;
	float *single_pulse_flux;
	float *greater_rms;
	float *phase_index_flux;
	float *phase_coor;
	float pulse_peak_max;
	float pulse_peak_min;
	int *phase_index_pulse;
	int *greater_rms_index;
};

struct pulse_draw
{
	struct data_alldm *da;
	struct pulse_num *pn;
};

void usage();
void *thread_function(void *args);
void free_data(struct data_alldm *da, struct smooth_data *sd, struct fold_buf *fb, struct pulse_num *pn);
double delay_from_dm(double dm, float freq_emitted);
int delay_data(struct psrfits *pf, float *data_in, float *data_in1, struct data_alldm *da, int flag_check_freq, int ndm);
int calc_dm(struct psrfits *pf, int **idelays, double *dm, int ndm);
int get_only_I(struct psrfits *pf, float *data_out);
int get_AABB_I(struct psrfits *pf, float *data_out);
int draw_flux_pic(float *data_adm, float *x_coor, float x_min, float x_max, float y_min, float y_max, double dt, int num);
int draw_dmpic(struct psrfits *pf, struct pulse_draw *pd, float dm_max, int ndm, int nr);
int calc_fold(struct psrfits *pf, struct fold_buf *fb, struct data_alldm *da, int nr);
int calc_smooth(struct psrfits *pf, struct smooth_data *sd, struct data_alldm *da, int nr);
int calc_smooth_fold(struct psrfits *pf, struct smooth_data *sd, struct fold_buf *fb);
int calc_fold_smooth(struct psrfits *pf, struct fold_buf *fb, struct smooth_data *sd);
int output_datatotxt(struct psrfits *pf, struct output_dmdata *opd, struct data_alldm *da);
int zero_dm(struct psrfits *pf, float *data, float *data_zerodm);
int draw_flux_dm_sample(struct psrfits *pf, struct pulse_draw *pd, float dm_max, int ndm, int nr);
float calc_dm_step(struct psrfits *pf);
float calc_rms(float *data_adm, float *rms, int nr);
