#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./pulse_dm.h"
#include "psrfits.h"
#include "fitsio.h"

#define PI M_PI

int output_function(FILE *fp, float *data, int size)
{
	int i = 0;
	fwrite(data, sizeof(float), size, fp);
	fclose(fp);
	return(0);
}

void *thread_function(void *args)
{
	struct output_buf *opb = (struct output_buf *)args;
	int rv = output_function(opb->fp, opb->odata, opb->size);
	pthread_exit(&rv);
}

int fold_phase(struct psrfits *pf, struct data_alldm *da, struct fold_buf *fb, int nr)
{
	struct hdrinfo *hdr = &(pf->hdr);
	float period = 1 / fb->pulsar_freq;
	
	int i = 0, j = 0;
	float isam_phase = 0.0f;
	float min = 0.0f, temp = 0.0f;
	int bin_index = 0;
	fb->phase_ibin = (float *)malloc(sizeof(float) * fb->nbin);
	for (i = 0; i < fb->nbin; i++) {
		fb->phase_ibin[i] = (i + 0.5) / fb->nbin;
		//printf("%f\n", fb->phase_ibin[i]);
	}
	fb->fold_data = (float *)malloc(sizeof(float) * fb->nbin);
	fb->ibin_count = (int *)malloc(sizeof(int) * fb->nbin);
	fb->fold_x_coor = (float *)malloc(sizeof(float) * fb->nbin);
	fb->freq_point = 0.0;
	for (i = 0; i < fb->nbin; i++) {
		fb->ibin_count[i] = 0;
		*(fb->fold_x_coor + i) = (float)i /fb->nbin;
		fb->freq_point++;
		*(fb->fold_data + i) = 0.0f;
	}
	for (i = 0; i < nr; i++) {
		isam_phase = (i + 1) * hdr->dt / period - (int)((i + 1) * hdr->dt / period);
		min = isam_phase - fb->phase_ibin[0];
		if (min < 0) min = -min;
		bin_index = 0;
		//printf("%f\n", isam_phase);
		for (j = 1; j < fb->nbin; j++) {
			temp = isam_phase - fb->phase_ibin[j];
			if (temp < 0) temp = -temp;
			
			if (min > temp) { 
				min = temp;
				bin_index = j;
			}
		}
		fb->fold_data[bin_index] += *(da->data_adm + i);
		fb->ibin_count[bin_index]++;
		//printf("%f %d\n", fb->fold_data[bin_index], fb->ibin_count[bin_index]);
	}
	for (i = 0; i < fb->nbin; i++) {
		fb->fold_data[i] /= fb->ibin_count[i];
		if (fb->fold_data_max < *(fb->fold_data + i))
			fb->fold_data_max = *(fb->fold_data + i);
		if (fb->fold_data_min > *(fb->fold_data + i))
			fb->fold_data_min = *(fb->fold_data + i);
	}
	

}
/* ********************************
 * set parameter for data_alldm   *
 * ********************************/
int set_rawdata(struct data_alldm *da, int nr, double dt)
{
	int i = 0;
	for (i = 0; i < nr; i++) {
		*(da->dm_x_coor + i) = i * dt;
		if (da->dm_data_max < *(da->data_adm + i))
			da->dm_data_max = *(da->data_adm + i);
		if (da->dm_data_min > *(da->data_adm + i))
			da->dm_data_min = *(da->data_adm + i);
	}
}

/* ********************************
 * Apply zero dm to data          *
 * ********************************/
int zero_dm(struct psrfits *pf, float *data, float *data_zerodm) 
{
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0, j = 0;
	float sum = 0.0f;
	float avg = 0.0f;
	for (i = 0; i < hdr->nsblk; i++) {
		for (j = 0; j < hdr->nchan; j++) {
			sum += *(data + j + i * hdr->nchan);
		}
		
		avg = sum / hdr->nchan;
		for (j = 0; j < hdr->nchan; j++) {
			*(data_zerodm + j + i * hdr->nchan) = *(data + j + i * hdr->nchan) - avg;
			//printf("%f\n", *(data_zerodm + j + i * hdr->nchan));
		}
		
		sum = 0.0f;
		avg = 0.0f;
	}
}
 
int calc_pulse_num(struct psrfits *pf, struct pulse_num *pn, struct smooth_data *sd, int nr)
{
	struct hdrinfo	*hdr = &(pf->hdr);
	int i = 0, greater_rms_count = 0, k = 0, j = 0, index = 0;
	float square_sum = 0.0f, data_max = 0.0f;
	float flux_sum = 0.0f;
	//pn->pulse_count = 0;
	//pn->count_threshold = 0;
	for (i = 0; i < nr; i++) {
			square_sum += *(sd->data_smooth + i) 
						* (*(sd->data_smooth + i));
	}
	pn->rms = pow(square_sum / nr, 0.5);
	
	pn->threshold_value = 6 * pn->rms;
	//printf("threshold_value: %f\n", pn->threshold_value);
	for (i = 0; i < nr; i++) {
			if (pn->threshold_value < *(sd->data_smooth + i))
				pn->count_threshold++;
	}
	//printf("1\n");
	pn->greater_rms = (float *)malloc(sizeof(float) * pn->count_threshold++);
	pn->greater_rms_index = (int *)malloc(sizeof(int) * pn->count_threshold++);
	pn->single_pulse_flux = (float *)malloc(sizeof(float) * pn->count_threshold++);
	j = 0;
	for (i = 0; i < nr; i++) {
		if (pn->threshold_value < *(sd->data_smooth + i)) {
			*(pn->greater_rms + greater_rms_count) = *(sd->data_smooth + i);
			*(pn->greater_rms_index + greater_rms_count) = i;
			greater_rms_count++;
			continue;
		}
		*(sd->data_smooth + i) = 0.0f;
	}
	
	for (i = 0; i < (greater_rms_count); i++) {
		if ((*(pn->greater_rms_index + i) + 1) == *(pn->greater_rms_index + i + 1)) {
			k++;
			flux_sum += *(pn->greater_rms + i);
			//printf("index: %d value :%f\n", *(pn->greater_rms_index + i), *(pn->greater_rms + i));
			
		} else {
			if (k >= 4) {
				//printf("k: %d fluxsum: %f\n", k, flux_sum);
				*(pn->single_pulse_flux + pn->pulse_count) = flux_sum;
				pn->pulse_count++;
				//printf("pulse_count: %d\n", pn->pulse_count);
			}
			k = 0;
			flux_sum = 0.0f;
		}
	}
	j = 0;
	k = 0;
	//printf("2\n");
	pn->phase_index_pulse = (int *)malloc(sizeof(int) * pn->pulse_count);
	pn->phase_index_flux = (float *)malloc(sizeof(float) * pn->pulse_count);
	for (i = 0; i < (greater_rms_count); i++) {
		if ((*(pn->greater_rms_index + i) + 1) == *(pn->greater_rms_index + i + 1)) {
			k++;
			if (data_max <= *(pn->greater_rms + i)) {
				data_max = *(pn->greater_rms + i);
				index = *(pn->greater_rms_index + i);
			}
			
		} else {
			if (k >= 4) {
				//printf("k: %d fluxsum: %f\n", k, flux_sum);
				*(pn->phase_index_pulse + j) = index;
				j++;
				data_max = 0.0f;
				//printf("j: %d\n", j);
			}
			k = 0;
		}
	}
	//printf("3\n");
	//printf("pulse_count:%d, j:%d\n", pn->pulse_count, j);
	pn->pulse_peak_max = *(sd->data_smooth + pn->phase_index_pulse[0]);
	pn->pulse_peak_min = *(sd->data_smooth + pn->phase_index_pulse[0]);
	for (i = 0; i < pn->pulse_count; i++) {
		*(pn->phase_index_flux + i) = *(sd->data_smooth + pn->phase_index_pulse[i]);
		
		if (pn->pulse_peak_max < *(pn->phase_index_flux + i))
			pn->pulse_peak_max = *(pn->phase_index_flux + i);
		if (pn->pulse_peak_min > *(pn->phase_index_flux + i))
			pn->pulse_peak_min = *(pn->phase_index_flux + i);
		
	}
	//for (i = 0; i < pn->pulse_count; i++) 
		//printf("%d\n", *(pn->phase_index_pulse + i));
	//printf("%d\n", pn->pulse_count);
	//printf("count:%d\n", j);
	//for (i = 0; i < pn->pulse_count; i++)
		//printf("%f\n", *(pn->phase_index_flux + i));
}


/* Output the dedipersion data if needed */
int output_datatotxt(struct psrfits *pf, struct output_dmdata *opd, struct data_alldm *da)
{
	struct hdrinfo *hdr = &(pf->hdr);
	
	int i = 0;
	int nr = hdr->nsblk * pf->rows_per_file;
	if ((opd->fp = fopen(opd->filename, "wb")) == NULL) {
		printf("Cannot create txt file\n");
		exit(1);
	}
	
	for (i = 0; i < nr; i++) {
		fprintf(opd->fp, "%f\t", *(da->data_adm + i));
		if ((i + 1) % 15 == 0) {
			fprintf(opd->fp, "\n");
		}
	}
	
	fclose(opd->fp);
}

/* ***************************
 * Calculate the smooth data *
 * ***************************/
int calc_smooth(struct psrfits *pf, struct smooth_data *sd, struct data_alldm *da,int nr)  
{	
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0, j = 0;
	for (i = 0; i < nr; i++) {
		sd->smooth_count++;
		if ((i + sd->dura_point) == nr) {
			break;
		}
	}
	
	sd->data_smooth = (float *)malloc(sizeof(float) * sd->smooth_count);
	sd->smooth_x_coor = (float *)malloc(sizeof(float) * sd->smooth_count);
	for (i = 0; i < sd->smooth_count; i++) {
		for (j = 0; j < sd->dura_point; j++) {
			sd->smooth_sum += *(da->data_adm + i + j);		
		}

		*(sd->data_smooth + i) = sd->smooth_sum / sd->dura_point;
		
		*(sd->smooth_x_coor + i) = (float)i * hdr->dt;
		sd->smooth_sum = 0.0f;
		if (sd->smooth_data_max < *(sd->data_smooth + i))
			sd->smooth_data_max = *(sd->data_smooth + i);
		if (sd->smooth_data_min > *(sd->data_smooth + i))
			sd->smooth_data_min = *(sd->data_smooth + i);
	}
}         

/* ******************************************
 * calculate floating average to smooth data*
 * ******************************************/
 int calc_float_avg_smooth(struct smooth_data *fltsd, struct data_alldm *da, double dt, int nr)
 {
	 int i = 0, j = 0;
	 int k = (fltsd->dura_point - fltsd->dura_point % 2) / 2;
	 for (i = 0; i < nr; i++) {
		 if ((i + fltsd->dura_point) > nr) break;
		 for (j = i; j < (i + fltsd->dura_point); j++) {
			 fltsd->smooth_sum += *(da->data_adm + j);
		 }
		 *(fltsd->data_smooth + i + k) = fltsd->smooth_sum / fltsd->dura_point;
		 fltsd->smooth_sum = 0.0f;
		
	 }
	 
	 for (i = 0; i < k; i++) {
		 *(fltsd->data_smooth + i) = *(da->data_adm + i); 
	 }
	 
	 i = nr - 1;
	 while (i > nr - 1 -k) {
		 *(fltsd->data_smooth + i) = *(da->data_adm + i); 
		 i--;
	 }
	 
	 for (i = 0; i < nr; i++) {
		 fltsd->smooth_count++;
		 if (fltsd->smooth_data_max < *(fltsd->data_smooth + i))
			 fltsd->smooth_data_max = *(fltsd->data_smooth + i);
		 if (fltsd->smooth_data_min > *(fltsd->data_smooth + i))
			 fltsd->smooth_data_min = *(fltsd->data_smooth + i);
		 
		 *(fltsd->smooth_x_coor + i) = (float)i * dt;
	 }
 }
 
 /* **************************************
  * use gaussian function to smooth data *
  * **************************************/
int calc_gauss_smooth(struct smooth_data *gaussd, struct data_alldm *da, double dt, int nr)
{
	 int i = 0, j = 0;
	 int k = (gaussd->dura_point - gaussd->dura_point % 2) / 2;
	 for (i = 0; i < gaussd->dura_point; i++) {
		 *(gaussd->gauss_coe + i) = (exp(-0.5 * pow((i - k), 2) / (gaussd->sigma * gaussd->sigma))) /
								     pow((2 * PI * gaussd->sigma * gaussd->sigma), 0.5);
	 }

	 for (i = 0; i < nr; i++) {
		 if ((i + gaussd->dura_point) > nr) break;
		 for (j = i; j < (i + gaussd->dura_point); j++) {
			 gaussd->smooth_sum += *(gaussd->gauss_coe + j - i) * (*(da->data_adm + j));
		 }
		 *(gaussd->data_smooth + i + k) = gaussd->smooth_sum;
		 gaussd->smooth_sum = 0.0f;
	 }
	 
	 for (i = 0; i < k; i++) {
		 *(gaussd->data_smooth + i) = *(da->data_adm + i); 
	 }
	 
	 i = nr - 1;
	 while (i > nr - 1 -k) {
		 *(gaussd->data_smooth + i) = *(da->data_adm + i); 
		 i--;
	 }
	 
	 for (i = 0; i < nr; i++) {
		 gaussd->smooth_count++;
		 if (gaussd->smooth_data_max < *(gaussd->data_smooth + i))
			 gaussd->smooth_data_max = *(gaussd->data_smooth + i);
		 if (gaussd->smooth_data_min > *(gaussd->data_smooth + i))
			 gaussd->smooth_data_min = *(gaussd->data_smooth + i);
		 
		 *(gaussd->smooth_x_coor + i) = (float)i * dt;
	 }
 }
 
/* **************************************
 * use exponential smooth to smooth data*
 * **************************************/
int calc_dbexp_smooth(struct smooth_data *double_expsd, struct data_alldm *da, int nr)
{
	int i = 0;
	*(double_expsd->data_smooth + 0) = *(da->data_adm + 0);
	*(double_expsd->data_trend + 0) = 0.0f;
	*(double_expsd->smooth_x_coor + 0) = 0.0f;
	double_expsd->smooth_data_max = *(double_expsd->data_smooth + 0);
	double_expsd->smooth_data_min = *(double_expsd->data_smooth + 0);
	double_expsd->smooth_count++;
	for (i = 1; i < nr - 20; i++) {
		*(double_expsd->data_smooth + i) = double_expsd->exp_alpha * (*(da->data_adm + i)) 
										+ (1 - double_expsd->exp_alpha) * (*(double_expsd->data_smooth + i - 1) 
										+ *(double_expsd->data_trend + i -1));
		*(double_expsd->data_trend + i) = double_expsd->exp_beta * (*(double_expsd->data_smooth + i) - 
										*(double_expsd->data_smooth + i - 1)) + (1 - double_expsd->exp_beta) 
										*(*(double_expsd->data_trend + i - 1));
		*(double_expsd->smooth_x_coor + i) = (float)i;
		if (double_expsd->smooth_data_max < *(double_expsd->data_smooth + i))
			double_expsd->smooth_data_max = *(double_expsd->data_smooth + i);
		if (double_expsd->smooth_data_min > *(double_expsd->data_smooth + i))
			double_expsd->smooth_data_min = *(double_expsd->data_smooth + i);
		double_expsd->smooth_count++;
 	}
}
 
/* *********************
 * Calculate fold data *
 * *********************/                                                                                                                                                                                    
int calc_fold(struct psrfits *pf, struct fold_buf *fb, struct data_alldm *da, int nr) 
{
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0, j = 0;
	int n = (int)(nr / fb->freq_point);
	for (i = 0; i < fb->freq_point; i++) {
		*(fb->fold_x_coor + i) = (float)i * hdr->dt;
	}
	for (j = 0; j < n; j++)
		for (i = 0; i < fb->freq_point; i++) {
			*(fb->fold_data + i) += *(da->data_adm + i + j * fb->freq_point) / n ; 
			if (j == n - 1) {
				if (fb->fold_data_max < *(fb->fold_data + i))
					fb->fold_data_max = *(fb->fold_data + i);
				if (fb->fold_data_min > *(fb->fold_data + i))
					fb->fold_data_min = *(fb->fold_data + i);
			}
		}	
}

/* **********************************
 * Calculate fold after smooth data *
 * **********************************/
int calc_smooth_fold(struct psrfits *pf, struct smooth_data *sd, struct fold_buf *fb)
{
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0, j = 0;
	int n = (int)(sd->smooth_count / fb->freq_point);
	for (i = 0; i < fb->freq_point; i++) {
		*(fb->fold_x_coor + i) = (float)i * hdr->dt;
	}
	for (j = 0; j < n; j++)
		for (i = 0; i < fb->freq_point; i++) {
			*(fb->fold_data + i) += *(sd->data_smooth + i + j * fb->freq_point) / n ; 
			if (j == n - 1) {
				if (fb->fold_data_max < *(fb->fold_data + i))
					fb->fold_data_max = *(fb->fold_data + i);
				if (fb->fold_data_min > *(fb->fold_data + i))
					fb->fold_data_min = *(fb->fold_data + i);
			}
		}	
}

/* *****************************
 * Calculate smooth after fold *
 * *****************************/
int calc_fold_smooth(struct psrfits *pf, struct fold_buf *fb, struct smooth_data *sd)
{
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0, j = 0; 
	sd->dura_point = 200;
	int n = (int)(fb->freq_point / sd->dura_point);
	for (i = 0; i < fb->freq_point; i++) {
		sd->smooth_count++;
		if ((i + sd->dura_point) > fb->freq_point) break;
	}

	sd->data_smooth = (float *)malloc(sizeof(float) * sd->smooth_count);
	sd->smooth_x_coor = (float *)malloc(sizeof(float) * sd->smooth_count);
	for (i = 0; i < sd->smooth_count; i++) {
		for (j = 0; j < sd->dura_point; j++) {
			sd->smooth_sum += *(fb->fold_data + i + j);		
		}
		//for (j = 0; j < sd->dura_point; j++)
		//	sd->smooth_sum += *(fb->fold_data + j + i * sd->dura_point);

		*(sd->data_smooth + i) = sd->smooth_sum;
		
		*(sd->smooth_x_coor + i) = (float)i * hdr->dt;
		sd->smooth_sum = 0.0f;
		if (sd->smooth_data_max < *(sd->data_smooth + i))
			sd->smooth_data_max = *(sd->data_smooth + i);
		if (sd->smooth_data_min > *(sd->data_smooth + i))
			sd->smooth_data_min = *(sd->data_smooth + i);
	}
}

/* **********************************************************
 * Return the delay in seconds caused by dispersion, given  *
 * a Dispersion Measure (dm) in cm-3 pc, and the emitted    *
 * frequency (freq_emitted) of the pulsar in MHz.           *
 * **********************************************************/
double delay_from_dm(double dm, float freq_emitted)
{
    return dm / (0.000241 * freq_emitted * freq_emitted);
}

/* **************************************
 * Use the data_freq to calculate delay *
 * *************************************/
int calc_dm(struct psrfits *pf, int **idelays, double *dm, int ndm) 
{
	struct subint *sub = &(pf->sub);
	struct hdrinfo *hdr = &(pf->hdr);
	double hifreq = 0.0f, dtmp = 0.0f, sub_delays = 0.0f;
	int flag_check_freq_seq = 0;
	int jj = 0, ii = 0;;
	double *chan_delays;
	chan_delays = (double *)malloc(sizeof(double) * hdr->nchan);
	for (ii = 0; ii < ndm; ii++) {
		if (sub->dat_freqs[0] > sub->dat_freqs[hdr->nchan - 1]) {
			hifreq = sub->dat_freqs[0] - hdr->df * 0.5;
			dtmp = hifreq + 0.5 * hdr->df;
			sub_delays = delay_from_dm(dm[ii], dtmp);
			for (jj = 0 ; jj < hdr->nchan ; jj++) {
				chan_delays[jj] = delay_from_dm(dm[ii], sub->dat_freqs[jj]);
				chan_delays[jj] -= sub_delays;
				*(idelays[ii] + jj) = (int)rint(chan_delays[jj] / hdr->dt);
				//printf("dm:%f\tidelays:%d\n", dm, idelays[jj]);
			}
		} else {
			flag_check_freq_seq = 1;
			hifreq = sub->dat_freqs[hdr->nchan - 1] - hdr->df * 0.5;
			dtmp = hifreq + 0.5 * hdr->df;
			sub_delays = delay_from_dm(dm[ii], dtmp);
			
			// Determine the dispersion delays and convert them
			// to offsets in units of sample times
			for (jj = hdr->nchan - 1; jj >= 0; jj--) {
				chan_delays[jj] = delay_from_dm(dm[ii], sub->dat_freqs[jj]);
				chan_delays[jj] -= sub_delays;
				*(idelays[ii] + jj) = (int)rint(chan_delays[jj] / hdr->dt);
				//printf("dm:%f\tidelays:%d\n", dm, idelays[jj]);
			}
		}
	}
	free(chan_delays);
	return flag_check_freq_seq;
}

/* **************************************************************
 * Get the stokes I of data if the poln_type is AABBCRCI or AABB*
 * **************************************************************/
int get_AABB_I(struct psrfits *pf, float *data_out)
{
	struct subint *sub = &(pf->sub);
	struct hdrinfo *hdr = &(pf->hdr);
	unsigned char *data, *data1;

	int i = 0, j = 0;
	for (i = 0; i < hdr->nsblk; i++) {
		data = sub->rawdata + i * hdr->nchan * hdr->npol;
		data1 = data + hdr->nchan;
		for (j = 0; j < hdr->nchan; j++, data++, data1++)
			*(data_out + j + i * hdr->nchan) = (float)*data + *data1;
		
	}
	
}

/********************************************************************
 * Apply zerodm to process the data 								*
 * ******************************************************************/
 int zerodm(struct psrfits *pf, float *data_in)
 {
	struct hdrinfo *hdr = &(pf->hdr);
	float sum = 0.0f;
	float avg = 0.0f;
	int i = 0, j = 0;
	for (i = 0; i < hdr->nsblk; i++) {
		for (j = 0; j < hdr->nchan; j++) {
			sum += *(data_in + j + i * hdr->nchan);
		}
		
		avg = sum / hdr->nchan;
		for (j = 0; j < hdr->nchan; j++) {
			*(data_in + j + i * hdr->nchan) -= avg;
			//printf("%f\n", *(data_zerodm + j + i * hdr->nchan));
		}
		sum = 0.0;
		avg = 0.0;
	}
 }
 
/* *******************************************************************
 * Get the stokes I of data if the poln_type is IQUV or AA or AA + BB*
 * *******************************************************************/
int get_only_I(struct psrfits *pf, float *data_out)
{
	struct subint *sub = &(pf->sub);
	struct hdrinfo	*hdr = &(pf->hdr);

	unsigned char *data;
	int i = 0, j = 0;
	for (i = 0; i < hdr->nsblk; i++) {
		data = sub->rawdata + i * hdr->nchan * hdr->npol;
		for (j = 0; j < hdr->nchan; j++, data++) {
			*(data_out + j + i * hdr->nchan) = (float)*data;
		}	
	}
}

/* *************************************
 *  Using the dm to delay the data     *
 * ************************************/
int delay_data(struct psrfits *pf, float *data_in, float *data_in1, struct data_alldm *da, int flag_check_freq, int ndm)
{
	struct hdrinfo	*hdr = &(pf->hdr);
	int i = 0, j = 0, k = 0;
	if (flag_check_freq) {
		for (i = 0; i < hdr->nsblk; i++) {
			for (j = hdr->nchan - 1; j >= 0; j--) {
				if ((i + *(da->idelays[ndm] + j)) >= hdr->nsblk) {
					break;
				}
				*(da->data_adm + i) += *(data_in + j + (i + *(da->idelays[ndm] + j)) * hdr->nchan);
			}
		}
		
		for (i = hdr->nchan - 1; i >= 0; i--)
			for (j = 0; j < *(da->idelays[ndm] + i); j++) {
				*(da->data_adm + hdr->nsblk - 1 -j) += *(data_in1 + i + j * hdr->nchan);
			}
	} else {
		
		for (i = 0; i < hdr->nsblk; i++) {
			for (j = 0; j < hdr->nchan; j++) {
				if ((i + *(da->idelays[ndm] + j)) >= hdr->nsblk) {
					break;
				}
				*(da->data_adm + i) += *(data_in + j + (i + *(da->idelays[ndm] + j)) * hdr->nchan);
				//printf("%d\n", *(da->data_adm + i));
			}
		}
		
		for (i = 0; i < hdr->nchan; i++)
			for (j = 0; j < *(da->idelays[ndm] + i); j++) {
				*(da->data_adm + hdr->nsblk - 1 -j) += *(data_in1 + i + j * hdr->nchan);
			}

	}
	
	/*if (flag_check_freq) {
		if (row != 0) {
			for (i = hdr->nchan - 1; i >= 0; i--)
				for (j = 0; j < *(da->idelays[ndm] + i); j++) {
					*(da->data_adm + hdr->nsblk - 1 - j + (row - 1) * hdr->nsblk) += *(data_in + i + j * hdr->nchan);
					
					//printf("%f\n", *(data_in + i + j * hdr->nchan));
				}
		}
		for (i = 0; i < hdr->nsblk; i++) {
			for (j = hdr->nchan - 1; j >= 0; j--) {
				if ((i + *(da->idelays[ndm] + j)) >= hdr->nsblk) {
					break;				
				}	
				*(da->data_adm + i + row * hdr->nsblk) += 
						*(data_in + j + (i + *(da->idelays[ndm] + j)) * hdr->nchan);
				//printf("%f\n", *(data_in + j + (i + da->idelays[j]) * hdr->nchan));
			}
		}
	} else {
		if (row != 0) {
			for (i = 0; i < hdr->nchan; i++) 
				for (j = 0; j < *(da->idelays[ndm] + i); j++) {
					*(da->data_adm + hdr->nsblk - 1 - j + (row - 1) * hdr->nsblk) += *(data_in + i + j * hdr->nchan);
					//printf("%d: %f\n", i, *(da->data_adm + hdr->nsblk - 1 - j + (row - 1) * hdr->nsblk));
				}
		}
		
		for (i = 0; i < hdr->nsblk; i++) {
			for (j = 0; j < hdr->nchan; j++) {

				if ((i + *(da->idelays[ndm] + j)) >= hdr->nsblk) {
					break;
					//if (ndm == 132)
					//printf("%d\n", i + da->idelays[j]);
				}

				*(da->data_adm + i + row * hdr->nsblk) += 
					*(data_in + j + (i + *(da->idelays[ndm] + j)) * hdr->nchan);

			}
		}
	}*/
}

/* ******************************************
 * To calculate the average of the data and *
 * use the data to substract the average    *
 * ******************************************/
 int calc_data_sub_avg(struct data_alldm *da, double dt, int nr)
 {
	 float avg_sum = 0.0f;
	 int i = 0;
   	 for (i = 0; i < nr; i++) {
		 avg_sum += *(da->data_adm + i);
	 }
	
	da->avg_flux = avg_sum / nr;
	
	for (i = 0; i < nr; i++) {
		*(da->data_adm + i) = (*(da->data_adm + i) - da->avg_flux) / 1000.0f;
		if (da->dm_data_max < *(da->data_adm + i)) {
			da->dm_data_max = *(da->data_adm + i);

		}
		if (da->dm_data_min > *(da->data_adm + i)) {
			da->dm_data_min = *(da->data_adm + i);
		}
		
		*(da->dm_x_coor + i) = (float)i * dt;
	}	
 }
 


float calc_dm_step(struct psrfits *pf)
/* Use the center frequency of the data */
/* and the bandwidth to calculate the 	*/
/* dm step 							    */
{
	struct hdrinfo *hdr = &(pf->hdr);
	return 1.205 * pow(10, -7) * hdr->dt * pow(hdr->fctr, 3) / hdr->BW;
}

/* **********************************************
 * Use after dedispersion data to calculate the *
 *  rms 										*
 * **********************************************/
float calc_rms(float *data_adm, float *rms, int nr)
{
	int i = 0;
	float square_sum = 0.0f, square_sum_avg = 0.0f;
	for (i = 0; i < nr; i++) {
		square_sum += *(data_adm + i) * (*(data_adm + i));
	}
	square_sum_avg = square_sum / nr;
	*rms = pow(square_sum_avg, 0.5);
}

