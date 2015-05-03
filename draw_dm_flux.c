#include "psrfits.h"
#include "pulse_dm.h"
#include "/home/igoab/pulsar_software/pgplot/cpgplot.h"

float max(float a, float b)
{
	if ( a < b) return (b);
	else return (a);
}

float min(float a, float b)
{
	if (a < b) return (a);
	else return (b);
}

/* ****************************************
 *  Use pgplot draw a dedispersion picture*
 * ****************************************/
/*int draw_dmpic(struct psrfits *pf, struct pulse_draw *pd, float dm_max, int ndm, int nr)
{
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0;
	if (pd->da[0].dm == 0.0) {
		if (cpgbeg(0, "/xs", 1, 1) != 1)
			return EXIT_FAILURE;
		cpgpage();
		cpgenv(-(pd->pn[0].pulse_peak_max * pd->pn[0].pulse_peak_max * 10), nr * hdr->dt, -(pd->pn[0].pulse_peak_max * pd->pn[0].pulse_peak_max * 10), dm_max, 1, 1);
		cpgsfs(2);
		
	}

	int j = 0;
	for (j = 0; j < ndm; j++)
		for (i = 0; i < pd->pn[j].pulse_count; i++) {
			
			cpgcirc(*(pd->pn[j].phase_index_pulse + i) * hdr->dt, pd->da[j].dm, *(pd->pn[j].phase_index_flux + i) * 10);
		}
	cpgbbuf();

}*/

/* use pgplot draw flux in different dm at different sample in gray plot */
int draw_flux_dm_grayplot(struct psrfits *pf, struct pulse_draw *pd, float dm_max, int ndm, int nr)
{
	
}

/* use pgplot draw flux in different dm at the same sample */
/*int draw_flux_dm_sample(struct psrfits *pf, struct pulse_draw *pd, float dm_max, int ndm, int nr)
{
	struct hdrinfo *hdr = &(pf->hdr);
	int i = 0, j = 0, key = 0;
	float *ddm, *flux, temp = 0.0f;
	ddm = (float *)malloc(sizeof(float) * ndm);
	flux = (float *)malloc(sizeof(float) * ndm);
	for (i = 0; i < pd->pn[0].pulse_count; i++) {
		if (temp < *(pd->pn[0].phase_index_flux + i)) {
			temp = *(pd->pn[0].phase_index_flux + i);
			key = i;
		} 
	}
	
	temp = 0.0f;
	
	for (i = 0; i < ndm; i++) {
		*(ddm + i) = pd->da[i].dm;
		*(flux + i) = *(pd->pn[i].phase_index_flux + key);	
		printf("ddm[%d]:%f flux[%d]:%f\n", i, *(ddm + i), i, *(flux + i));
		if (temp < *(flux + i)) {
			temp = *(flux + i);
		}
	}
	
	if (cpgbeg(0, "/xs", 1, 1) != 1)	
		return EXIT_FAILURE;

    cpgenv(0.0, dm_max, 0.0, temp, 0, 1);
    //cpgenv(0.0, 256.0, 0.0, 0.5, 0, 1);
	cpgline(ndm, ddm, flux);
	cpgend();
	//cpgsubp(-2, 1);
	
	//cpgpanl(1, 1);
	//cpgpage();
	//cpgenv(0.0, fb->freq_point / fb->nbin, 0.0, 0.5, 0, 1);
	//cpgpt(pn->pulse_count, pn->phase_coor, pn->phase_index_flux, -5);
	
}*/
int draw_fb_pulse_dis(struct psrfits *pf, struct fold_buf *fb, struct pulse_num *pn)
{
	struct hdrinfo *hdr = &(pf->hdr);
	float period = 1 / fb->pulsar_freq;
	
	int i = 0, j = 0;
	float isam_phase = 0.0f;
	float min = 0.0f, temp = 0.0f, max = 0.0f;
	int bin_index = 0;
	
	pn->phase_coor = (float *)malloc(sizeof(float) * pn->pulse_count);
	for (i = 0; i < pn->pulse_count; i++) {
		isam_phase = (*(pn->phase_index_pulse + i) + 1) * hdr->dt / period - (int)((*(pn->phase_index_pulse + i) + 1) * hdr->dt / period);
		min = isam_phase - fb->phase_ibin[0];
		if (min < 0) min = -min;
		bin_index = 0;
		//printf("%f\n", isam_phase);
		for (j = 1; j < fb->nbin; j++) {
			temp = isam_phase -fb->phase_ibin[j];
			if (temp < 0) temp = -temp;
			
			if (min > temp) { 
				min = temp;
				bin_index = j;
			}
		}
		*(pn->phase_index_pulse + i) = bin_index;
		//printf("%f %d\n", fb->fold_data[bin_index], fb->ibin_count[bin_index]);
		//printf("%d\n", *(pn->phase_index_pulse + i));
	}
	min = *(pn->phase_index_flux + 0);
	max = *(pn->phase_index_flux + 0);
	for (i = 1; i < pn->pulse_count; i++) {
		if (min > *(pn->phase_index_flux + i))
			min = *(pn->phase_index_flux + i);
		if (max < *(pn->phase_index_flux + i))
			max = *(pn->phase_index_flux + i);
	}
	
	for (i = 0; i < pn->pulse_count; i++) {
		*(pn->phase_coor + i) = (float)*(pn->phase_index_pulse + i) / fb->nbin;
		*(pn->phase_index_flux + i) -= min;
	}
	
	if (cpgbeg(0, "/xs", 1, 2) != 1)	
		return EXIT_FAILURE;
	//cpgsubp(1, 2);
	//cpgpanl(1, 1);
	//cpgpage();
	//cpgbbuf();
    cpgenv(0.0, fb->freq_point / fb->nbin, fb->fold_data_min, fb->fold_data_max, 0, 1);
    //cpgenv(0.0, 256.0, 0.0, 0.5, 0, 1);
	cpgline(fb->freq_point, fb->fold_x_coor, fb->fold_data);
	//cpgsubp(-2, 1);
	
	//cpgpanl(1, 1);
	//cpgpage();
	cpgenv(0.0, fb->freq_point / fb->nbin, 0.0, 0.5, 0, 1);
	cpgpt(pn->pulse_count, pn->phase_coor, pn->phase_index_flux, -5);
	
}

/* ******************************************
 *  Use pgplot to draw a flux picture       *
 * *****************************************/
int draw_flux_pic(float *data_adm, float *x_coor, float x_min, float x_max, float y_min ,float y_max, double dt, int num) 
/* data_adm is the data after dedispersion */
/* x_coor is the x point array */
/* x_min is the small x_coordinate */
/* x_max is the great y_coordinate */
/* y_min and y_max is the same to x_coordinate */
{
	float xl_co = 0.0f, xr_co = 0.0f;
	float yl_co = 0.0f, yr_co = 0.0f;
	char ch_scalar;				//the character typed by the user
	float x_cu[2] = {0.0f, 0.0f};   // the world of x coordinate of the cursor
	float y_cu[2] = {0.0f, 0.0f};   // the world of y coordinate of the cursor
	float xref = 0.0f, yref = 0.0f; // xref the x anchor point yref the y anchor point
	float xmin = 0.0f;
	float inc = (float)num;
	float width = 0.0f;
	float xmax = 0.0f + inc;
	/* Call PGBEG to initiate PGPLOT and open the output device; PGBEG
	 * will prompt the user to supply the device name and type. */
	if (cpgbeg(0, "/xs", 1, 1) != 1)	
		return EXIT_FAILURE;
		
	cpgpage();
	cpgenv(x_min, x_max * dt, y_min, y_max, 0, 1);
   
	cpgline(x_max, x_coor, data_adm);

	/*
	 * Finally, call PGEND to terminate things properly.
	 */
	while(TRUE) {
		cpgband(0, 0, xref, yref, x_cu, y_cu, &ch_scalar);
		
		switch(ch_scalar) {
			case '>':
					  cpgenv(xmin * dt, xmax * dt, y_min, y_max, 0, 1);
					  cpgline(x_max, x_coor, data_adm);
					  xmin += inc;
			          xmax = xmin + inc;
					  if (xmax > x_max) {
						cpgend();
				 		return EXIT_SUCCESS;
				  	 }
					  break;
			case '<':					 
					 xmin -= inc;
			         xmax -= inc;
					 if (xmin < 0.0f) {
						cpgend();
				 		return EXIT_SUCCESS;
				  	 }
				  	 cpgenv(xmin * dt, xmax * dt, y_min, y_max, 0, 1);
					 cpgline(x_max, x_coor, data_adm);
					  break;
			case 'r':
					 xmin = 0.0f;
				     xmax = xmin + inc;
			         cpgenv(x_min * dt, x_max * dt, y_min, y_max, 0, 1);
			         cpgline(x_max, x_coor, data_adm);
			         break;
			case '+':
					printf("Please input the minimum and maximum of x axis\n");
					printf("the minimum point is %f, the maximum point is %f\n", x_min, x_max);
					scanf("%f %f", &xmin, &xmax);
					cpgenv(xmin * dt, xmax * dt, y_min, y_max, 0, 1);
					cpgline(x_max, x_coor, data_adm);
					break;
			case 'h':
					puts("> right shift the point");
					puts("< left shift the point");
					puts("r reset the axis to the minimum and the maximum");
					puts("+ set the axis to between two point");
					puts("q quit pgplot");
					break;
			case 'q': cpgend();
					  return EXIT_SUCCESS;
					  break;
		}
	}
}

/******************************************
 * Use pgplot to draw a gray-scale map    *
 * ****************************************/
int draw_gray_map(float *arr, int nr, int ndm, int dt)
{
	int idim = 0, jdim = 0, i1 = 0, i2 = 0, j1 = 0, j2 = 0;
	float hi, lo, scale, x1, x2, xleft, xright, xscale;
	float y1, y2, ybottom, yscale, ytop;
	float max1, min1;
	int nx, ny;
	float fg = 0.0f, bg = 0.0f;
	int i = 0, j = 0;
	float r, g, b;
	float tr[6]; //tr[2] & tr[4] usually 0
	/*float **arr1 = (float **)malloc(sizeof(float *) * ndm);
	for (i = 0; i < ndm; i++) arr1[i] = (float *)malloc(sizeof(float) * nr);
	for (i = 0; i < ndm; i++) 
		for (j = 0; j < nr; j++) {
			arr1[i][j] = 0.1 + i * 1.5;
		}*/
	tr[0] = 10000.0;
	tr[1] = 10000.0;
	tr[2] = 0.0;
	tr[3] = 0.0;
	tr[4] = 0.0;
	tr[5] = 1.5;
	
	max1 = *(arr + 0);
	min1 = *(arr + 0);
	for (i = 1; i < ndm; i++) 
		for (j = 0; j < nr; j++) //*(arr + j + i * nr) = *(arr + j + i * nr) * 100;
		{
			max1 = max(*(arr + j + i * nr), max1);
			min1 = min(*(arr + j + i * nr), min1);
		}
	printf("%f %f\n", max1, min1);
	nx = nr;
	ny = ndm;
	hi = max1;
	lo = min1;
	if (cpgbeg(0, "/XS", 1, 1) != 1)
		return EXIT_FAILURE;
	
	//cpgwnad(0.0f, 1.0f, 0.0f, 1.0f);
	//cpgqwin(&x1, &x2, &y1, &y2);
	/*xscale = (x2 - x1)/nx;
	yscale = (y2 - y1)/ny;
	printf("%f %f %f %f\n", x1, x2, y1, y2);
	scale = (xscale < yscale) ? xscale : yscale;
	xleft = 0.5f * (x1 + x2 - nx * scale);
	xright = 0.5f * (x1 + x2 + nx * scale);
	ybottom = 0.5f * (y1 + y2 - ny * scale);
	ytop = 0.5f * (y1 + y2 + ny * scale);
	tr[0] = xleft - 0.5f * scale;
	tr[1] = scale;
	tr[2] = 0.0;
	tr[3] = ybottom - 0.5f * scale;
	tr[4] = 0.0;
	tr[5] = scale;*/
	//printf("%f %f %f %f %f %f\n", tr[0], tr[1], tr[2], tr[3], tr[4], tr[5]);
	/*cpgpage();
	cpgenv(0.0, nr, 0.0, ndm, 0, 0);
	cpgwnad(0.0f, 1.0f, 0.0f, 1.0f);
	cpgqwin(&x1, &x2, &y1, &y2);
	xscale = (x2 - x1)/nx;
	yscale = (y2 - y1)/ny;
	scale = (xscale < yscale) ? xscale : yscale;
	
	xleft = 0.5f * (x1 + x2 - nx * scale);
	xright = 0.5f * (x1 + x2 + nx * scale);
	ybottom = 0.5f * (y1 + y2 - ny * scale);
	ytop = 0.5f * (y1 + y2 + ny * scale);
	float tr[] = {xleft - 0.5f * scale, scale, 0.0f,
					ybottom - 0.5f * scale, 0.0f, scale};
	cpgimag(arr, nx, ny, 1, nx, 1, ny, hi, lo, tr);*/
	//cpgscir(1, 0);
	//cpgqcir(&cr1, &cr2);
	//printf("cr1: %d cr2: %d\n", cr1, cr2);
	//cpgpage();
	//cpgscr(0, 0.0, 0.3, 0.2);
	//cpgenv(0.0, nr, 0.0, ndm, 0, 0);
	//cpgsvp(0.05, 0.95, 0.05, 0.95);
	//cpgwnad(0.0, 1.0 , 0.0, 1.0);
	//cpgsci(1);
	//cpgsitf(0);
	//cpgmtxt("t", 1.0, 0.0, 0.0, "11");
	//cpgbox("bcnts", 0.0, 0, "bcnts", 0.0, 0);
	//cpgpixl(arr1[0], idim, jdim, i1, i2, j1, j2, 0.0, 1.0, 0.0, 1.0);
	cpgenv(1 - 0.5, nr+0.5, 1 - 0.5, ndm + 0.5, 0, 0);
	cpgqwin(&x1, &x2, &y1, &y2);
	xscale = (x2 - x1)/nx;
	yscale = (y2 - y1)/ny;
	printf("%f %f %f %f\n", x1, x2, y1, y2);
	scale = (xscale < yscale) ? xscale : yscale;
	xleft = 0.5f * (x1 + x2 - nx * scale);
	xright = 0.5f * (x1 + x2 + nx * scale);
	ybottom = 0.5f * (y1 + y2 - ny * scale);
	ytop = 0.5f * (y1 + y2 + ny * scale);
	tr[0] = xleft - 0.5f * scale;
	tr[1] = scale;
	tr[2] = 0.0;
	tr[3] = ybottom - 0.5f * scale;
	tr[4] = 0.0;
	tr[5] = scale;
	cpgscir(1, 50);
	cpgsitf(0);
	cpgimag(arr, nx, ny, 1, nx, 1, ny, hi, lo, tr);
	cpgbbuf();
	cpgend();
	return EXIT_SUCCESS;
}
/* ****************************************
 * Use pgplot to draw a histogram picture *
 * ****************************************/
int draw_histo_pic(struct pulse_num *pn)
{
	int i = 0;
	float datmax = 0.0f;
	for (i = 0; i < pn->pulse_count; i++) {
		if (datmax < *(pn->single_pulse_flux + i))
			datmax = *(pn->single_pulse_flux + i);
	}

	if (cpgbeg(0, "/xs", 1, 1) != 1)
		return EXIT_FAILURE;
	
	cpgenv(0.0, datmax + 1, 0.0, 100.0, 0, 1);
	cpghist(pn->pulse_count, pn->single_pulse_flux, 0.0, datmax + 1, 15, 1);
	
	cpgbbuf();
	cpgend();
	return EXIT_SUCCESS;
}
