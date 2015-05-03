#ifndef _FOLD_H
#define _FOLD_H
#include "polyco.h"

typedef struct foldbuf {
	int nbin;
	int nchan;
	int npol;
	float *data;
	unsigned *count;
}foldbuf;

void malloc_foldbuf(foldbuf *f);

void free_foldbuf(foldbuf *f);

void clear_foldbuf(foldbuf *f);

size_t foldbuf_data_size(const foldbuf *f);
size_t foldbuf_count_size(const foldbuf *f);

int normalize_transpose_folds(float *out, const foldbuf *f);

typedef struct fold_args {
    polyco *pc;
    int imjd;
    double fmjd;
    char *data;
    int nsamp;
    double tsamp;
    int raw_signed;
    foldbuf *fb;
    float *scale;
    float *offset;	
}fold_args;

void *fold_8bit_power_thread(void *_args);

int fold_8bit_power(const polyco *pc, int imjd, double fmjd,
		const char *data, int nsamp, double tsamp, int raw_signed,
		foldbuf *f);
		
void *fold_16bit_power_thread(void *_args);

int fold_16bit_power(const polyco *pc, int imjd, double fmjd,
		const int16_t *data, int nsamp, double tsamp, int raw_signed,
		foldbuf *f);
		
int scale_offset_folds(foldbuf *f,
		const float *scale, const float *offset);

int accumulate_folds(foldbuf *ftot, const foldbuf *f);

#endif
