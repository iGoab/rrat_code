#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
#include <math.h>

#define CVALLEN 10
#define SIGNUM 90000
#define DMNUM 10
#define RMSBLOCK 11

struct Fitsinfo {
	char pol_type[CVALLEN];
	int nCHAN;
	int nPOL;
	int nSBLK;
	int nROW;
	float tBIN;
	float dmArr[DMNUM];
	int lowF;
	int highF;
};

struct FreqData {
	float freq;
	int* data;
	float dm;
};

struct DMData {
	float dm[10];
	int value[10];
	int index;
	int dmNum;
};

struct BigSignal {
	int value[SIGNUM];
	int index[SIGNUM];
};

void printerror(int status);

void initDmArr(float dmArr[],float realDm);

void initFitsInfo(fitsfile* fptr,struct Fitsinfo* fi,float realDm,int* status);

void initTotalData(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** ptr,int* status);
void initSuitData(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** ptr,int* status);

void initTotalI(fitsfile* fptr,struct Fitsinfo* fi, unsigned char** iptr,int* status);
void initSuitI(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** iptr,int* status);

void initTotalAABB(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** abptr,int* status);
void initSuitAABB(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** abptr,int* status);

void initFreqData(fitsfile* fptr,struct FreqData* fd,unsigned char** ptr,int freqIndex,long samNum,int* status);

int calcZeroDm(unsigned char** ptr,int chanNum,int samIndex);

void initFreqDataWithZeroDm(fitsfile* fptr,struct FreqData* fd,unsigned char** ptr,int freqIndex,long samNum,int chanNum,int* status);

int getMoveDist(float tbin,float dm,float freqO,float freqM);

void addFreqData(struct FreqData* fd,struct FreqData* fdm,int dist,long samNum);

void addAllFreqData(fitsfile* fptr,struct FreqData* fd,struct FreqData* fdm,unsigned char** ptr,float dm,struct Fitsinfo* fi,int* status);

void addSuitFreqData(fitsfile* fptr,struct FreqData* fd,struct FreqData* fdm,unsigned char** ptr,float dm,struct Fitsinfo* fi,int* status);

void addSuitFreqDataWithZeroDm(fitsfile* fptr,struct FreqData* fd,struct FreqData* fdm,unsigned char** ptr,float dm,struct Fitsinfo* fi,int* status);

void writeIntData(char* filename,float dm,int* data,int dataNum);

void writeFloatData(char* filename,float dm,float* data,int dataNum);

void fold(int foldData[],int* data,float tbin,double psrFreq,long samNum);

int calcAveArr(float aveData[],int* data,int totalNum,int aveNum);

void getDataFromSamps(int resSamp[],int* data, int num, int startIndex);

int getIndexofMaxfromArr(int* arr,int num,int* maxVal);

int getRms(int* val, int size);

int getBigSignalNumFromData(int* data,int size,int rms,float rate);

void writeBigSignal(struct BigSignal bs,struct FreqData** fd,int signum);

void compTwoFreq(fitsfile* fptr,int freqlow,int freqhigh,int dist,struct Fitsinfo* fi,unsigned char** ptr,int* status);

void writeHistogramData(struct FreqData** fd,struct Fitsinfo* fi);

int main(int argc, char *argv[])
{
	fitsfile* fptr;
	//char filename[] = "guppi_56590_B1133+16_0003_0001.fits";
	char filename[] = "~/pulsar_data/B0329.fits";
	int status = 0;
	float realDm = 26.8;
	double psrFreq = 1.399541539;
	int i,j; //for loop
	int sampNum;


	fits_open_file(&fptr,filename,READONLY,&status);
	
	//初始化fits文件信息
	struct Fitsinfo* fi = (struct Fitsinfo*)malloc(sizeof(struct Fitsinfo));
	initFitsInfo(fptr,fi,realDm,&status);
	sampNum = fi->nSBLK*fi->nROW;
	//初始化保存所有流量数据的二维数组
	/*
	unsigned char** ptr = (unsigned char**)malloc(fi->nROW*fi->nSBLK*sizeof(unsigned char*));
	for (i=0; i<fi->nROW*fi->nSBLK; i++) {
		ptr[i] = (unsigned char*)malloc(fi->nCHAN*sizeof(unsigned char));
	}
	initTotalData(fptr,fi,ptr,&status);
	*/

	//初始化合适的频道
	unsigned char** ptr = (unsigned char**)malloc(fi->nROW*fi->nSBLK*sizeof(unsigned char*));
	for (i=0; i<fi->nROW*fi->nSBLK; i++) {
		ptr[i] = (unsigned char*)malloc((fi->highF-fi->lowF+1)*sizeof(unsigned char));
	}
	initSuitData(fptr,fi,ptr,&status);
	//给需要处理的某两个通道分配空间	
    	struct FreqData* fdm = (struct FreqData*)malloc(sizeof(struct FreqData));
	fdm->data = (int*)malloc(sizeof(int)*fi->nSBLK*fi->nROW);
	struct FreqData** fd = (struct FreqData**)malloc(DMNUM*sizeof(struct FreqData*));	
	for (i=0; i<DMNUM; i++) {
		fd[i] = (struct FreqData*)malloc(sizeof(struct FreqData));
		fd[i]->data = (int*)malloc(sizeof(int)*fi->nSBLK*fi->nROW);
		//addSuitFreqData(fptr,fd[i],fdm,ptr,fi->dmArr[i],fi,&status);
		addSuitFreqDataWithZeroDm(fptr,fd[i],fdm,ptr,fi->dmArr[i],fi,&status);
	}
	for (i=0; i<fi->nROW*fi->nSBLK; i++) {
		free(ptr[i]);
	}
	free(ptr);

	/* draw fd[DMNUM/2]->data */
	/* draw fd[DMNUM/2]->data */
	if (cpgbeg(0, "/xs", 1, 1) != 1) {
		return EXIT_FAILURE;
	}
	float ns = fi->nROW * fi->nSBLK;
	float *x_coor, *data_adm;
	x_coor = (float *)malloc(sizeof(float) * fi->nROW * fi->nSBLK);
	data_adm = (float *)malloc(sizeof(float) * ns);
	for (i = 0; i < fi->nROW * fi->nSBLK; i++) {
		*(x_coor + i) = i;
		*(data_adm  + i) = (float)fd[DMNUM/2]->data[i];
	}
	
	cpgpage();
	cpgenv(0.0, fi->nROW * fi->nSBLK, -5000.0f, 5000.0f, 0, 1);
	
	cpgline(ns, x_coor, data_adm);
	cpgbbuf();
	cpgend();
	return EXIT_SUCCESS;
	
	float aveData[20000] = {0};
	writeIntData("26.8.dat",26.8,fd[DMNUM/2]->data,fi->nSBLK*fi->nROW);
	calcAveArr(aveData,fd[DMNUM/2]->data,2000000,100);
	writeFloatData("ave26.8.dat",26.8,aveData,20000);
	int foldData[20000] = {0};
	fold(foldData,fd[DMNUM/2]->data,fi->tBIN,psrFreq,fi->nSBLK*fi->nROW);

	writeHistogramData(fd,fi);
	
	

	fits_close_file(fptr,&status);

}

//打印错误信息
void printerror(int status) {
    if (status) {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* stop the program, returning error status */
    }
    return;
}

/*
初始化实验的dm数组
*/
void initDmArr(float dmArr[],float realDm) {
	int i;
	int delta = DMNUM/2;
	for (i=0; i<DMNUM; i++) {
		dmArr[i] = realDm-delta+i;
	}
	return;
}

/*
初始化fits信息的结构体
*/
void initFitsInfo(fitsfile* fptr,struct Fitsinfo* fi,float realDm,int* status) {
	int value;
	float tbin;
	char strVal[CVALLEN];
	//fits_read_key(fptr,TSTRING,"OBS_MODE",strVal,NULL,status);
	//strcpy(fi->obs_Mode,strVal);
	fits_movabs_hdu(fptr, 2, NULL, status);
	fits_read_key(fptr,TSTRING,"POL_TYPE",strVal,NULL,status);
	strcpy(fi->pol_type,strVal);
	fits_read_key(fptr,TINT,"NCHAN",&value,NULL,status);
	fi->nCHAN = value;
	fits_read_key(fptr,TINT,"NPOL",&value,NULL,status);
	fi->nPOL = value;
	fits_read_key(fptr,TINT,"NSBLK",&value,NULL,status);
	fi->nSBLK = value;
	fits_read_key(fptr,TINT,"NAXIS2",&value,NULL,status);
	fi->nROW = value;
	fits_read_key(fptr,TFLOAT,"TBIN",&tbin,NULL,status);
	fi->tBIN = tbin;
	initDmArr(fi->dmArr,realDm);
	fi->lowF = 280;
	fi->highF = 343;
	return;
}
/*
判断pol的模式，选取合适的函数来初始化所有强度值
*/
void initTotalData(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** ptr,int* status) {
	if (!strcmp(fi->pol_type,"AABBCRCI")) {
		initTotalAABB(fptr,fi,ptr,status);
	} else if(!strcmp(fi->pol_type,"IQUV")) {
		initTotalI(fptr,fi,ptr,status);
	}
}

void initSuitData(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** ptr,int* status) {
	if (!strcmp(fi->pol_type,"AABBCRCI")) {
		initSuitAABB(fptr,fi,ptr,status);
	} else if(!strcmp(fi->pol_type,"IQUV")) {
		initSuitI(fptr,fi,ptr,status);
	}
}
/*
如果是IQUV模式，则按此函数来初始化一个二维unsigned char型数组
*/
void initTotalI(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** iptr,int* status) {
	int n = fi->nROW; int k = fi->nSBLK; int m = fi->nCHAN;
	int step = fi->nCHAN*fi->nPOL; int start;
	unsigned char fluxI[m];
	int dataCol; int nulval = 0; int anynull; int i,j,p;
	fits_get_colnum(fptr,CASEINSEN,"DATA",&dataCol,status);
	for (i=0; i<n; i++) {
		for (j=0; j<k; j++) {
			start = 1+j*step;
			fits_read_col(fptr,TBYTE,dataCol,i+1,start,m,&nulval,fluxI,&anynull,status);
			for (p=0; p<m; p++) {
				iptr[j+i*k][p] = fluxI[p];
			}
		}
		printf("All Freq I init:Totalrow %d,row completed %d\n",n,i);
	}
}

void initSuitI(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** iptr,int* status) {
	int n = fi->nROW; int k = fi->nSBLK; int lf = fi->lowF; int hf = fi->highF;
	int chanNum = hf-lf+1;
	int step = fi->nCHAN*fi->nPOL; int start;
	unsigned char fluxI[chanNum];
	int dataCol; int nulval = 0; int anynull; int i,j,p;
	fits_get_colnum(fptr,CASEINSEN,"DATA",&dataCol,status);
	for (i=0; i<n; i++) {
		for (j=0; j<k; j++) {
			start = lf+1+j*step;
			fits_read_col(fptr,TBYTE,dataCol,i+1,start,chanNum,&nulval,fluxI,&anynull,status);
			for (p=0; p<chanNum; p++) {
				iptr[j+i*k][p] = fluxI[p];
			}
		}
		printf("Suit Freq I init:Totalrow %d,row completed %d\n",n,i);
	}
}

/*
初始化整个aa+bb序列，为了省内存，用unsigned char类型存储，并除2使之不超过范围
*/
void initTotalAABB(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** abptr,int* status) {
	int n = fi->nROW;int k = fi->nSBLK; int m = fi->nCHAN;
	int step = fi->nCHAN*fi->nPOL; int start;
	unsigned char aa[m]; unsigned char bb[m];
	int dataCol; int nulval = 0; int anynull; int i,j,p;
	fits_get_colnum(fptr,CASEINSEN,"DATA",&dataCol,status);
	for (i=0; i<n; i++) {
		for (j=0; j<k; j++) {
			start = 1+j*step;
			fits_read_col(fptr,TBYTE,dataCol,i+1,start,m,&nulval,aa,&anynull,status);
			fits_read_col(fptr,TBYTE,dataCol,i+1,start+m,m,&nulval,bb,&anynull,status);
			for (p=0; p<m; p++) {
				abptr[j+i*k][p] = (aa[p]+bb[p])/2;
			}
		}
		printf("All Freq AABB init:Totalrow %d,row completed %d\n",n,i);
	}
	return;
}

void initSuitAABB(fitsfile* fptr,struct Fitsinfo* fi,unsigned char** abptr,int* status) {
	int n = fi->nROW;int k = fi->nSBLK; int lf = fi->lowF; int hf = fi->highF; int m = fi->nCHAN; 
	int chanNum = hf-lf+1;
	int step = chanNum*fi->nPOL; int start;
	unsigned char aa[chanNum]; unsigned char bb[chanNum];
	int dataCol; int nulval = 0; int anynull; int i,j,p;
	fits_get_colnum(fptr,CASEINSEN,"DATA",&dataCol,status);
	for (i=0; i<n; i++) {
		for (j=0; j<k; j++) {
			start = lf+1+j*step;
			fits_read_col(fptr,TBYTE,dataCol,i+1,start,chanNum,&nulval,aa,&anynull,status);
			fits_read_col(fptr,TBYTE,dataCol,i+1,start+m,chanNum,&nulval,bb,&anynull,status);
			for (p=0; p<chanNum; p++) {
				abptr[j+i*k][p] = (aa[p]+bb[p])/2;
			}
		}
		printf("Suit Freq AABB init:Totalrow %d,row completed %d\n",n,i);
	}
	return;
}

/**
@param freqIndex 从0开始
*/
void initFreqData(fitsfile* fptr,struct FreqData* fd,unsigned char** ptr,int freqIndex,long samNum,int* status) {
	int freqCol; int nulval = 0; int anynull; float tmpFreq; int i;
	fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&freqCol,status);
	fits_read_col(fptr,TFLOAT,freqCol,1,freqIndex+1,1,&nulval,&tmpFreq,&anynull,status);
	fd->freq = tmpFreq;
	for (i=0; i<samNum; i++) {
		fd->data[i] = ptr[i][freqIndex];
	}
	return;
}

int calcZeroDm(unsigned char** ptr,int chanNum,int samIndex) {
	int zeroDmVal = 0;
	int sum = 0;
	int i;
	for (i=0; i<chanNum; i++) {
		sum += ptr[samIndex][i];
	}
	zeroDmVal = sum/chanNum;
	return zeroDmVal;
}
	

void initFreqDataWithZeroDm(fitsfile* fptr,struct FreqData* fd,unsigned char** ptr,int freqIndex,long samNum,int chanNum,int* status) {
	int freqCol; int nulval = 0; int anynull; float tmpFreq; int i;
	fits_get_colnum(fptr,CASEINSEN,"DAT_FREQ",&freqCol,status);
	fits_read_col(fptr,TFLOAT,freqCol,1,freqIndex+1,1,&nulval,&tmpFreq,&anynull,status);
	fd->freq = tmpFreq;
	for (i=0; i<samNum; i++) {
		fd->data[i] = ptr[i][freqIndex]-calcZeroDm(ptr,chanNum,i);
	}
	return;
}

/*
根据色散公式，求各个频道分别偏移点数
@param freqO为目标频道
@param freqM为需要移动的频道
@return retDist为偏移点数
*/

int getMoveDist(float tbin,float dm,float freqO,float freqM) {
	int retDist;
	if (freqO <= freqM) {
		retDist = dm*(1/(freqO*freqO)-1/(freqM*freqM))*4.15e+3/tbin;
	} else {
		retDist = dm*(1/(freqM*freqM)-1/(freqO*freqO))*4.15e+3/tbin;
	}
	return retDist;
}
/*
根据偏移点数将某两个频道的采样点叠加
*/
void addFreqData(struct FreqData* fd, struct FreqData* fdm,int dist,long samNum) {
	int i;
	if (fd->freq > fdm->freq) {
		for (i=0; i+dist<samNum; i++) {
			fd->data[i] = fd->data[i]+fdm->data[i+dist];
		}
	} else {
		for (i=dist; i<samNum; i++) {
			fd->data[i] = fd->data[i]+fdm->data[i-dist];
		}
	}
	return;
}

/*
将所有的频道叠加，并记录dm在最终的FreqData中
*/
void addAllFreqData(fitsfile* fptr,struct FreqData* fd,struct FreqData* fdm,unsigned char** ptr,float dm,struct Fitsinfo* fi,int* status) {
	int i,dist;
	int samNum = fi->nROW*fi->nSBLK;
     	initFreqData(fptr,fd,ptr,0,samNum,status);
	fd->dm = dm;
     	for (i=1; i<fi->nCHAN; i++) {
        	initFreqData(fptr,fdm,ptr,i,samNum,status);
        	dist = getMoveDist(fi->tBIN,dm,fd->freq,fdm->freq);
        	addFreqData(fd,fdm,dist,samNum);
        	printf("Total:%d chan,Now:%d,dist==%d,dm==%.1f\n",fi->nCHAN,i,dist,fd->dm);
        }
}

void addSuitFreqData(fitsfile* fptr,struct FreqData* fd,struct FreqData* fdm,unsigned char** ptr,float dm,struct Fitsinfo* fi,int* status) {
	int i,dist,ave;
	int chanNum = fi->highF-fi->lowF+1;
	int samNum = fi->nROW*fi->nSBLK; 
        initFreqData(fptr,fd,ptr,0,samNum,status);
	fd->dm = dm;
	long long sum =0;
        for (i=1; i<chanNum; i++) {
        	initFreqData(fptr,fdm,ptr,i,samNum,status);
        	dist = getMoveDist(fi->tBIN,dm,fd->freq,fdm->freq);
        	addFreqData(fd,fdm,dist,samNum);
        	printf("Suit:%d chan,Now:%d,dist==%d,dm==%.1f\n",chanNum,i,dist,fd->dm);
        }
	for (i=0; i<samNum; i++) {
		sum += fd->data[i];
	}
	ave = sum/samNum;	
	for (i=0; i<samNum; i++) {
		fd->data[i] = fd->data[i]-ave;
	}
}

void addSuitFreqDataWithZeroDm(fitsfile* fptr,struct FreqData* fd,struct FreqData* fdm,unsigned char** ptr,float dm,struct Fitsinfo* fi,int* status) {
	int i,dist,ave;
	int chanNum = fi->highF-fi->lowF+1;
	int samNum = fi->nROW*fi->nSBLK; 
        initFreqDataWithZeroDm(fptr,fd,ptr,0,samNum,chanNum,status);
	fd->dm = dm;
	long long sum =0;
        for (i=1; i<chanNum; i++) {
        	initFreqDataWithZeroDm(fptr,fdm,ptr,i,samNum,chanNum,status);
        	dist = getMoveDist(fi->tBIN,dm,fd->freq,fdm->freq);
        	addFreqData(fd,fdm,dist,samNum);
        	printf("Suit WithZeroDM:%d chan,Now:%d,dist==%d,dm==%.1f\n",chanNum,i,dist,fd->dm);
        }
	for (i=0; i<samNum; i++) {
		sum += fd->data[i];
	}
	ave = sum/samNum;	
	for (i=0; i<samNum; i++) {
		fd->data[i] = fd->data[i]-ave;
	}
}

/*
写入数据
*/
void writeIntData(char* filename,float dm,int* data,int dataNum) {
	FILE* fp;
	fp = fopen(filename,"w");
	int i;
	for (i=0; i<dataNum; i++) {
		fprintf(fp,"%d\n",data[i]);
		//printf("written in %d\n",data[i]);
	}
	fclose(fp);
}

void writeFloatData(char* filename,float dm,float* data,int dataNum) {
	FILE* fp;
	fp = fopen(filename,"w");
	int i;
	for (i=0; i<dataNum; i++) {
		fprintf(fp,"%f\n",data[i]);
	}
	fclose(fp);
} 

/*
将最终的消色散的数据折叠（还是有问题啊，看不到轮廓）
*/
void fold(int foldData[],int* data,float tbin,double psrFreq,long samNum) {
	int foldProNum = 1/tbin/psrFreq;
	//int foldData[foldProNum] = {0};
	int i,j;
	for (i=0,j=0; j<samNum; i++,j++) {
		foldData[i] += data[j];
		if (i == foldProNum-1) {
			i = 0;
		}
	}
	writeIntData("fold26.8.dat",26.8,foldData,foldProNum);
}		

/*
计算平滑数据
@param aveData 储存平滑数据的数组
@param data 需要平滑的数组
@param totalNum 待平滑的数组元素个数
@param  avenum  以 avenum 来做平滑
@return resNum  平滑后的数组元素个数（画图用）
*/
int calcAveArr(float aveData[],int* data,int totalNum,int aveNum) {
	float tmp;
	int i,j;
	int resNum = totalNum/aveNum;
	for (i=0,j=0; i<resNum; ) {
		tmp += data[j];
		j++;
		if (j%aveNum == 0) {
			aveData[i] = tmp/aveNum;
			tmp = 0;
			i++;
		}
	}
	return resNum;
}
/*
从各个频段消色散叠加后的数据中取走一部分，因为原始数据前一部分不合理
@param resSamp 取走的数据
@param data 原始数据
@param num 取出的数据个数
@param startIndex 从 startIndex 开始取  一般是最大偏移点数
*/
void getDataFromSamps(int resSamp[],int* data, int num, int startIndex) {
	int i;
	int index = startIndex;
	for (i=0; i<num; i++,index++) {
		resSamp[i] = data[index];
	}
}


int getIndexofMaxfromArr(int* arr,int num,int* maxVal) {
	int max = arr[0];
	int maxIndex;
	int i;
	for (i=0; i<num; i++) {
		if (arr[i]>max) {
			max = arr[i];
			maxIndex = i;
		}
	}
	*maxVal = max;
	return maxIndex;
}


/*
返回rms（均方根），用整型取近似，方便计算
*/
int getRms(int* val,int size) {
	double tmp = 0;
	int ret;
	int i;
	for (i=0; i<size; i++) {
		tmp += val[i]*val[i]/size;
	}
	ret = (int)sqrt(tmp);
	return ret;
}	
	

/*
@param data 频道强度序列指针
@param size   采样点总个数
@return j  大数据总个数
*/	
int getBigSignalNumFromData(int* data,int size,int rms,float rate) { 
	int i,j;
	for (i=0,j=0; i<size; i++) {
		if(data[i]>rate*rms) {
			j++;
		}
	}
	return j;
}


void writeBigSignal(struct BigSignal bs,struct FreqData** fd,int signum) {
	FILE* fp;
	fp = fopen("BigSignal.dat","w");
	int i,j;
	fprintf(fp,"0   ");
	for (j=0; j<DMNUM; j++) {
		fprintf(fp,"%f  ",fd[j]->dm);
	}
	fprintf(fp,"\n");
	for (i=0; i<signum; i++) {
		fprintf(fp,"%d  ",bs.index[i]);
		for (j=0; j<DMNUM; j++) {
			fprintf(fp,"%d  ",fd[j]->data[bs.index[i]]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void compTwoFreq(fitsfile* fptr,int freqlow,int freqhigh,int dist,struct Fitsinfo* fi,unsigned char** ptr,int* status) {
	int i;
    	struct FreqData* fd = (struct FreqData*)malloc(sizeof(struct FreqData));
	fd->data = (int*)malloc(sizeof(int)*fi->nSBLK*fi->nROW);	
	struct FreqData* fdtmp = (struct FreqData*)malloc(sizeof(struct FreqData));
	fdtmp->data = (int*)malloc(sizeof(int)*fi->nSBLK*fi->nROW);

	initFreqData(fptr,fd,ptr,freqlow,fi->nSBLK*fi->nROW,status);
	initFreqData(fptr,fdtmp,ptr,freqhigh,fi->nSBLK*fi->nROW,status);

	addFreqData(fd,fdtmp,dist,fi->nSBLK*fi->nROW);
	FILE* fp;
	fp = fopen("addtest.dat","w");
	fprintf(fp,"%5d%5d%5d\n",freqlow,freqhigh,0);
	for (i=0; i<20000; i++) {
		fprintf(fp,"%5d%5d%5d\n",ptr[i][freqlow],ptr[i][freqhigh],fd->data[i]);
	}
	fclose(fp);
}

void setArrToZero(int count[],int size) {
	int i;
	for (i=0; i<size; i++) {
		count[i] = 0;
	}
	return;
}

void writeHistogramData(struct FreqData** fd,struct Fitsinfo* fi) {
	int count[RMSBLOCK] = {0};
	int beyond[RMSBLOCK] = {0};
	int i,j,k;
	int tmp;
	long size = fi->nROW*fi->nSBLK;
	int rms;
	FILE* fp;
	fp = fopen("histoGram.dat","w");
	for (i=0; i<DMNUM; i++) {
		fprintf(fp,"dm %.1f : \n",fd[i]->dm);
		rms = getRms(fd[i]->data,size);
		fprintf(fp,"rms = %d\n",rms);
		for (j=0; j<RMSBLOCK-1; j++) {
			beyond[j] = getBigSignalNumFromData(fd[i]->data,size,rms,0.5+0.1*j);
		}
		for (k=0; k<RMSBLOCK; k++) {
			if (0 == k) {
				count[k] = size-beyond[k];
				fprintf(fp,"0-0.5:%d   ",count[k]);
			} else {
				count[k] = beyond[k-1]-beyond[k];
				fprintf(fp,"%.1f~:%d   ",0.5+0.1*(k-1),count[k]);
			}
		}
		fprintf(fp,"\n-------------------\n");
	}
}







