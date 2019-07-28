/*
 * psd.c
 *
 *  Created on: Apr 29, 2019
 *      Author: mturbin
 */
#include <math.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double psd_function(double x)
{
//	return 1./(1.+pow(M_PI*x, 2));
	return exp(-M_PI*pow(x,2));
}

double psd_correction_log(double x, void *p)
{
	double q = *((double *)p);
	return log(1.+q*psd_function(x));
}

void psd_fill(double *psd_sample, int N, int M)
{
	/*
	 * It is assumed that observation time T=1, so sampling rate=1/N and, hence, frequency step=1.
	 */
	int i;
	for(i=0;i<N;i++)
		psd_sample[i] = psd_function(((double)i)/M);
}
