/*
 * sample_gen.c
 *
 *  Created on: May 1, 2019
 *      Author: mturbin
 */
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include "defs.h"

void SampleGauss(gsl_rng *gen, double *x, int length)
{
	double f, x1, x2, r2;
	while(1)
	{
        do {
            x1 = 2.0*gsl_rng_uniform(gen) - 1.0;
            x2 = 2.0*gsl_rng_uniform(gen) - 1.0;
            r2 = x1*x1 + x2*x2;
        }
        while (r2 >= 1.0 || r2 == 0.0);
        f = sqrt(-2.0*log(r2)/r2);

        *(x++) = x1*f;
        if(--length==0)
        	return;
        *(x++) = x2*f;
        if(--length==0)
        	return;
	}
}

void GenSquaredCoefficients(gsl_rng *gen, double *psd_sample, double *c2)
{
	int i=0;
	double sqrt_psd_sample;
	double *x1, *x2;

	x1 = (double*) malloc(N*sizeof(double));
	x2 = (double*) malloc(N*sizeof(double));
	SampleGauss(gen, x1, N);
	gsl_fft_real_radix2_transform(x1, 1, N);
	SampleGauss(gen, x2, N);
	gsl_fft_real_radix2_transform(x2, 1, N);

	sqrt_psd_sample = sqrt(Q*psd_sample[0]);
	c2[0] = gsl_pow_2(x1[0]*sqrt_psd_sample+x2[0]);
	for(i = 1; i < N - i; i++)
	{
		sqrt_psd_sample = sqrt(Q*psd_sample[i]);
		c2[i] = gsl_pow_2(x1[i]*sqrt_psd_sample+x2[i]) + gsl_pow_2(x1[N-i]*sqrt_psd_sample+x2[N-i]);
	}
	if(i == N - i)
	{
		sqrt_psd_sample = sqrt(Q*psd_sample[i]);
		c2[i] = gsl_pow_2(x1[i]*sqrt_psd_sample+x2[i]);
	}
	free(x1);
	free(x2);
}
