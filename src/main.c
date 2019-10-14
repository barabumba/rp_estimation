/*
 * main.c
 *
 *  Created on: Apr 29, 2019
 *      Author: mturbin
 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>
//#include <omp.h>
#include "sample_gen.h"
#include "psd.h"
#include "defs.h"
#include "processing.h"
#include "misc.h"
#include "gsl/gsl_statistics_double.h"

double q_search_values[Q_SRCH_RANGE_LEN];
double bw_search_values[M_SRCH_RANGE_LEN];
double psd_sample_true[N/2+1];


double calculate_correction(double q)
{
	gsl_integration_workspace *giw = gsl_integration_workspace_alloc(100);
	double result, abserror;
	gsl_function F;
	F.function = &psd_correction_log;
	F.params = &q;
	gsl_integration_qagi(&F, 1e-6, 1e-6, 100, giw, &result, &abserror);
	gsl_integration_workspace_free(giw);
	return result;
}





int main(void)
{
	psd_fill(psd_sample_true, N/2+1, M);

	for(int i=0;i<Q_SRCH_RANGE_LEN;i++)
		q_search_values[i] = (Q_SRCH_RANGE_LEN>1) ? Q/2.+i*(1.*Q/Q_SRCH_RANGE_LEN) : Q;

	for(int i=0;i<M_SRCH_RANGE_LEN;i++)
		bw_search_values[i] = (M_SRCH_RANGE_LEN>1) ? M/2+i*(M/M_SRCH_RANGE_LEN) : M;
	assert(bw_search_values[M_SRCH_RANGE_LEN/2] == M);

	for(int q_indx=0; q_indx<Q_SRCH_RANGE_LEN; q_indx++)
	{
		double correction_base = calculate_correction(q_search_values[q_indx]);
//#pragma omp parallel
		{
			double *c2 = (double*) malloc((N/2+1)*sizeof(double));
			double *statistic = (double*) malloc((NUMBER_OF_TESTS)*sizeof(double));
			gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);

			gsl_function F;
			F.function = &process_channel_gsl;
//#pragma omp critical
		{
				gsl_rng_set(gen, GetRandomSeed_64b());
		}
//#pragma omp for
			for(int i=0; i<NUMBER_OF_TESTS; i++)
			{
				double bw;
				GenSquaredCoefficients(gen, psd_sample_true, c2);
				struct my_f_params parameters = {c2, &correction_base, q_search_values[q_indx]};
				F.params = &parameters;

				bw = find_min(&F);

			}
		}
	}
	return 0;
}
