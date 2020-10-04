/*
 * main.c
 *
 *  Created on: Apr 29, 2019
 *      Author: mturbin
 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <omp.h>
#include "sample_gen.h"
#include "psd.h"
#include "defs.h"
#include "processing.h"
#include "misc.h"

double psd_array[M_SRCH_RANGE_LEN][N/2+1];
double correction_base[Q_SRCH_RANGE_LEN];
int64_t bw_search_values[M_SRCH_RANGE_LEN];
double q_search_values[Q_SRCH_RANGE_LEN];
double *psd_sample_true;

int64_t SaveGSLMatrtix(gsl_matrix *m, char *name);
int64_t SaveGSLMatrtixInt(gsl_matrix_int *m, char *name);

int main(void)
{
	gsl_matrix_int *mutual_results = gsl_matrix_int_calloc((size_t) Q_SRCH_RANGE_LEN, (size_t) M_SRCH_RANGE_LEN);

	for(int i=0;i<M_SRCH_RANGE_LEN;i++)
	{
		bw_search_values[i] = (M_SRCH_RANGE_LEN>1) ? M_SRCH_MIN+i*((M_SRCH_MAX-M_SRCH_MIN)/M_SRCH_RANGE_LEN) : M;
		printf("%d, ", M/2+i*(M/M_SRCH_RANGE_LEN));
		psd_fill(psd_array[i], N/2+1, bw_search_values[i]);
	}
	printf("\n");
	assert(bw_search_values[M_SRCH_RANGE_LEN/2] == M);
	psd_sample_true = psd_array[M_SRCH_RANGE_LEN/2];

	gsl_integration_workspace *giw = gsl_integration_workspace_alloc(100);
	for(int i=0;i<Q_SRCH_RANGE_LEN;i++)
	{
		double result, abserror;
		gsl_function F;
		q_search_values[i] = (Q_SRCH_RANGE_LEN>1) ? Q_SRCH_MIN+i*((Q_SRCH_MAX-Q_SRCH_MIN)/Q_SRCH_RANGE_LEN) : Q;
		F.function = &psd_correction_log;
		F.params = &(q_search_values[i]);
		gsl_integration_qagi(&F, 1e-6, 1e-6, 100, giw, &result, &abserror);
		correction_base[i] = result;
	}
	gsl_integration_workspace_free(giw);
	assert(q_search_values[Q_SRCH_RANGE_LEN/2] == Q);

	//omp_set_num_threads(4);
#pragma omp parallel
{
	double *c2 = (double*) malloc((N/2+1)*sizeof(double));

	gsl_matrix *search_field = gsl_matrix_alloc((size_t) Q_SRCH_RANGE_LEN, (size_t) M_SRCH_RANGE_LEN);
	gsl_matrix_int *mutual_results_aux = gsl_matrix_int_calloc((size_t) Q_SRCH_RANGE_LEN, (size_t) M_SRCH_RANGE_LEN);

	gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);

#pragma omp critical
	{
	gsl_rng_set(gen, GetRandomSeed_64b());
	}

#pragma omp for
	for(int i=0;i<NUMBER_OF_TESTS;i++)
	{
		size_t q_indx, m_indx;

		GenSquaredCoefficients(gen, psd_sample_true, c2);
		for(q_indx=0;q_indx<Q_SRCH_RANGE_LEN;q_indx++)
		{
			for(m_indx=0;m_indx<M_SRCH_RANGE_LEN;m_indx++)
			{
				double channel_output = process_channel(c2, psd_array[m_indx], &correction_base[q_indx], q_search_values[q_indx], bw_search_values[m_indx]);
				gsl_matrix_set(search_field, q_indx, m_indx, channel_output);
			}
		}
		if(i%10000==0)
			printf("%d %d \n", omp_get_thread_num(), i);
		gsl_matrix_max_index(search_field, &q_indx, &m_indx);
		(*(gsl_matrix_int_ptr(mutual_results_aux, q_indx, m_indx)))++;
	}
#pragma omp critical
	{
		gsl_matrix_int_add(mutual_results, mutual_results_aux);
	}
	free(c2);
	gsl_matrix_free(search_field);
	gsl_matrix_int_free(mutual_results_aux);
	gsl_rng_free(gen);
}

	for(int i=0; i<Q_SRCH_RANGE_LEN;i++)
	{
		for(int j=0; j<M_SRCH_RANGE_LEN;j++)
		{
			printf("%03d ", gsl_matrix_int_get(mutual_results, i, j));
		}
		printf("\n");
	}
	SaveGSLMatrtixInt(mutual_results, "result_matrix.dat");
	gsl_matrix_int_free(mutual_results);
	return 0;
}

int64_t SaveGSLMatrtix(gsl_matrix *m, char *name)
{
	FILE *fp;
	fp = fopen(name, "w");

    if(fprintf(fp, "%f;%d:%f:%f\n", Q, Q_SRCH_RANGE_LEN, q_search_values[0], q_search_values[m->size1-1]) < 0)
    	return -1;

    if(fprintf(fp, "%d;%d:%ld:%ld\n", M, M_SRCH_RANGE_LEN, bw_search_values[0], bw_search_values[m->size2-1]) < 0)
    	return -1;

    for (size_t i = 0; i < m->size1; i++) {
    	for (size_t j = 0; j < m->size2; j++) {
			if ((fprintf(fp, "%f ", gsl_matrix_get(m, i, j))) < 0)
				return -1;
        }
		if ((fprintf(fp, "\n")) < 0)
			return -1;
    }
    fclose(fp);
    return 0;
}

int64_t SaveGSLMatrtixInt(gsl_matrix_int *m, char *name)
{
	FILE *fp;
	fp = fopen(name, "w");

    if(fprintf(fp, "%f;%d:%f:%f\n", Q, Q_SRCH_RANGE_LEN, q_search_values[0], q_search_values[m->size1-1]) < 0)
    	return -1;

    if(fprintf(fp, "%d;%d:%ld:%ld\n", M, M_SRCH_RANGE_LEN, bw_search_values[0], bw_search_values[m->size2-1]) < 0)
    	return -1;

    for (size_t i = 0; i < m->size1; i++) {
    	for (size_t j = 0; j < m->size2; j++) {
			if ((fprintf(fp, "%d ", gsl_matrix_int_get(m, i, j))) < 0)
				return -1;
        }
		if ((fprintf(fp, "\n")) < 0)
			return -1;
    }
    fclose(fp);
    return 0;
}
