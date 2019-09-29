/*
 * misc.c
 *
 *  Created on: Sep 29, 2019
 *      Author: mturbin
 */
#include "defs.h"


uint64_t GetRandomSeed_64b(void)
{
	uint64_t data;
	FILE *fp;
	fp = fopen("/dev/urandom", "r");
	fread(&data, sizeof(data), 1, fp);
	fclose(fp);
	return data;
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
