/*
 * sample_gen.h
 *
 *  Created on: May 1, 2019
 *      Author: mturbin
 */

#ifndef SRC_SAMPLE_GEN_H_
#define SRC_SAMPLE_GEN_H_

void SampleGauss(gsl_rng *gen, double *x, int length);
void GenSquaredCoefficients(gsl_rng *gen, double *psd_sample, double *c2);

#endif /* SRC_SAMPLE_GEN_H_ */
