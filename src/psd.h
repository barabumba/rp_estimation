/*
 * psd.h
 *
 *  Created on: May 2, 2019
 *      Author: mturbin
 */

#ifndef SRC_PSD_H_
#define SRC_PSD_H_

double psd_function(double x);
double psd_correction_log(double x, void *p);
void psd_fill(double *psd_sample, int N, int M);

#endif /* SRC_PSD_H_ */
