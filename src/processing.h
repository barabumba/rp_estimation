/*
 * processing.h
 *
 *  Created on: May 4, 2019
 *      Author: mturbin
 */

#ifndef SRC_PROCESSING_H_
#define SRC_PROCESSING_H_

struct my_f_params { double *c2; double (*correction_base)(); double q_ch;};

double process_channel(double *c2, double *psd_ch_sample, double (*correction_base)(), double q_ch, double m_ch);
double process_channel_gsl(double x, void *p);
double find_min(gsl_function *F);

#endif /* SRC_PROCESSING_H_ */
