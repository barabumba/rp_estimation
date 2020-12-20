/*
 * processing.c
 *
 *  Created on: May 2, 2019
 *      Author: mturbin
 */
#include "defs.h"


double process_channel(double *c2, double *psd_ch_sample, double *correction_base, double q_ch, double m_ch)
{
	int64_t i;
	double sum=0;
	q_ch *= M/m_ch;
	for(i=0;i<N/2+1;i++)
	{
		sum += c2[i]*psd_ch_sample[i]/(q_ch*psd_ch_sample[i] + 1);
	}
	sum -= 0.5*c2[0]*psd_ch_sample[0]/(q_ch*psd_ch_sample[0] + 1) + 0.5*c2[N/2]*psd_ch_sample[N/2]/(q_ch*psd_ch_sample[N/2] + 1);
	sum *= q_ch/N;
	sum -= m_ch * (*correction_base) / 2;
	return sum;
}

