/*
 * processing.c
 *
 *  Created on: May 2, 2019
 *      Author: mturbin
 */
#include "defs.h"
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "processing.h"
#include "psd.h"

double process_channel(double *c2, double *psd_ch_sample, double (*correction_base)(), double q_ch, double m_ch)
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
	sum -= m_ch * (*correction_base)(q_ch) / 2;
	return sum;
}

double process_channel_gsl(double x, void *p)
{
	double output;
	struct my_f_params *params = (struct my_f_params *)p;

	double *psd = (double*) malloc((N/2+1)*sizeof(double));
	psd_fill(psd, N/2+1, x);

	output = process_channel(params->c2, psd, params->correction_base, params->q_ch, x);
	free(psd);
	return -output;
}

double find_min(gsl_function *F)
{
	int iter = 0, status;
	const int max_iter = 100;

	double m = M;
	double a = m/2, b = 2*m;

	double f = process_channel_gsl(a, F->params);
	double f_m = process_channel_gsl(m, F->params);
	if(f < f_m)
	{
		do {
			a /= 2;
			m /= 2;
			b /= 2;
			f = process_channel_gsl(a, F->params);
			f_m = process_channel_gsl(m, F->params);;
		} while(f < f_m);
	}
	else if((f = process_channel_gsl(b, F->params)) < f_m)
	{
		do {
			a *= 2;
			m *= 2;
			b *= 2;
			f = process_channel_gsl(b, F->params);
			f_m = process_channel_gsl(m, F->params);;
		} while(f < f_m);
	}

	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;

	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
//	printf("%f %f %f \n", process_channel_gsl(a, F->params), process_channel_gsl(b, F->params), process_channel_gsl(m, F->params));
	gsl_min_fminimizer_set (s, F, m, a, b);

//	printf ("using %s method\n",
//			gsl_min_fminimizer_name (s));
//
//	printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//		  "iter", "lower", "upper", "min",
//		  "err", "err(est)");
//
//	printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
//		  iter, a, b,
//		  m, b - a);

	do
	{
	  iter++;
	  status = gsl_min_fminimizer_iterate (s);

	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);

	  status
		= gsl_min_test_interval (a, b, 0.5, 0.0);

//	  if (status == GSL_SUCCESS)
//		printf ("Converged:\n");
//
//	  printf ("%5d [%.7f, %.7f] "
//			  "%.7f %.7f\n",
//			  iter, a, b,
//			  m, b - a);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_min_fminimizer_free (s);

	if (status != GSL_SUCCESS)
		printf ("Do not converged\n");
	return m;
}
