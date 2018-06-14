#ifndef _DZETA_FUNCTION_UTILS_H
#define _DZETA_FUNCTION_UTILS_H


#define DZETA_EPSILON 1.e-14

inline double _DSQR( double const x[3] ) {
  return ( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
}

inline double _DSCP( double x[3], double y[3] ) {
  return( x[0] * y[0] + x[1] * y[1] + x[2] * y[2] );
}

#define _PLM_L_MAX 6

extern double plm_norm[(_PLM_L_MAX+1)*(_PLM_L_MAX+1)];

int dinit_plm_norm (void);

double dplm( unsigned int const l, int const m, double const x );

double dintegral_I12 ( gsl_function F, int limit, int key, double epsabs, double epsrel, gsl_integration_workspace * workspace );

double dintegral_I32 ( gsl_function F, int limit, int key, double epsabs, double epsrel, gsl_integration_workspace * workspace );
 
double dintegrand_I12 (double t, void*params);

double dintegrand_I32 (double t, void*params);

double dintegrand_Imlm32 (double t, void*params);

double dintegrand_K12 ( double t, void * params);

double dintegrand_lp32 (double t, void * params);

#endif
