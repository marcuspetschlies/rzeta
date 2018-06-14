#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include "init_zeta_function.h"
#include "dzeta_function_utils.h"
/*******************************************************************************
 * dvec = integer vector, P = 2 pi d / L
 *
 * P^m_l(cos(theta)) 
 *  double gsl_sf_legendre_sphPlm (int l, int m, double x) 
 *
 * check convergence with epsAbs and epsRel
 *
 * typedef for integration function double (*intfunction)    (double x, void*params)
 *
 * We characterize the n-vectors into 8 clases; we use 0 < h < k < l
 *     repres.  0-type <=-type
 * (0) (0,0,0)  3      0
 * ======================
 * (1) (0,0,h)  2      1
 * ======================
 * (2) (0,h,h)  1      2
 * ----------------------
 * (3) (0,h,k)  1      3
 * ======================
 * (4) (h,h,h)  0      0
 * ----------------------
 * (5) (h,h,k)  0      1
 * ----------------------
 * (6) (h,k,k)  0      2
 * ----------------------
 * (7) (h,k,l)  0      3
 *
 *  normalized associated Legendre polynomicals with negative m
 *
 *  sqrt{(2l+1)/(4 pi)}  x sqrt{(l-m)/(l+m) } P_l^m =
 *    sqrt{(2l+1)/(4 pi)}  x sqrt{(l+|m|)/(l-|m|) } P_l^{-|m|} =
 *    sqrt{(2l+1)/(4 pi)}  x sqrt{(l+|m|)/(l-|m|) } x (-1)^m  (l-|m|)/(l+|m|)  P_l^{|m|}
 *    sqrt{(2l+1)/(4 pi)}  x sqrt{(l-|m|)/(l+|m|) } x (-1)^m  P_l^{|m|}
 *    (-1)^m x associated Legendre polynomial with (l, |m|)
 *******************************************************************************/

/***********************************************************
 * angular momentum functions
 ***********************************************************/

double plm_norm[(_PLM_L_MAX+1)*(_PLM_L_MAX+1)];

int dinit_plm_norm (void) {

  int k, l, m;
  int idxp, idxm, idx0;
  double const qone_over_sqrt4pi = 0.25 * M_2_SQRTPI;

  for( int l=0; l <= _PLM_L_MAX; l++ ) {

    int const idx0 = l * (l + 1);
    plm_norm[idx0] = sqrt(2. * l + 1.) * qone_over_sqrt4pi;

    for ( int m = 1; m<=l; m++) {

      int const idxm = idx0 - m;
      int const idxp = idx0 + m;

      plm_norm[idxp] = plm_norm[idx0];

      double z = (double)(l+m);
      for ( int k=l+m-1; k>l-m; k--) {
        z *= (double)k;
      }
      plm_norm[idxp] = plm_norm[idx0] / sqrt( z );

      plm_norm[idxm] =  (double)(1 - 2*(m%2)) * plm_norm[idxp];
    }
  }

  /* TEST */
/*
  for ( int l=0; l <= _PLM_L_MAX; l++) {
    for( int m = -l; m <= l; m++) {
      int const idx0 = l*l + (l+m);
      fprintf(stdout, "# [dinit_plm_norm] l=%3d m=%3d idx=%3d norm=%25.16e\n", l, m, idx0, plm_norm[idx0]);
    }
    fprintf(stdout, "# [dinit_plm_norm] ----------------------------------------------------------------------------\n");
  }
*/
  /* END OF TEST */

  return(0);
}

/***************************************************************************/
/***************************************************************************/

/***************************************************************************
 * normalized associated Legendre polynomials
 * 
 * - input: l and m; l >= 0 and -l <= m <= l
 *          x with  -1 <= x <= 1; x = cos(theta)
 * - output: sqrt{(2l+1)/(4 pi)} x sqrt{(l-m)/(l+m) } P_l^m (x)
 * - explicit form of P_l^m from
 *     https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 ***************************************************************************/
double dplm( unsigned int const l, int const m, double const x ) {

  int const idx = l * ( l + 1 ) + m;
  int const mabs = m<0 ? -m : m;

  double const norm = plm_norm[idx];
  double const xx = x * x;


  if(l == 0) {
    return( norm );
  } else if(l == 1 ) {
    if(mabs == 0) {
      return( x * norm ) ;
    } else if (mabs == 1) {
      return( -sqrt( 1. - xx ) * norm ) ;
    }
  } else if( l == 2 ) {
    if( mabs == 0 ) {
      return(  0.5 * (3.*xx - 1. ) * norm );
    } else if (mabs == 1) {
      return( -3. * x * sqrt( 1. - xx ) * norm );
    } else if (mabs == 2) {
      return(  3. * ( 1. - xx ) * norm );
    }
  } else if( l == 3 ) {
    double const s = sqrt( 1. - xx );
    if( mabs == 0 ) {
      return ( 0.5 * x * ( 5.*xx - 3. ) * norm );
    } else if(mabs == 1) {
      return( -1.5 * (5.*xx - 1.) *  s * norm );
    } else if(mabs == 2) {
      return( -15. * x * (-1. + xx ) * norm ) ;
    } else if(mabs == 3) {
      return( -15. * s * s * s * norm );
    }
  } else if(l == 4) {
    double s = sqrt(1. - xx);
    if(mabs == 0) {
      return ( 0.125 * ( ( 35.*xx -30. ) * xx +3. ) * norm );
    } else if (mabs == 1) {
      return ( -2.5 * x * (7. * xx - 3.) * s * norm );
    } else if (mabs == 2) {
      return ( -7.5 * (7. * xx - 1.) * (-1. + xx) * norm );
    } else if (mabs == 3) {
      return(  -105. * x * (1. - xx) * s * norm );
    } else if (mabs == 4) {
      return ( 105. * (-1. + xx) * (-1. + xx) * norm );
    }
  } else if ( l == 5 ) {
    double s = sqrt ( 1. - xx );
    if(mabs == 0) {
      return ( 0.125 * x * ( ( 63*xx - 70 ) * xx +15 ) * norm );
    } if(mabs == 1) {
      return ( -1.875 * s * ( ( 21*xx - 14 ) * xx + 1 ) * norm );
    } if(mabs == 2) {
      return ( -52.5 * x * ( -1 + xx ) * ( -1 + 3*xx ) * norm );
    } if(mabs == 3) {
      return ( 52.5 * s * ( -1 + xx ) * ( -1 + 9*xx ) * norm );
    } if(mabs == 4) {
      return ( 945 * x * (1 - xx) * (1 - xx) * norm );
    } if(mabs == 5) {
      return ( -945 * ( 1 - xx )*(1 - xx) * s * norm );
    }
  } else if ( l == 6 ) {
    double s = sqrt ( 1. - xx );
    if(mabs == 0) {
      return ( 0.0625 * ( -5 + xx * ( 105 + xx * (-315 + 231 * xx )  )) * norm );
    } if(mabs == 1 ) {
      return ( -2.625 * s * x * ( 5 + xx * ( -30. + 33. * xx  )) * norm );
    } if(mabs == 2 ) {
      return ( -13.125 * ( - 1 + xx ) * ( 1 + xx * ( -18. + 33. * xx  )) * norm );
    } if(mabs == 3 ) {
      return ( 157.5 * s * x * ( - 1 + xx ) * ( -3. + 11. * xx ) * norm );
    } if(mabs == 4 ) {
      return ( 472.5 * ( - 1 + xx ) * ( -1 + xx ) * ( -1. + 11. * xx ) * norm );
    } if(mabs == 5 ) {
      return ( -10395 * ( - 1 + xx ) * ( -1 + xx ) * s * norm );
    } if(mabs == 6 ) {
      return ( -10395 * ( - 1 + xx ) * ( -1 + xx ) * ( -1 + xx ) * norm );
    }
  } else {
    return( sqrt(-1.) );
  }
}

/**************************************************************************
 * integrands
 **************************************************************************/

double dintegrand_I12 (double t, void*params) {
  double q2   = (double)((double*)params)[0];
  double a    = (double)((double*)params)[1];
  return( exp( q2 * t - a / t) *sqrt(t) ); 
}  /* end of dintegrand_I12 */

double dintegrand_I32 (double t, void*params) {
  double q2   = (double)((double*)params)[0];
  double a    = (double)((double*)params)[1];
  return( exp( q2 * t - a / t) *sqrt(t) * t ); 
}  /* end of dintegrand_I32 */

double dintegrand_Imlm32 (double t, void*params) {
  double q2   = (double)((double*)params)[0];
  double a    = (double)((double*)params)[1];
  double l    = (double)((double*)params)[2];
  return( exp( q2 * t - a / t)  / ( sqrt(t) * pow(t,l+1) ) ); 
}  /* end of dintegrand_mlm32 */


double dintegrand_K12 (double t, void*params) {
  double q2   = (double)((double*)params)[0];
  return( exp(t*q2) * sqrt( t ) );
}  /* end of dintegrand_K12 */


double dintegrand_Jm32 (double t, void*params) {
  double q2  = (double)((double*)params)[0];
  return(  ( exp(t*q2) - 1. ) / ( t * sqrt(t) ) );
}  /* end of dintegrand_Jm32 */


/**************************************************************************/
/**************************************************************************/


/**************************************************************************
 * I12 integral values
 **************************************************************************/
double dintegral_I12 ( gsl_function F, int limit, int key, double epsabs, double epsrel, gsl_integration_workspace * workspace ) {

  const gsl_function F;
  F.function = dintegrand_I12;
  F.params   = (void*)parameters;

  int exitstatus;
  double value, error;

  exitstatus = gsl_integration_qag (&F, 0., 1., epsabs, epsrel, limit, key, workspace, &value, &error);
  if( exitstatus != 0 ) {
    fprintf(stderr, "[dintegral_I12] Error from qag, status was %d\n", status);
    return(sqrt(-1.));
  }
  /* fprintf(stdout, "# [dintegral_I12] integration 12 result   = %25.16e%25.16e\n", value, error); */

  return(value);

}  // end of function dintegral_I12

/**************************************************************************/
/**************************************************************************/

/**************************************************************************
 * I32 integral values
 **************************************************************************/
double dintegral_I32 ( gsl_function F, int limit, int key, double epsabs, double epsrel, gsl_integration_workspace * workspace ) {

  int exitstatus;
  double value, error;

  exitstatus = gsl_integration_qag ( &F, 0., 1., epsabs, epsrel, limit, key, workspace, &value, &error);
  if( exitstatus != 0 ) {
    fprintf(stderr, "[dintegral_I32] Error from qag, status was %d\n", status);
    return(sqrt(-1.));
  }
  /* fprintf(stdout, "# [dintegral_I32] integration 12 result   = %25.16e%25.16e\n", value, error); */

  return(value);

}  // end of function dintegral_I32


/**************************************************************************/
/**************************************************************************/

/**************************************************************************
 * Imlm32 integral values
 **************************************************************************/
double dintegral_Imlm32 ( gsl_function F, int limit, int key, double epsabs, double epsrel, gsl_integration_workspace * workspace ) {

  double const q2    = ((double*)F.params)[0];
  double const a2    = ((double*)F.params)[1];
  double const a2inv = 1. / a2;
  double const dtmp3 = exp( -( a2 - q2 ) );
  double dtmp1 = dintegral_I32 ( F, limit, key, epsabs, epsrel, workspace );
  double dtmp2 = dintegral_I12 ( F, limit, key, epsabs, epsrel, workspace );
    

  for( int i = 0; i < l+2; i++) {
    double dtmp = ( dtmp3  - q2 * dtmp1 + ((double)i - 1.5)  * dtmp2 ) * a2inv;
    dtmp1 = dtmp2;
    dtmp2 = dtmp;
  }
  /* fprintf(stdout, "# [dintegral_Imlm32] integral value = %25.16e\n", dtmp2 ); */

  return( dtmp2 );
}  // end of function dintegral_Imlm32
