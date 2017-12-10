#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include "init_zeta_function.h"
#include "dzeta_function_deriv.h"

/*****************************************************************************************************************
 * Zeta function
 * input:
 *   q2 - the value for q^2
 *   l, m - angular momentum quantum numbers
 *   dvec - vector d related to total linear momentum P = 2 pi/L x d
 *   gamma - gamma from Lorentz transformation, gamma = 1 / sqrt( 1 - v^2 )
 *   A - numerical factor for case of different masses
 *   epsAbs, epsRel - convergence criteria
 *****************************************************************************************************************/

int dzeta_function_deriv (double z[2], double q2, int l, int m, int*dvec, double gamma, double A, double epsAbs, double epsRel, int Lmax) {

  const double int_epsabs = 1.e-08;
  const double int_epsrel = 1.e-07;
  const size_t int_limit  = 500000;
  const size_t int_key    = 6;

  static int rotref_list[48][2][3];
  static int rotref_selection_permutation[4][4][48][3];
  static int rotref_selection_sign[4][4][48][3];
  static int rotref_number[4][4];

  static int refrot_list_is_initialized = 0;
  static int plm_norm_is_initialized = 0;

  int k1, k2, k3, status, i, irotref;
  int nvec[3];
  double qgamma    = (double)gamma;
  double qgammaInv = 1.0 / qgamma;
  double qA        = (double)A;

  gsl_integration_workspace *int_workspace = NULL;
  gsl_function F;

  int parity_factor =  1 - (l % 2);
  int lessequals_type, zeros_type;
  double int_parameters[3], int_l0m0_parameters;
  double int_integral_l0m0_value, int_integral_l0m0_error;
  double int_integral_12_value, int_integral_12_error, int_integral_32_value, int_integral_32_error;

  double qdvec[3], qnd, qddInv, qr_rrInv, qdd;
  double qterm1[2], qterm2[2], qterm1c[2], qterm2c[2], qterm3[3];
  double qdvecGammaInvMinusOne[3];
  double qint_norm_const[2], qint_norm_var[2];
  double qint_iterate_const=0.0, qint_integral_value=0.0;
  double qint_l0m0_norm_const, qint_l0m0_add;
  double qshift[3], qn[3];
  double qtmp, qtmp1, qtmp2, qw[2], qw2[2];
  double qr[3], qr_r, qr_rr, qr_costheta, qr_phi;

  if ( refrot_list_is_initialized == 0 ) {
    fprintf(stdout, "# [dzeta_function] initialize refrot_list\n");
    init_refrot_list (rotref_list, rotref_selection_permutation, rotref_selection_sign, rotref_number);
    refrot_list_is_initialized = 1;
  }
  if ( plm_norm_is_initialized == 0 ) {
    fprintf(stdout, "# [dzeta_function] initialize plm_norm\n");
    dinit_plm_norm();
    plm_norm_is_initialized = 1;
  }

  qdvec[0] = (double)dvec[0];
  qdvec[1] = (double)dvec[1];
  qdvec[2] = (double)dvec[2];
  qdd      = _DSQR(qdvec);
  qddInv   = (qdd == 0.0) ? 0.0 : 1.0 / qdd;

  qshift[0] = -0.5 * (qgammaInv * qA * qdvec[0]);
  qshift[1] = -0.5 * (qgammaInv * qA * qdvec[1]);
  qshift[2] = -0.5 * (qgammaInv * qA * qdvec[2]);

  qdvecGammaInvMinusOne[0] = qdvec[0] * ( qgammaInv - 1.0 );
  qdvecGammaInvMinusOne[1] = qdvec[1] * ( qgammaInv - 1.0 );
  qdvecGammaInvMinusOne[2] = qdvec[2] * ( qgammaInv - 1.0 );


  /* initialize first integration */
  if( (int_workspace = gsl_integration_workspace_alloc ( (size_t)(int_limit+2))) == NULL ) {
    fprintf(stderr, "# [dzeta_function] Error from gsl_integration_workspace_alloc\n");
    return(1);
  }

  qtmp = qgamma * M_PI * sqrt(M_PI);
  i = l % 4;
  if(i == 0) {
    qint_norm_const[0] = qtmp;
    qint_norm_const[1] = 0.0;
  } else if (i == 1) {
    qint_norm_const[0] = 0.0;
    qint_norm_const[1] = -qtmp;
  } else if (i == 2) {
    qint_norm_const[0] = -qtmp;
    qint_norm_const[1] = 0.0;
  } else {
    qint_norm_const[0] = 0.0;
    qint_norm_const[1] = qtmp;
  }

  /* fprintf(stdout, "# [dzeta_function] int_norm_const = %e + I %e\n", qint_norm_const[0], qint_norm_const[1]); */


  /* q2 value */
  int_parameters[0] = q2;
  
  /* initialize second integration */
  qint_l0m0_norm_const = -qgamma * M_PI * 2.0 * q2 * q2;
  qint_l0m0_add        =  qgamma * M_PI * ( 2.0 * q2 - 1.0 ) * exp(q2);
  int_l0m0_parameters  = q2;

  /* initialize real and imaginary part of term 1 and 2 */
  _DCO_EQ_ZERO( qterm1  );
  _DCO_EQ_ZERO( qterm1c );
  _DCO_EQ_ZERO( qterm2  );
  _DCO_EQ_ZERO( qterm2c );
  _DCO_EQ_ZERO( qterm3  );


  for(k3=0; k3 <= Lmax; k3++) {
    nvec[2] = k3;
  for(k2=0; k2 <= k3;    k2++) {
    nvec[1] = k2;
  for(k1=0; k1 <= k2;    k1++) {
    nvec[0] = k1;

    lessequals_type = 2*(int)(k1 < k2) + (int)(k2 < k3);
    zeros_type  = (int)(k1 ==  0) + (int)(k2 ==  0) + (int)(k3 ==  0);

    /* fprintf(stdout, "# [dzeta_function] %3d%3d%3d\t%3d%3d\t%3d\n", k1, k2, k3, lessequals_type, zeros_type, rotref_number[lessequals_type][zeros_type]); */

    for(irotref = 0; irotref < rotref_number[lessequals_type][zeros_type]; irotref++) {

      qn[0] = (double)( nvec[ rotref_selection_permutation[lessequals_type][zeros_type][irotref][0] ] * rotref_selection_sign[lessequals_type][zeros_type][irotref][0] );
      qn[1] = (double)( nvec[ rotref_selection_permutation[lessequals_type][zeros_type][irotref][1] ] * rotref_selection_sign[lessequals_type][zeros_type][irotref][1] );
      qn[2] = (double)( nvec[ rotref_selection_permutation[lessequals_type][zeros_type][irotref][2] ] * rotref_selection_sign[lessequals_type][zeros_type][irotref][2] );


      /* fprintf(stdout, "# [dzeta_function] %3d%3d%3d\t%3d\t%3.0f %3.0f %3.0f\n", k1, k2, k3, irotref, qn[0], qn[1], qn[2]); */

      qnd   = _DSCP(qn, qdvec);

      /***************
       * SECOND TERM *
       ***************/
      qtmp  = qnd * qddInv;

      qr[0] = qn[0] + qtmp * qdvecGammaInvMinusOne[0] + qshift[0];
      qr[1] = qn[1] + qtmp * qdvecGammaInvMinusOne[1] + qshift[1];
      qr[2] = qn[2] + qtmp * qdvecGammaInvMinusOne[2] + qshift[2];

      qr_rr  = _DSQR(qr);
      qr_r   = sqrt(qr_rr);

      if(qr_r < DZETA_EPSILON) {
        qr_costheta = 1.0;
        qr_phi      = 0.0;
      } else {
        qr_costheta = qr[2] / qr_r;
        qtmp = qr[0] * qr[0] + qr[1] * qr[1];

        if(qtmp < DZETA_EPSILON) {
          qr_phi = 0.0;
        } else {
          qr_phi = atan2(qr[1], qr[0]);
        }
      }

      /* fprintf(stdout, "# [dzeta_function] %e %e %e \t %e %e\t %e %e \n", qr[0], qr[1], qr[2], qr_r, qr_rr, qr_costheta, qr_phi); */

      qtmp  = (double)m * qr_phi;
      qtmp1 = pow(qr_r, (double)l ) * dplm(l, m, qr_costheta);
      qtmp2 = qr_rr - q2;
      /* fprintf(stdout, "# [dzeta_function]  %e %e  %e\n", qtmp, qtmp1, qtmp2); */

      qw[0] = cos( qtmp ) * qtmp1;
      qw[1] = sin( qtmp ) * qtmp1;

      qtmp1 = exp( -qr_rr ) / qtmp2;
      _DCO_TI_EQ_RE( qw, qtmp1 );
      /* fprintf(stdout, "# [dzeta_function] qw = %e   %e\n", qw[0], qw[1]); */

      _DKAHAN_SUM_CO_PL_CO(qterm2, qterm2c, qw);
      /* fprintf(stdout, "# [dzeta_function] term2 %3.0f %3.0f %3.0f \t %e \t %e\n", qn[0], qn[1], qn[2], qterm2[0], qterm2[1]); */



      /***************
       * FIRST TERM  *
       ***************/
      if(zeros_type < 3) {
        qtmp = qnd * qddInv * ( qgamma - 1.0 );

        qr[0] = -M_PI * ( qn[0] + qtmp * qdvec[0] ); 
        qr[1] = -M_PI * ( qn[1] + qtmp * qdvec[1] ); 
        qr[2] = -M_PI * ( qn[2] + qtmp * qdvec[2] ); 

        qr_rr    = _DSQR(qr);
        qr_rrInv = 1.0 / qr_rr;
        qr_r     = sqrt(qr_rr);

        if(qr_r < DZETA_EPSILON) {
          qr_costheta = 1.0;
          qr_phi      = 0.0;
        } else {
          qr_costheta = qr[2] / qr_r;

          qtmp = qr[0] * qr[0] + qr[1] * qr[1];
          if(qtmp < DZETA_EPSILON) {
            qr_phi = 0;
          } else {
            qr_phi = atan2( qr[1], qr[0]);
          }
        }

        /* fprintf(stdout, "# [dzeta_function] qr=(%e, %e, %e) qr_r = %e qr_costheta = %e qr_phi = %e qr_rr = %e\n", qr[0], qr[1], qr[2], qr_r, qr_costheta, qr_phi, qr_rr); */

        int_parameters[1] = (double)qr_rr;
        F.params          = (void*)int_parameters;
  
        F.function = dintegrand_12;
        status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_12_value, &int_integral_12_error);
        if(status != 0) {
          fprintf(stderr, "[dzeta_function] Error from qag, status was %d\n", status);
          return(1);
        }
        /* fprintf(stdout, "# [dzeta_function] integration 12 result   = %25.16e%25.16e\n", int_integral_12_value, int_integral_12_error); */

        F.function = dintegrand_32;
        status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_32_value, &int_integral_32_error);
        if(status != 0) {
          fprintf(stderr, "[dzeta_function] Error from qag, status was %d\n", status);
          return(1);
        }
        /* fprintf(stdout, "# [dzeta_function] integration 32 result   = %25.16e%25.16e\n", int_integral_32_value, int_integral_32_error); */

        qint_iterate_const = exp(-(qr_rr - q2));

        qtmp1 = (double)int_integral_32_value;
        qtmp2 = (double)int_integral_12_value;

        for(i=0; i<l+2; i++) {
          qtmp = ( qint_iterate_const - q2 * qtmp1 + ((double)i - 1.5)  * qtmp2 ) * qr_rrInv;
          qtmp1 = qtmp2;
          qtmp2 = qtmp;
        }
        qint_integral_value = qtmp2;
        /* fprintf(stdout, "# [dzeta_function] integral value          = %25.16e\n", qint_integral_value); */

        /* TEST */
/*
        F.function = dintegrand_lp32;
        int_parameters[2] = (double)l;
        status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_32_value, &int_integral_32_error);
        if(status != 0) {
          fprintf(stderr, "[dzeta_function] Error from qag, status was %d\n", status);
          return(1);
        }
        fprintf(stdout, "# [dzeta_function] integration lp32 result = %25.16e%25.16e\n", int_integral_32_value, int_integral_32_error);
*/

        qtmp  = (double)m * qr_phi;
        qtmp1 = pow(qr_r, l) * dplm (l, m, qr_costheta);

        qw[0] = qtmp1 * qint_integral_value;
        qw[1] = sin( qtmp ) * qw[0];
        qw[0] *= cos( qtmp );
        /* fprintf(stdout, "# [dzeta_function] qw = %25.16e + I %25.16e\n", qw[0], qw[1]); */

        qtmp = M_PI * qA * qnd;
        qint_norm_var[0] =  cos(qtmp) * (double)      parity_factor;
        qint_norm_var[1] =  sin(qtmp) * (double)( 1 - parity_factor );
        /* fprintf(stdout, "# [dzeta_function] qint_norm_var = %25.16e + I %25.16e\n", qint_norm_var[0], qint_norm_var[1]); */

            
        _DCO_EQ_CO_TI_CO(qw2, qw, qint_norm_var);
        /* fprintf(stdout, "# [dzeta_function] qw2 = %25.16e + I %25.16e\n", qw2[0], qw2[1]); */

        _DKAHAN_SUM_CO_PL_CO(qterm1, qterm1c, qw2);


        /* fprintf(stdout, "# [dzeta_function] term1 %3.0f %3.0f %3.0f \t %25.16e \t %25.16e\n", qr[0], qr[1], qr[2], qterm1[0], qterm1[1]); */
      }  /* end of if zeros_type < 3 */

    }  /* end of loop on rotations and reflections */

  }  /* end of loop on k1 */
  }  /* end of loop on k2 */

    /* fprintf(stdout, "# [dzeta_function] qconvergence %3d \t %25.16e \t %25.16e\n", k3, qterm2[0], qterm2[1]); */

  }  /* end of loop on k3 */

  /* subtraction for l = 0 and m = 0 */
  if(l == 0 && m == 0) {

    F.function = dint_l0m0_kernelFunction;
    F.params   = (void*)(&int_l0m0_parameters);
    status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_l0m0_value, &int_integral_l0m0_error);
    if(status != 0) {
      fprintf(stderr, "[dzeta_function] Error from qag, status was %d\n", status);
      return(1);
    }
    qtmp = qint_l0m0_norm_const * (double)int_integral_l0m0_value + qint_l0m0_add;
    /* TEST */
    /* fprintf(stdout, "# [dzeta_function] (l=0 m=0) modified integral value %16.7e  %25.16e\n", q2, qtmp); */
    qterm3[0] = qtmp;
    qterm3[1] = 0.;
  }

  if(int_workspace != NULL) {
    gsl_integration_workspace_free(int_workspace);
  }

  /* multiply qterm2 with exp( q2 ) */
  qtmp = exp(  q2 );
  _DCO_TI_EQ_RE(qterm2, qtmp);

  /* multiply qterm1 with constant integral normalization */
  _DCO_EQ_CO(qw, qterm1);
  _DCO_EQ_CO_TI_CO(qterm1, qw, qint_norm_const);

  z[0] = qterm1[0] + qterm2[0] + qterm3[0];  
  z[1] = qterm1[1] + qterm2[1] + qterm3[1];
  
  /* TEST */
/*
  fprintf(stdout, "# [zeta_function] term1 %e + I %e\n", qterm1[0], qterm1[1]);
  fprintf(stdout, "# [zeta_function] term2 %e + I %e\n", qterm2[0], qterm2[1]);
  fprintf(stdout, "# [zeta_function] term3 %e + I %e\n", qterm3[0], qterm3[1]);
*/

#if 0
#endif  /* of if 0 */
  return(0);
}  /* end of dzeta_function_deriv */
