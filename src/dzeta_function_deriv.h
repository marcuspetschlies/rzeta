#ifndef _DZETA_FUNCTION_DERIV_H
#define _DZETA_FUNCTION_DERIV_H


#ifdef DZETA_USE_KAHAN

#define _DKAHAN_SUM_CO_PL_CO(_a, _b, _c) { \
  double _ks_y, _ks_t;                 \
  _ks_y  = (_c)[0] - (_b)[0];              \
  _ks_t  = (_a)[0] + _ks_y;                \
  (_b)[0] = (_ks_t - (_a)[0] ) - _ks_y;    \
  (_a)[0] = _ks_t;                         \
  _ks_y  = (_c)[1] - (_b)[1];              \
  _ks_t  = (_a)[1] + _ks_y;                \
  (_b)[1] = (_ks_t - (_a)[1] ) - _ks_y;    \
  (_a)[1] = _ks_t;                         \
}
    
#define _DKAHAN_SUM_CO_PL_RE(_a, _b, _r) { \
  double _ks_y, _ks_t;                 \
  _ks_y  = (_r) - (_b)[0];                 \
  _ks_t  = (_a)[0]  + _ks_y;               \
  (_b)[0] = (_ks_t - (_a)[0] ) - _ks_y;    \
  (_a)[0] = _ks_t;                         \
}

#else

#define _DKAHAN_SUM_CO_PL_CO(_a, _b, _c) { \
  (_a)[0] += (_c)[0];\
  (_a)[1] += (_c)[1];\
}
    
#define _DKAHAN_SUM_CO_PL_RE(_a, _b, _r) { \
  (_a)[0] += (_r);\
}

#endif  /* of if def DZETA_USE_KAHAN */


#define DZETA_EPSILON 1.e-14

#define _DSQR(_x) ((_x)[0]*(_x)[0]+ (_x)[1]*(_x)[1]+ (_x)[2]*(_x)[2])

#define _DSCP(_x, _y) ((_x)[0]*(_y)[0] + (_x)[1]*(_y)[1] + (_x)[2]*(_y)[2])

#define _DCO_EQ_CO_TI_CO(_c, _a, _b) { \
  (_c)[0] = (_a)[0] * (_b)[0] - (_a)[1] * (_b)[1];\
  (_c)[1] = (_a)[0] * (_b)[1] + (_a)[1] * (_b)[0];\
}

#define _DCO_EQ_CO_TI_RE(_c, _a, _r) {\
  (_c)[0] = (_a)[0] * (_r);           \
  (_c)[1] = (_a)[1] * (_r);           \
}

#define _DCO_PL_EQ_CO_TI_RE(_c, _a, _r) {\
  (_c)[0] += (_a)[0] * (_r);             \
  (_c)[1] += (_a)[1] * (_r);             \
}

#define _DCO_TI_EQ_RE(_c, _r) { \
  (_c)[0] *= (_r);              \
  (_c)[1] *= (_r);              \
}

#define _DCO_EQ_ZERO(_c)  {\
  (_c)[0] = 0.0;          \
  (_c)[1] = 0.0;          \
}

#define _DCO_EQ_CO(_a, _b) {\
  (_a)[0] = (_b)[0];        \
  (_a)[1] = (_b)[1];        \
}

#define _PLM_L_MAX 4

extern double plm_norm[(_PLM_L_MAX+1)*(_PLM_L_MAX+1)];

double dplm( unsigned int l, int m, double x);

#endif
