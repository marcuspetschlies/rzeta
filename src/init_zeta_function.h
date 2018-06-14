#ifndef _INIT_ZETA_FUNCTION_H
#define _INIT_ZETA_FUNCTION_H

inline int _V3_EQ_V3( int const a[3], int const b[3] ) {
  return ( ( a[0]==b[0] ) && ( a[1]==b[1] ) && ( a[2]==b[2] ) );
}


int init_refrot_list (int rotref_list[48][2][3], int rotref_selection_permutation[4][4][48][3], int rotref_selection_sign[4][4][48][3], int rotref_number[4][4]);
#endif
