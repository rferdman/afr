/* Straight summation of an array */

int ISum (int *Array, int ArrLen) {

  int i_arr;
  int ArrSum;

  ArrSum=0;
  for (i_arr=0; i_arr<ArrLen; i_arr++){
    ArrSum += Array[i_arr];
  }

  return ArrSum;
}

float FSum (float *Array, int ArrLen) {

  int i_arr;
  float ArrSum;

  ArrSum=0.;
  for (i_arr=0; i_arr<ArrLen; i_arr++){
    ArrSum += Array[i_arr];
  }

  return ArrSum;
}

double DSum (double *Array, int ArrLen) {

  int i_arr;
  double ArrSum;

  ArrSum=0.;
  for (i_arr=0; i_arr<ArrLen; i_arr++){
    ArrSum += Array[i_arr];
  }

  return ArrSum;
}
