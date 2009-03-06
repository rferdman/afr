/* Straight summation of an array */

int ISum (int *Array, int ArrLen, int *ArrSum) {

  int i_arr;

  ArrSum=0;
  for (i_arr=0; i_arr<ArrLen; i_arr++){
    *ArrSum += Array[i_arr];
  }

  return 0;
}

int FSum (float *Array, int ArrLen, float *ArrSum) {

  int i_arr;

  ArrSum=0.;
  for (i_arr=0; i_arr<ArrLen; i_arr++){
    *ArrSum += Array[i_arr];
  }

  return 0;
}

int DSum (double *Array, int ArrLen, double *ArrSum) {

  int i_arr;

  ArrSum=0.;
  for (i_arr=0; i_arr<ArrLen; i_arr++){
    *ArrSum += Array[i_arr];
  }

  return 0;
}
