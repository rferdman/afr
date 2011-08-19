/* Little routine to test for file existence */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int FileExists(char *testfile)
{
  
  FILE *Ftest;
  
  if((Ftest=fopen(testfile, "r"))==NULL){
    /* File cannot be opened */
    return 0;
  }
  else{
    /* File exists! */
    fclose(Ftest);
    return 1;
  }
  
}

