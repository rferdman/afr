/* Basic shell sort; obtained from 
   http://linux.wku.edu/~lamonml/algor/sort/shell.html */

/* Performs in-place sorting of input array */

void ISort(int n_array, int array[])
{
  int   i, j, increment;
  int   temp;

  increment = 3;
  while (increment > 0)
  {
    for (i=0; i < n_array; i++)
    {
      j = i;
      temp = array[i];
      while ((j >= increment) && (array[j-increment] > temp))
      {
        array[j] = array[j - increment];
        j = j - increment;
      }
      array[j] = temp;
    }
    if (increment/2 != 0)
      increment = increment/2;
    else if (increment == 1)
      increment = 0;
    else
      increment = 1;
  }
}


void FSort(int n_array, float array[])
{
  int   i, j, increment;
  float temp;

  increment = 3;
  while (increment > 0)
  {
    for (i=0; i < n_array; i++)
    {
      j = i;
      temp = array[i];
      while ((j >= increment) && (array[j-increment] > temp))
      {
        array[j] = array[j - increment];
        j = j - increment;
      }
      array[j] = temp;
    }
    if (increment/2 != 0)
      increment = increment/2;
    else if (increment == 1)
      increment = 0;
    else
      increment = 1;
  }
}


void DSort(int n_array, double array[])
{
  int    i, j, increment;
  double temp;

  increment = 3;
  while (increment > 0)
  {
    for (i=0; i < n_array; i++)
    {
      j = i;
      temp = array[i];
      while ((j >= increment) && (array[j-increment] > temp))
      {
        array[j] = array[j - increment];
        j = j - increment;
      }
      array[j] = temp;
    }
    if (increment/2 != 0)
      increment = increment/2;
    else if (increment == 1)
      increment = 0;
    else
      increment = 1;
  }
}


