#include <stdio.h>
#include <stdlib.h>

typedef struct ThisStruct{
  int nx;
  int ny;
  double *x;
  double *y;
} ThisStruct;

int thisfunction(int nstruct, ThisStruct *thesestructs)
{

  for(int i=0; i<nstruct; i++){

    ThisStruct *thisstruct = thesestructs+i;
    
    printf("nx %d\n", thisstruct->nx);
    printf("ny %d\n", thisstruct->ny);
    
    for(int i=0; i<10; i++) printf("x %lf\n", thisstruct->x[i]);
    for(int i=0; i<10; i++) printf("y %lf\n", thisstruct->y[i]);
  }
  
  return EXIT_SUCCESS;
}
