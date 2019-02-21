
// Printlist for MATLAB
#include <stdio.h>
#include "hycon3d.h"
#include "printlist.h"

void printlist(point_avl list[],int n) 
{ 
  int i;    
  printf("L = [...\n");
  for(i=0;i<n;i++)       
  {
    printf("   %f, %f, %f",
       list[i].x, list[i].y, list[i].z); 
    if (i < n-1) printf(";...\n");
  }
  printf("]; drawZeroBoxes(L);");
} 
