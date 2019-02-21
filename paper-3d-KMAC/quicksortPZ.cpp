// Quicksort
#include "quicksortPZ.h"
#include "hycon3d.h"


void swap(point_avl *x,point_avl *y) 
{ 
  point_avl temp; 
  temp = *x; 
  *x = *y;    
  *y = temp; 
}  
 
int choose_pivot(int i,int j ) 
{    
   return((i+j) /2); 
} 

 
void quicksortPZ(point_avl list[],int m,int n) 
{ 
   point_avl key;
   int i,j,k;    
   if( m < n)    
   { 
      k = choose_pivot(m,n);     
      swap(&list[m],&list[k]);    
      key = list[m];    
      i = m+1;    
      j = n;    
      while(i <= j)    
      {       
           while((i <= n) && (list[i].z >= key.z))         
           i++;  
           while((j >= m) && (list[j].z < key.z))         
           j--;  
      if( i < j) swap(&list[i],&list[j]);
     }     
     // swap two elements     
     swap(&list[m],&list[j]); 

     // recursively sort the lesser list 
     quicksortPZ(list,m,j-1);     
     quicksortPZ(list,j+1,n); 
   } 
} 

