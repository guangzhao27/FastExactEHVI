#include "volBox.h"

double volBox(box b)
{ 
      if ((b.ux - b.lx)*(b.uy-b.ly)*(b.uz-b.lz) < 0)
      {
        printf("Achtung! Achtung! Box Volume is negative!\n");
        printf("Faulty coordinates:");
        printf("%f, %f, %f, %f, %f, %f", 
                      b.lx, b.ly, b.lz, 
                      b.ux, b.uy, b.uz);
        getchar();
      }; 
      return (b.ux - b.lx)*(b.uy-b.ly)*(b.uz-b.lz);  
           
}
