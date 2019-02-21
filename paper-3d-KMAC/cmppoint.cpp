
#include "cmppoint.h"
int cmppoint(void* a,void* b){
    if (((point_avl*)a)->x > ((point_avl*)b)->x)
     {return 1;}
     else
     if (((point_avl*)a)->x < ((point_avl*)b)->x)
     {return -1;}
     else 
     return 0;
}
