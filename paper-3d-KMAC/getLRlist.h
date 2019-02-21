#include "hycon3d.h"

#ifndef _GETLRLIST_H
#define _GETLRLIST_H 1

int getLRlist(int *inner, int* left, int* right, 
              point_avl* a, 
              int c, double y_up, double x_up, double *best_right, double *best_left);

void test_getLRList();
#endif
