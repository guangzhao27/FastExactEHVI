#include "avl.h"

#ifndef _HYCON3D_H
#define _HYCON3D_H 1
#define MAXDOUBLE 100000000000.0
#define MAXPOINTS 20000
struct point_avl:public avl{
        double x;
	double y;
	double z;
	int myindex;
	int originalindex;
};



struct box {
       double lx;
       double ly;
       double lz;
       double ux;
       double uy;
       double uz;
};

int contributions(double* results, point_avl* mypoint, int k, double* ref);
// Standard includes
#include <cstdlib>
#include <iostream>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include "hycon3d.h"
#include "boxlist.h"
#include "getLRlist.h"
#include "printlist.h"
#include "listPointTree.h"
#include "quicksortPZ.h"
#include "volBox.h"
#include "cmppoint.h"

#define DEBUGMODUS 1
#define PRINTCHECK(A) printf("\n__________ CHECKPOINT __________\n"); printf(A); printf("\n--------------------------------\n")
#include "boxlist.h"

#endif
