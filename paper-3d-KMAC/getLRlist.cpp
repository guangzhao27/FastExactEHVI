#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include"cmppoint.h"
#include"listPointTree.h"
#include "avl.h"
#include "boxlist.h"

// Computes a list of point indexes with c elements, starting from the point 
// that are just bigger in the x-coordinate than x_up ending in the last point
// in the y-coordinate that is smaller than y_up. 
// It uses the avl routines. The initial node 'a' must be the root when being
// called. The index must be c. Make sure the array inner is big enough to 
// comprise the number of nodes in the tree.

// int* inner (array): Record of all points inside the range. 
// int* right: The element that is just bigger in the x coordinate
// int* left: The element that is just bigger in the y coordinate

int getLRlist(int *inner, int* left, int* right, 
              point_avl* a, 
              int c, double y_up, double x_up, double* left_best, double* right_best)
{         
  if (a==0) return c;

  // conditional update of the rightmost point
  if (a->x > x_up && a->x <= *right_best)  
  {
           *right_best = a->x;
           *right = a->myindex;
           }
           
  if (a->right)
  {
    point_avl* aright = (point_avl*)(a->right);
                
    // go to the right subtree rooted in a->right
    if (a->x <= x_up)
      {c=getLRlist(inner, left, right, aright, c, y_up, x_up, left_best, right_best);}
  }
  // record a if it is in the feasible range
  if ((a->y <= y_up) && (a->x <= x_up)) 
  {   
      inner[c]=a->myindex; c++;
  } 
      
  // conditional update of the leftmost point
  if (a->y > y_up && a->y <= *left_best)
  {
           *left_best = a->y;
           *left= a->myindex;
           }
  
  //printf("Visiting: %f %f index=%i ...\n",a->x,a->y, a->myindex);
  if (a->left)
  {
    point_avl* aleft = (point_avl*)(a->left);      
              
    // go to the left subtree rooted in a->left
    if (a->y <= y_up)
      {c=getLRlist(inner, left, right, aleft, c, y_up, x_up, left_best, right_best);}
  } 
  return c;
}



/*
main()
{
    int i;
    avl_tree points;
    point_avl mypoint[MAXPOINTS];
    points.compar=cmppoint;
	points.root=0;
	
	//mypoint[3].x=1;
	//mypoint[3].y=3;
	mypoint[1].x=1;
	mypoint[1].y=3;
	mypoint[2].x=1.25;
	mypoint[2].y=2.75;
	mypoint[0].x=4;
	mypoint[0].y=2;
	//mypoint[4].x=0;
	//mypoint[4].y=10;

	for(i=0;i<3;i++) mypoint[i].myindex=i;
	
	// Build a balanced search tree
	for(i=0;i<3;i++)
    {
		printf("-------------\n");
		avl_insert(&points,(avl*)&mypoint[i]);
		listPointTree((point_avl*)points.root,0);
	}
	
	int inner[5];
	int left=-1;
	int right=-1;
	int counter=0;
	double x_up=1.5;
	double y_up=2.85;
    printf("TESTING getLRlist: ...\n");

    printf("upper y=%f bounds left point.\n",y_up);
    printf("upper x=%f bounds right point.\n",x_up);  
    double right_best=MAXDOUBLE; double left_best=MAXDOUBLE;
    counter = getLRlist(inner, &left, &right, (point_avl*)points.root, counter, y_up, x_up, &left_best, &right_best);
    
    printf("Right beyond boundary point: %i with coordinates (%f, %f)\n", 
      right, mypoint[right].x, mypoint[right].y);

    printf("Inner points \n");
    for (i=0; i < counter; i++)
        printf("Inner point: %i with coordinates (%f, %f)\n", 
          inner[i], mypoint[inner[i]].x, mypoint[inner[i]].y);
    printf("\n");
    printf("Left beyond boundary point: %i with coordinates (%f, %f)\n", 
      left, mypoint[left].x, mypoint[left].y);
    
    getchar();
    exit(0);
    
}
*/



