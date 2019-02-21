// The main file for Calculation of EHVI for 3-D case, using nlogn method
// C++ version
// (C) by Kaifeng Yang, Michael T. M. Emmerich, 
// k.yang@liacs.leidenuniv.nl
// m.t.m.emmerich@liacs.leidenuniv.nl
// Date: Nov. 21, 2015
// For details, please refer to README.txt in the current folder
#include <vector>
#include <deque>
#include <algorithm>
#include <math.h>
#include "ehvi_hvol.h"
#include "ehvi_multi.h"
#include "avl.h"
#include <iostream>
#include "helper.h"
#include "ehvi_consts.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmppoint.h"
#include "listPointTree.h"
#include "avl.h"
#include "boxlist.h"
//#include <omp.h>


// The parameter stores the boxes' information
static box layer_box[1001][1001];
static int test;

// Get the lower and upper bound information for each boxes
struct box box_info(point_avl* crt_p_l, point_avl* crt_p_r, point_avl* crt_p, double r[], int iter){
  box temp_box;
  if (iter==0){
    temp_box.lx = r[0];
    temp_box.ly = r[1];
    temp_box.lz = crt_p->z;
    temp_box.ux = crt_p->x;
    temp_box.uy = crt_p->y;
    temp_box.uz = INFINITY;
  }
  else{
    // The new point is added to the right 
    if (crt_p_r->x==0 && crt_p_l->x!=0){
      temp_box.lx = crt_p_l->x;
      temp_box.ly = r[1];
      temp_box.lz = crt_p->z;
      temp_box.ux = crt_p->x;
      temp_box.uy = crt_p->y;
      temp_box.uz = INFINITY;
    }
    // The new point is added to the left
    else if (crt_p_l->x==0 && crt_p_r->x!=0){
      temp_box.lx = r[0];
      temp_box.ly = crt_p_r->y;
      temp_box.lz = crt_p->z;
      temp_box.ux = crt_p->x;
      temp_box.uy = crt_p->y;
      temp_box.uz = INFINITY;
    }
    // The new point is added between 2 points
    else {
      temp_box.lx = crt_p_l->x;
      temp_box.ly = crt_p_r->y;
      temp_box.lz = crt_p->z;
      temp_box.ux = crt_p->x;
      temp_box.uy = crt_p->y;
      temp_box.uz = INFINITY;
    }
  }
  return temp_box;
}

// Calculation for EHVI for each hyper boxes
vector<double> calculation_box(box a, vector<mus *> & pdf, double r[], int j, int n_p){
  double eq1,eq2,eq3,eq4;
  double ep_u1_u1,ep_u1_l1,ep_u2_u2,ep_u2_l2,ep1,ep2,ep3;
  int n=pdf.size();
  int i;
  
  //double asw;
  vector<double> asw;
  while (asw.size() < n)
        asw.push_back(0);
   
  for(i=0;i<n;i++){
    ep_u1_u1 = exipsi(a.lx,a.ux,pdf[i]->mu[0],pdf[i]->s[0]);
    ep_u1_l1 = exipsi(a.lx,a.lx,pdf[i]->mu[0],pdf[i]->s[0]);
    ep_u2_u2 = exipsi(a.ly,a.uy,pdf[i]->mu[1],pdf[i]->s[1]);
    ep_u2_l2 = exipsi(a.ly,a.ly,pdf[i]->mu[1],pdf[i]->s[1]);
    
    ep3 = exipsi(a.lz,a.lz,pdf[i]->mu[2],pdf[i]->s[2]) - exipsi(a.lz,a.uz,pdf[i]->mu[2],pdf[i]->s[2]);
    eq1 = (ep_u1_u1-ep_u1_l1)*(ep_u2_u2-ep_u2_l2)*ep3;
    eq2 = -(ep_u1_u1-ep_u1_l1)*(a.uy-a.ly)*(1-gausscdf(a.uy, pdf[i]->mu[1], pdf[i]->s[1]))*ep3;
    eq3 = -(a.ux-a.lx)*(1-gausscdf(a.ux, pdf[i]->mu[0], pdf[i]->s[0]))*(ep_u2_u2-ep_u2_l2)*ep3;
    eq4 = (a.ux-a.lx)*(1-gausscdf(a.ux, pdf[i]->mu[0], pdf[i]->s[0]))*(a.uy-a.ly)*(1-gausscdf(a.uy, pdf[i]->mu[1], pdf[i]->s[1]))*ep3;
    if (j<n_p)
	    asw[i] = eq1+eq2+eq3+eq4+asw[i];
    else{ 
      if(a.ux==INFINITY)
        ep1 = exipsi(a.lx,a.lx,pdf[i]->mu[0],pdf[i]->s[0]) - exipsi(a.lx,a.ux,pdf[i]->mu[0],pdf[i]->s[0]);
      else
        ep1 = exipsi(a.lx,a.lx,pdf[i]->mu[0],pdf[i]->s[0]) - exipsi(a.ux,a.ux,pdf[i]->mu[0],pdf[i]->s[0]);
      if(a.uy==INFINITY)
        ep2 = exipsi(a.ly,a.ly,pdf[i]->mu[1],pdf[i]->s[1]) - exipsi(a.ly,a.uy,pdf[i]->mu[1],pdf[i]->s[1]);
      else
        ep2 = exipsi(a.ly,a.ly,pdf[i]->mu[1],pdf[i]->s[1]) - exipsi(a.uy,a.uy,pdf[i]->mu[1],pdf[i]->s[1]);

	    asw[i] = asw[i] + ep1*ep2*ep3;
    }
  }
 return asw;
}

/* This function is nlogn method for 3-d case
 * The input parameters are:
 *                P:      The current Pareto front
 *                r:      reference point
 *                pdf:    The evaluated points
*/
vector<double> ehvi3d_nlogn(deque<individual*> P, double r[], vector<mus *> & pdf){
  vector<double> answer,temp_answer;                //The answer
  int n_eval=pdf.size();                //The number of evaluated points
  int n_p=P.size();			// The number of current pf
  avl_tree points;
  point_avl mypoint[n_p];
  point_avl temp_point;
  points.compar=cmppoint;
  points.root=0;
  avl_tree find_node;
  int i,j,k;
  int inner[n_p];
  int left=-1;
  int right=-1;
  int counter=0;
  
  double right_best=MAXDOUBLE; double left_best=MAXDOUBLE;
  double box_test[n_p+1]; 
  double box_index[n_p+1];      // the value is the number of the box
  
  for(int i=0;i<n_eval;i++)   {answer.push_back(0); temp_answer.push_back(0);}
  // Sort by z value
  sort(P.begin(), P.end(), zcomparator);
  // Initilize the nodes
  for (i=0;i<n_p;i++){
    mypoint[i].x = P[i]->f[0];
    mypoint[i].y = P[i]->f[1];
    mypoint[i].z = P[i]->f[2];
    mypoint[i].myindex = i;
  }

  // Get lower bound and upper bound 
  // Three cases:
  // 1. a point with the largest z value
  // 2. a point with the z value of 0
  // 3. points with the z value between largest and 0
  counter = 0;
  avl_insert(&points,(avl*)&mypoint[0]);
  box_index[0] = 1;
  layer_box[0][0] = box_info(&mypoint[left], &mypoint[right], &mypoint[0], r, 0);
  for(i=1;i<n_p+1;i++)
  {
    counter = 0;
    left = -1;
    right = -1;
    left_best = MAXDOUBLE; right_best = MAXDOUBLE;
    if(i==n_p){ // case 2
      counter = getLRlist(inner, &left, &right, (point_avl*)points.root, counter, INFINITY, INFINITY, &left_best, &right_best);
      box_index[i] = counter + 1;
      for(j=0;j<counter+1;j++){
          if(j==0){
            temp_point.x = INFINITY;
            temp_point.y = INFINITY;
            temp_point.z = r[2];
            layer_box[i][j] = box_info(&mypoint[inner[j]], &mypoint[right], &temp_point, r, i);
          }
          else if(j==counter){
            temp_point.x = mypoint[inner[j-1]].x;
            temp_point.y = INFINITY;
            temp_point.z = r[2];
            layer_box[i][j] = box_info(&mypoint[left], &mypoint[inner[j-1]], &temp_point, r, i);}
          else{
            temp_point.x = mypoint[inner[j-1]].x;
            temp_point.y = INFINITY;
            temp_point.z = r[2];
            layer_box[i][j] = box_info(&mypoint[inner[j]], &mypoint[inner[j-1]], &temp_point, r, i);}
        }
    }
    else{ // case 3
      counter = getLRlist(inner, &left, &right, (point_avl*)points.root, counter, mypoint[i].y, mypoint[i].x, &left_best, &right_best);
      box_index[i] = counter + 1;
      if(counter==0){
          layer_box[i][0] = box_info(&mypoint[left], &mypoint[right], &mypoint[i], r, i);
      }
      // CASE: at least one point needs to be removed
      else{
        for(j=0;j<counter+1;j++){
          if(j==0)
            layer_box[i][j] = box_info(&mypoint[inner[j]], &mypoint[right], &mypoint[i], r, i);
          else if(j==counter){
            temp_point.x = mypoint[inner[j-1]].x;
            temp_point.y = mypoint[i].y;
            temp_point.z = mypoint[i].z;
            layer_box[i][j] = box_info(&mypoint[left], &mypoint[inner[j-1]], &temp_point, r, i);}
          else{
            temp_point.x = mypoint[inner[j-1]].x;
            temp_point.y = mypoint[i].y;
            temp_point.z = mypoint[i].z;
            layer_box[i][j] = box_info(&mypoint[inner[j]], &mypoint[inner[j-1]], &temp_point, r, i);}
        }
      }
      // Remove unused points
      for(j=0;j<counter;j++)
      	avl_remove(&points,(avl*)&mypoint[inner[j]]);
      // Insert the new point
      avl_insert(&points,(avl*)&mypoint[i]);
      //listPointTree((point_avl*)points.root,0);

    }
  }
  
  // Calculation for EHVI for each hyper-box
  for(i=0;i<n_p+1;i++){
    for(j=0;j<box_index[i];j++){
      temp_answer = calculation_box(layer_box[i][j], pdf, r, i, n_p);
      for(int k=0; k<pdf.size();k++){
        answer[k] = answer[k] + temp_answer[k];}
    }
  }
  
  return answer;
}


