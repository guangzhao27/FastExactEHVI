#include "ehvi_consts.h"
#include <deque>
#include <algorithm>
#include <iostream> //For error on exception only.
#include <cmath> //INFINITY macro
#include "ehvi_hvol.h"
#include <math.h>

using namespace std;

#ifndef EHVI_SLICEUPDATE
#define EHVI_SLICEUPDATE

//static inline int point_dominate_set(double t[], deque<individual*> P){}

double box_integral(double mu, double sigma, double upper, double lower){
    double integral = 1, chi1, chi2;
    chi1 = chi(upper, mu, sigma);
    if (lower == -INFINITY){
        return chi1;
    }
    chi2 = chi(lower, mu, sigma);
    integral = chi1 - chi2;
    return integral;
}

double box_integral_max(double mu, double sigma, double upper, double lower){
    double chi1, chi2;
    chi1 = chi(-lower, -mu, sigma);
    if (upper == INFINITY){
        return chi1;
    }
    chi2 = chi(-upper, -mu, sigma);
    return chi1 - chi2;
}



double ehvi3d_sliceprob(deque<individual*> P, double r[], double mu[], double s[]){
//EHVI calculation algorithm with time complexity O(n^3).
  double answer = 0, cellength[3]; //The eventual answer.
    
    int d_num = 3;//num of dimension
    int odemeter[d_num], o_bound[d_num];//odometer index initialize odemeter = 0;
    for (int d = 0; d<d_num;d++) odemeter[d] = 0;
    
    
    
  specialind *newind;
  int n = P.size(); //Holds amount of points.
  thingy *Pstruct; //2D array with information about the shape of the dominated hypervolume
  deque<specialind*> o_value[d_num]; //P sorted by x/y/z coordinate with extra information.    //define o_value here

    for (int d = 0; d<d_num; d++) o_bound[d] = n;
  try{
    //Create sorted arrays which also contain extra information allowing the location in
    //the other sorting orders to be ascertained in O(1).
    sort(P.begin(), P.end(), xcomparator);
    for (unsigned int i=0;i<n;i++){
      newind = new specialind;
      newind->point = P[i];
      newind->xorder = i;
      o_value[0].push_back(newind);
      o_value[1].push_back(newind);
      o_value[2].push_back(newind);
    }
    sort(o_value[1].begin(), o_value[1].end(), specialycomparator);
    for (unsigned int i=0;i<n;i++){
      o_value[1][i]->yorder = i;
    }
    sort(o_value[2].begin(), o_value[2].end(), specialzcomparator);
    for (unsigned int i=0;i<n;i++){
      o_value[2][i]->zorder = i;
    }
    //Then also reserve memory for the structure array.
      
  }
  catch (...) {
    cout << "An exception was thrown. There probably isn't enough memory available." << endl;
    cout << "-1 will be returned." << endl;
    return -1;
  }
  //Now we establish dominance in the 2-dimensional slices. Note: it is assumed that
  //P is mutually nondominated. This implementation of that step is O(n^3).
    int dominate_zvalue[n][n];
    for (int j=0; j<n; j++){
        for (int k=0; k<n; k++){
            dominate_zvalue[j][k] = -1;
        }
    }
    for (int i=0;i<n;i++){
        for (int j=o_value[2][i]->yorder;j>=0;j--)
            for (int k=o_value[2][i]->xorder;k>=0;k--){
                dominate_zvalue[k][j] = i;
          //if p dominate (z, x, y) then it also dominate (z, x-1, y-1)
            }
    /*for (int j=o_value[0][i]->zorder;j>=0;j--)
      for (int k=o_value[0][i]->yorder;k>=0;k--){
        Pstruct[k+j*n].xlim = o_value[0][i]->point->f[0] - r[0];
      }
    for (int j=o_value[1][i]->zorder;j>=0;j--)
      for (int k=o_value[1][i]->xorder;k>=0;k--){
        Pstruct[k+j*n].ylim = o_value[1][i]->point->f[1] - r[1];
      }*/
    }
    

    int digit_idx = d_num-1;
    double cl[d_num], cu[d_num];
    while (digit_idx>=0) {
        
        /*for (int d=0;d<d_num;d++){
            cl[d] = (odemeter[d]==0? r[d] :o_value[d][odemeter[d]-1]->point->f[d]);
            cu[d] = (odemeter[d]==n? INFINITY :o_value[d][odemeter[d]]->point->f[d]);
            cellength[d] = cu[d] - cl[d];
        }*/
        cl[0] = (odemeter[0] == 0 ? r[0] : o_value[0][odemeter[0]-1]->point->f[0]);
        cl[1] = (odemeter[1] == 0 ? r[1] : o_value[1][odemeter[1]-1]->point->f[1]);
        cl[2] = (odemeter[2] == 0 ? r[2] : o_value[2][odemeter[2]-1]->point->f[2]);
        cu[0] = (odemeter[0] == n ? INFINITY : o_value[0][odemeter[0]]->point->f[0]);
        cu[1] = (odemeter[1] == n ? INFINITY : o_value[1][odemeter[1]]->point->f[1]);
        cu[2] = (odemeter[2] == n ? INFINITY : o_value[2][odemeter[2]]->point->f[2]);
        cellength[0] = cu[0] - cl[0];
        cellength[1] = cu[1] - cl[1];
        cellength[2] = cu[2] - cl[2];

        if (cellength[0] >0 && cellength[1] >0 && cellength[2]> 0){
            //if (odemeter[0] == n || odemeter[1] == n || Pstruct[odemeter[0]+odemeter[1]*n].highestdominator < odemeter[2]){
            if (odemeter[0] == n || odemeter[1] == n || dominate_zvalue[odemeter[0]][odemeter[1]]< odemeter[2]){
                
                double sum = 1;
                for (int i = 0; i<3; i++){
                    sum *= box_integral_max(mu[i], s[i], cu[i], cl[i]);
                }
                if (sum > 0)
                    answer += sum;
            }
        }
        
        digit_idx = d_num - 1;
        odemeter[digit_idx]++;
        while (odemeter[digit_idx]>o_bound[digit_idx]) {
            odemeter[digit_idx] = 0;
            digit_idx--;
            if (digit_idx >= 0){
                odemeter[digit_idx]++;
            }
        }
        
    }

  return answer;
}
#endif

