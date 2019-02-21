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



class Upper
{
public:
    Upper(int min = 0):m_min(min){}
    bool operator() (specialind *A, specialind *B)
    {
        return A->point->f[m_min] < B->point->f[m_min];
    }
private:
    int m_min;
};



double ehvi3d_sliceprob(deque<individual*> P, double *r, double *mu, double *s, int dimension){
//EHVI calculation algorithm with time complexity O(n^3).
  
    
    int d_num = dimension;//sizeof(mu)/sizeof(mu[0]);//num of dimension
    
    double answer = 0, cellength[d_num]; //The eventual answer.
    
    
    
    
  specialind *newind;
  unsigned long n = P.size(); //Holds amount of points, points of fronts.
  deque<specialind*> o_value[d_num]; //P sorted by x/y/z coordinate with extra information.    //define o_value here

  try{
    //Create sorted arrays which also contain extra information allowing the location in
    //the other sorting orders to be ascertained in O(1).
    sort(P.begin(), P.end(), xcomparator);
    for (unsigned int i=0;i<n;i++){
      newind = new specialind;
      newind->point = P[i];
      newind->sorder[0] = i;
        
        for (unsigned int m = 0; m<d_num; m++){
            o_value[m].push_back(newind);   //here we didn't copy the individual points, we just copy an array of address of individual points, so each list o_value[0, 1, 2] shares the same sorder[3], they just store the address in different orders.
        }
      /*o_value[0].push_back(newind);
      o_value[1].push_back(newind);
      o_value[2].push_back(newind);*/
    }
    //sort(o_value[1].begin(), o_value[1].end(), specialycomparator);
      for (unsigned int m = 1; m<d_num; m++){
          sort(o_value[m].begin(), o_value[m].end(), Upper(m));
          for (unsigned int i = 0; i<n; i++){
              o_value[m][i] -> sorder[m] = i;
          }
      }
   /* sort(o_value[1].begin(), o_value[1].end(), Upper(1));
    for (unsigned int i=0;i<n;i++){
      o_value[1][i]->sorder[1] = i;
    }
    sort(o_value[2].begin(), o_value[2].end(), Upper(2));
    for (unsigned int i=0;i<n;i++){
      o_value[2][i]->sorder[2] = i;
    }*/
      /*sort(o_value[3].begin(), o_value[3].end(), specialkcomparator);
      for (unsigned int i=0;i<n;i++){
          o_value[3][i]->korder = i;
      }*/
    
    //Then also reserve memory for the structure array.
      
  }
  catch (...) {
    cout << "An exception was thrown. There probably isn't enough memory available." << endl;
    cout << "-1 will be returned." << endl;
    return -1;
  }
  
    
    
    //Now we establish dominance in the 2-dimensional slices. Note: it is assumed that
  //P is mutually nondominated. This implementation of that step is O(n^3).
    int digit_idx;
    int odemeter[d_num], o_bound[d_num];//odometer index initialize odemeter = 0, odometer bound index the carry in each digit;
    for (int d = 0; d<d_num;d++) odemeter[d] = 0;
    
    
    
    int * dominate_value = new int [long(pow(n, d_num-1))];//dominate value is an array used to store the highest value of final coordiantes for which the cell is dominated.
    for (int j=0; j<pow(n, d_num-1); j++){
        dominate_value[j] = -1;
    }
    /*for (int i=0;i<n;i++){
       // for (int j = o_value[0][i] -> sorder[1]; j>=0; j--)
         //   for (int k = o_value[0][i] -> sorder[2]; k>=0; k--)
        for (int j=o_value[0][i]->sorder[2];j>=0;j--)
            for (int k=o_value[0][i]->sorder[1];k>=0;k--){
                dominate_value[k+j*n] = i;
            }
    }*/
    o_bound[0] = n-1;
    
    int tempidx;
    digit_idx = 0;
    do{
        if (digit_idx == 0 ){
            for (int d = 1; d<d_num; d++){
                o_bound[d] = o_value[0][odemeter[0]]->sorder[d];//for each dimension, the cell with value smaller than sorder[odemeter[d]] is dominated by front.
            }
        }
        
        tempidx = 0;
        for (int m = 1; m<d_num; m++){
            tempidx += odemeter[m]*pow(n, m-1);
        }
        dominate_value[tempidx] = odemeter[0];
        digit_idx = d_num -1;
        odemeter[digit_idx]++;
        
        while (odemeter[digit_idx]>o_bound[digit_idx] && digit_idx > 0) {
            odemeter[digit_idx] = 0;
            digit_idx--;
            odemeter[digit_idx]++;
        }
    }while (odemeter[0]<=o_bound[0]);


    for (int d = 0; d<d_num;d++) odemeter[d] = 0;
    for (int d = 0; d<d_num; d++) o_bound[d] = n;
    double cl[d_num], cu[d_num];
    do{
        bool activeflag = true, flag2 = false;// flag = true stands for an active cell
        tempidx = 0;
        for (int d=0;d<d_num;d++){
            cl[d] = (odemeter[d]==0? r[d] :o_value[d][odemeter[d]-1]->point->f[d]);
            cu[d] = (odemeter[d]==n? INFINITY :o_value[d][odemeter[d]]->point->f[d]);
            cellength[d] = cu[d] - cl[d];
            activeflag = activeflag && (cu[d]-cl[d]>0);//length must be larger than 0
        }
        if (activeflag){
            for (int m = 1; m<d_num; m++){
                flag2 |= (odemeter[m] == n);
                if (flag2) break;
                tempidx += odemeter[m]*pow(n, m-1);
            }
                flag2 = flag2 || (dominate_value[tempidx]<odemeter[0]);
        }
        
        
       // if (cellength[0] >0 && cellength[1] >0 && cellength[2]> 0){
            //if (odemeter[0] == n || odemeter[1] == n || Pstruct[odemeter[0]+odemeter[1]*n].highestdominator < odemeter[2]){
         //   if (odemeter[1] == n || odemeter[2] == n || dominate_value[odemeter[1]+odemeter[2]*n]< odemeter[0]){
            if (activeflag & flag2){
                double sum = 1;
                for (int i = 0; i<d_num; i++){
                    sum *= box_integral_max(mu[i], s[i], cu[i], cl[i]);
                }
                if (sum > 0)
                    answer += sum;
            }
       // }
        
        digit_idx = d_num - 1;
        odemeter[digit_idx]++;
        while (odemeter[digit_idx]>o_bound[digit_idx] && digit_idx > 0) {
            odemeter[digit_idx] = 0;
            digit_idx--;
            odemeter[digit_idx]++;
        }
    }while (odemeter[0]<=o_bound[0]);
    
  return answer;
}
#endif

