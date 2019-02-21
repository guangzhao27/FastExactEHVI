//This is the implementation of variants of the 5term and slice-update scheme in which
//a vector of Gaussian PDFs is used instead of just one.
#include <vector>
#include <deque>
#include <algorithm>
#include <math.h>
#include "ehvi_hvol.h"
#include "ehvi_multi.h"
#include <iostream>
#include "newhelper.h"
#include "ehvi_consts.h"
#include <math.h>
double gausscdf1[1001][1001];
double gausscdf2[1001][1001];
double gausscdf3[1001][1001];
char flag1[1001][1001];
char flag2[1001][1001];
char flag3[1001][1001];
double exipsi1[1001][1001];
double exipsi2[1001][1001];
double exipsi3[1001][1001];
char fpsi1[1001][1001];
char fpsi2[1001][1001];
char fpsi3[1001][1001];




using namespace std;

vector<double> ehvi3d_5term(deque<individual*> P, double r[], vector<mus *> & pdf){
//5-term 3-dimensional ehvi calculation scheme. Subtracts 4 quantities off a rectangular volume.
  vector<double> answer; //The eventual answer
  int n = P.size(); //Holds amount of points.
  double Sminus; //Correction term for the integral.
  deque<individual*> Py, Pz; //P sorted by y/z coordinate
  sort(P.begin(), P.end(), ycomparator);
  for (int i=0;i<P.size();i++){
    Py.push_back(P[i]);
  }
  sort(P.begin(), P.end(), zcomparator);
  for (unsigned int i=0;i<P.size();i++){
    Pz.push_back(P[i]);
  }
  sort(P.begin(), P.end(), xcomparator);
  for (int i=0;i<pdf.size();i++){
    answer.push_back(0);
  }
  for (int z=0;z<=n;z++){
    for (int y=0;y<=n;y++){
      for (int x=0;x<=n;x++){
        double v[DIMENSIONS]; //upper corner of hypervolume improvement box
        double cl[DIMENSIONS], cu[DIMENSIONS]; //Boundaries of grid cells
        cl[0] = (x == 0 ? r[0] : P[x-1]->f[0]);
        cl[1] = (y == 0 ? r[1] : Py[y-1]->f[1]);
        cl[2] = (z == 0 ? r[2] : Pz[z-1]->f[2]);
        cu[0] = (x == n ? INFINITY : P[x]->f[0]);
        cu[1] = (y == n ? INFINITY : Py[y]->f[1]);
        cu[2] = (z == n ? INFINITY : Pz[z]->f[2]);
        //We have to find v. This naive implementation is O(n) per iteration.
        v[0] = r[0];
        v[1] = r[1];
        v[2] = r[2];
        bool dominated = false;
        for (unsigned int i=0;i<P.size();i++){
          if (P[i]->f[0] >= cu[0] && P[i]->f[1] >= cu[1] && P[i]->f[2] >= cu[2]){
            dominated = true;
            break;
          }
          else if (P[i]->f[0] <= cu[0] && P[i]->f[1] >= cu[1] && P[i]->f[2] >= cu[2]){
            if (P[i]->f[0] > v[0])
              v[0] = P[i]->f[0];
          }
          else if (P[i]->f[0] >= cu[0] && P[i]->f[1] <= cu[1] && P[i]->f[2] >= cu[2]){
            if (P[i]->f[1] > v[1])
              v[1] = P[i]->f[1];
          }
          else if (P[i]->f[0] >= cu[0] && P[i]->f[1] >= cu[1] && P[i]->f[2] <= cu[2]){
            if (P[i]->f[2] > v[2])
              v[2] = P[i]->f[2];
          }
        }
        if (dominated)
          continue; //Cell's contribution is 0.
        Sminus = hvol3d(Pz, v, cl);
        double xslice = calculateslice(P, v, cl, 0);
        double yslice = calculateslice(Py, v, cl, 1);
        double zslice = calculateslice(Pz, v, cl, 2);
      //And then we integrate.
        for (int i=0;i<pdf.size();i++){
            double psi1 = exipsi(v[0],cl[0],pdf[i]->mu[0],pdf[i]->s[0]) - exipsi(v[0],cu[0],pdf[i]->mu[0],pdf[i]->s[0]);
            double psi2 = exipsi(v[1],cl[1],pdf[i]->mu[1],pdf[i]->s[1]) - exipsi(v[1],cu[1],pdf[i]->mu[1],pdf[i]->s[1]);
            double psi3 = exipsi(v[2],cl[2],pdf[i]->mu[2],pdf[i]->s[2]) - exipsi(v[2],cu[2],pdf[i]->mu[2],pdf[i]->s[2]);

            double gausscdf1 = gausscdf((cu[0]-pdf[i]->mu[0])/pdf[i]->s[0]) - gausscdf((cl[0]-pdf[i]->mu[0])/pdf[i]->s[0]);
            double gausscdf2 = gausscdf((cu[1]-pdf[i]->mu[1])/pdf[i]->s[1]) - gausscdf((cl[1]-pdf[i]->mu[1])/pdf[i]->s[1]);
            double gausscdf3 = gausscdf((cu[2]-pdf[i]->mu[2])/pdf[i]->s[2]) - gausscdf((cl[2]-pdf[i]->mu[2])/pdf[i]->s[2]);
            double sum = (psi1*psi2*psi3) - (Sminus*gausscdf1*gausscdf2*gausscdf3);
            //gausscdf represents chance of a point falling within the range [cl,cu)
            //psi = partial expected improvement
            //so psi - (gausscdf * (cl - v)) = p's expected distance from cl
            sum -= (xslice * gausscdf2 * gausscdf3 * (psi1 - (gausscdf1 * (cl[0]-v[0]))));
            sum -= (yslice * gausscdf1 * gausscdf3 * (psi2 - (gausscdf2 * (cl[1]-v[1]))));
            sum -= (zslice * gausscdf1 * gausscdf2 * (psi3 - (gausscdf3 * (cl[2]-v[2]))));
            if (sum > 0)
              answer[i] += sum;
        }
      }
    }
  }
  return answer;
}


 
