/*

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------

*/

// To do: 
// - can we sort less often or reduce/optimise dominance checks? 
// - should we use FPL's data structure? 
// - two changes in read.c 
// - heuristics 

// opt:  0 = basic, 1 = sorting, 2 = slicing to 2D, 3 = slicing to 3D 

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "wfg.h"
#include "avl.h"
#include "newhelper.h"
#include <stdlib.h>

#define MAXIMISING true

#if MAXIMISING
#define BEATS(x,y)   (x >  y) 
#define BEATSEQ(x,y) (x >= y) 
#else
#define BEATS(x,y)   (x <  y) 
#define BEATSEQ(x,y) (x <= y) 
#endif

#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y)) 
#define BETTER(x,y)  (BEATS(y,x) ? (y) : (x))
#define MU 10
#define SIGMA 2.5
#define REF 0

double box_integral(double mu, double sigma, double upper, double lower){
	double integral = 1, chi1, chi2;
	chi1 = chi(upper, mu, sigma);
  if (lower == -INFINITY) return chi1;
  chi2 = chi(lower, mu, sigma);
  integral = chi1 - chi2;
  return integral;
}



double g_function(double l, double u, double a, double b, double mu, double sigma){
    // a function used to caculate dimension factor of cell contribution, i.e. cellcontribution = \prod g_function_j, j<= m, refer eq 16 in "Fast calculation of multiobjective probability of improvement and expected improvement criteria for Pareto optimization"
  
  double dimension_factor = 0;
  dimension_factor = BETTER(1, 2);
  if (a >= u){
    dimension_factor = (b-a)*(gausscdf((u-mu)/sigma)-gausscdf((l-mu)/sigma));
    return dimension_factor;
  }
  if (b > l){
    double normal1, normal2;
    
    normal2 = (WORSE(b, u)-mu)/sigma;
    if (a == -INFINITY){
      dimension_factor = (b-mu)*(gausscdf(normal2)-gausscdf((l-mu)/sigma))+\
      sigma*(gausspdf(normal2)-gausspdf((l-mu)/sigma));
      return dimension_factor;
    }
    normal1 = (BETTER(a, l)-mu)/sigma;
    dimension_factor = (b-a)*(gausscdf(normal1)-gausscdf((l-mu)/sigma))+\
    (b-mu)*(gausscdf(normal2)-gausscdf(normal1) )+\
    sigma*(gausspdf(normal2)-gausspdf(normal1));
    return dimension_factor;
  }
  
  return dimension_factor;
}







POINT mu, sigma;


int n;     // the number of objectives 
POINT ref; // the reference point 

FRONT *fs;      // memory management stuff 
int fr = 0;     // current depth 
int frmax = -1; // max depth malloced so far (for opt = 0) 
int maxm = 0;   // identify the biggest fronts in the file 
int maxn = 0;


double hv(FRONT);
static avl_tree_t *tree;

static int compare_tree_asc( const void *p1, const void *p2)
{
  const double x1= *((const double *)p1+1);
  const double x2= *((const double *)p2+1);
  
  if (x1 != x2) return (x1 > x2) ? -1 : 1;
  else          return 0;
}



int greater(const void *v1, const void *v2)
// this sorts points improving in the last objective
{
  POINT p = *(POINT*)v1;
  POINT q = *(POINT*)v2;
  #if opt == 1
  for (int i = n - fr - 1; i >= 0; i--)
  #else
    for (int i = n - 1; i >= 0; i--)
  #endif
      if BEATS(p.objectives[i],q.objectives[i]) return  1;
    else
      if BEATS(q.objectives[i],p.objectives[i]) return -1;
    return 0;
  }


  int dominates2way(POINT p, POINT q)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise 
  {
  // domination could be checked in either order 
  #if opt == 1
    for (int i = n - fr - 1; i >= 0; i--)
  #else
      for (int i = n - 1; i >= 0; i--)
  #endif
        if BEATS(p.objectives[i],q.objectives[i]) 
          {for (int j = i - 1; j >= 0; j--) 
           if BEATS(q.objectives[j],p.objectives[j]) return 0; 
           return -1;}
           else
            if BEATS(q.objectives[i],p.objectives[i]) 
              {for (int j = i - 1; j >= 0; j--) 
               if BEATS(p.objectives[j],q.objectives[j]) return 0; 
               return  1;}
               return 2;
             }


             void makeDominatedBit(FRONT ps, int p)
// creates the front ps[p+1 ..] in fs[fr], with each point bounded by ps[p] and dominated points removed 
             {
  // when opt = 0 each new frame is allocated as needed, because the worst-case needs #frames = #points 
  #if opt == 0
              if (fr > frmax)
                {frmax = fr;
                 fs[fr].points = malloc(sizeof(POINT) * maxm);
                 for (int j = 0; j < maxm; j++) 
                 {
                   fs[fr].points[j].objectives = malloc(sizeof(OBJECTIVE) * maxn);
                 }
               }
  #endif

               int z = ps.nPoints - 1 - p;
               for (int i = 0; i < z; i++)
                for (int j = 0; j < n; j++) 
                  fs[fr].points[i].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[p + 1 + i].objectives[j]); 
                POINT t;
                fs[fr].nPoints = 1;
                for (int i = 1; i < z; i++)
                  {int j = 0;
                   bool keep = true;
                   while (j < fs[fr].nPoints && keep)
                     switch (dominates2way(fs[fr].points[i], fs[fr].points[j]))
                   {case -1: t = fs[fr].points[j];
                     fs[fr].nPoints--; 
                     fs[fr].points[j] = fs[fr].points[fs[fr].nPoints]; 
                     fs[fr].points[fs[fr].nPoints] = t; 
                     break;
                     case  0: j++; break;
          // case  2: printf("Identical points!\n");
                     default: keep = false;
                   }
                   if (keep) {t = fs[fr].points[fs[fr].nPoints]; 
                    fs[fr].points[fs[fr].nPoints] = fs[fr].points[i]; 
                    fs[fr].points[i] = t; 
                    fs[fr].nPoints++;}
                  }
                  fr++;
                }



                double inclhv(POINT p)
// returns the inclusive hypervolume of p
                {
                  double volume = 1;
                  for (int i = 0; i < n; i++) {
    //volume *= fabs(p.objectives[i] - ref.objectives[i]);
                   volume *= box_integral(-(mu.objectives[i]), sigma.objectives[i], -ref.objectives[i], -p.objectives[i]);
                 }
                 return volume;
               }

               double exclhv(FRONT ps, int p)
// returns the exclusive hypervolume of ps[p] relative to ps[p+1 ..] 
               {
                double volume = inclhv(ps.points[p]);
                if (ps.nPoints > p + 1) 
                {
                 makeDominatedBit(ps, p);
                 volume -= hv(fs[fr - 1]);
                 fr--;
               }
               return volume;
             }

             double hv(FRONT ps){
// returns the ehvi of ps[0 ..]
              qsort(ps.points, ps.nPoints, sizeof(POINT), greater);
              double volume = 0;
              for (int i = 0; i < ps.nPoints; i++) volume += exclhv(ps, i);
                return volume;
            }

            double ehvi_wfg(FRONT ps){
	//ehvi equal to the integral over region = (half region) - (dominated region)
             double halfregion = 1;
             for(int i = 0; i<n; i++){
              halfregion *= chi(-ref.objectives[i], -mu.objectives[i], sigma.objectives[i]);
            }
            return halfregion - hv(ps);
          }






          double hv3_AVL(FRONT ps)
/* hv3_AVL: 3D algorithm code taken from version hv-1.2 available at
 http://iridia.ulb.ac.be/~manuel/hypervolume and proposed by:
 
 Carlos M. Fonseca, Luís Paquete, and Manuel López-Ibáñez.  An improved
 dimension-sweep algorithm for the hypervolume indicator. In IEEE
 Congress on Evolutionary Computation, pages 1157-1163, Vancouver,
 Canada, July 2006.
 
 Copyright (c) 2009
 Carlos M. Fonseca <cmfonsec@ualg.pt>
 Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
 Luis Paquete <paquete@dei.uc.pt>
 */

// returns the hypervolume of ps[0 ..] in 3D
// assumes that ps is sorted improving
          {
    qsort(ps.points, ps.nPoints, sizeof(POINT), greater);//sort from small to large, objective is minimizing
    avl_init_node(ps.points[0].tnode,ps.points[0].objectives);
    avl_insert_top(tree,ps.points[0].tnode);
    
    //double hypera = (ref.objectives[0] - ps.points[0].objectives[0]) *
    //(ref.objectives[1] - ps.points[0].objectives[1]);
    
    double integral0 = box_integral(-(mu.objectives[0]), sigma.objectives[0], ref.objectives[0], ps.points[0].objectives[0]), integral1 =box_integral(-(mu.objectives[1]), sigma.objectives[1], ref.objectives[1], ps.points[0].objectives[1]);
    
    double hypera = integral0*integral1;
    
    double height;
    if (ps.nPoints == 1)
      height = box_integral(-(mu.objectives[2]), sigma.objectives[2], ref.objectives[2], ps.points[0].objectives[2]);
      //  height = ref.objectives[2] - ps.points[0].objectives[2];
    // height = -(ref.objectives[2] - ps.points[ps.nPoints-1].objectives[2]);
    else
      height = box_integral(-(mu.objectives[2]), sigma.objectives[2], ps.points[1].objectives[2], ps.points[0].objectives[2]);
    //    height = ps.points[1].objectives[2] - ps.points[0].objectives[2];
    //    height = ps.points[ps.nPoints-2].objectives[2] - ps.points[ps.nPoints-1].objectives[2];
    
    double hyperv = hypera * height;
    
    for (int i = 1; i < ps.nPoints; i++)
    {
      if (i == ps.nPoints - 1)
            //height = ref.objectives[2] - ps.points[i].objectives[2];
        height = box_integral(-(mu.objectives[2]), sigma.objectives[2], ref.objectives[2], ps.points[i].objectives[2]);
      else
            //height = ps.points[i-1].objectives[2] - ps.points[i].objectives[2];
        height = box_integral(-(mu.objectives[2]), sigma.objectives[2], ps.points[i+1].objectives[2], ps.points[i].objectives[2]);
            //height = ps.points[i+1].objectives[2]- ps.points[i].objectives[2];
      
        // search tree for point q to the right of current point
      const double * prv_ip, * nxt_ip;
      avl_node_t *tnode;
      
      avl_init_node(ps.points[i].tnode, ps.points[i].objectives);
      
      if (avl_search_closest(tree, ps.points[i].objectives, &tnode) <= 0) {
        nxt_ip = (double *)(tnode->item);
        tnode = tnode->prev;
      } else {
        nxt_ip = (tnode->next!=NULL)
        ? (double *)(tnode->next->item)
        : ref.objectives;
      }
        // if p is not dominated
      if (nxt_ip[0] > ps.points[i].objectives[0]) {
        
            // insert p in tree
        avl_insert_after(tree, tnode, ps.points[i].tnode);
        
        if (tnode !=NULL) {
          prv_ip = (double *)(tnode->item);
          
          if (prv_ip[0] > ps.points[i].objectives[0]) {
            const double * cur_ip;
            
            tnode = ps.points[i].tnode->prev;
                    // cur_ip = point dominated by pp with highest [0]-coordinate
            cur_ip = (double *)(tnode->item);
            
                    // for each point in s in tree dominated by p
            while (tnode->prev) {
              prv_ip = (double *)(tnode->prev->item);
                        // decrease area by contribution of s
                        //hypera -= (prv_ip[1] - cur_ip[1])*(nxt_ip[0] - cur_ip[0]);
              integral0 = box_integral(-(mu.objectives[0]), sigma.objectives[0], nxt_ip[0], cur_ip[0]);
              integral1 =box_integral(-(mu.objectives[1]), sigma.objectives[1], prv_ip[1], cur_ip[1]);
              hypera -= integral0*integral1;
              if (prv_ip[0] < ps.points[i].objectives[0])
                            break; // prv is not dominated by pp
                          cur_ip = prv_ip;
                        // remove s from tree
                          avl_unlink_node(tree,tnode);
                          tnode = tnode->prev;
                        }
                        
                    // remove s from tree
                        avl_unlink_node(tree,tnode);
                        
                        if (!tnode->prev) {
                        // decrease area by contribution of s
                        //hypera -= (ref.objectives[1] - cur_ip[1])*(nxt_ip[0] - cur_ip[0]);
                          integral0 = box_integral(-(mu.objectives[0]), sigma.objectives[0], nxt_ip[0], cur_ip[0]);
                          integral1 =box_integral(-(mu.objectives[1]), sigma.objectives[1], ref.objectives[1], cur_ip[1]);
                          hypera -= integral0*integral1;
                          prv_ip = ref.objectives;
                        }
                      }
                    } else
                    prv_ip = ref.objectives;
                    
            // increase area by contribution of p
            //hypera += (prv_ip[1] -ps.points[i].objectives[1])*(nxt_ip[0] -ps.points[i].objectives[0]);
                    integral0 = box_integral(-(mu.objectives[0]), sigma.objectives[0], nxt_ip[0], ps.points[i].objectives[0]);
                    integral1 =box_integral(-(mu.objectives[1]), sigma.objectives[1], prv_ip[1], ps.points[i].objectives[1]);
                    hypera += integral0*integral1;
                    
                  }
                  
                  if (height > 0)
                    hyperv += hypera * height;
                }
                avl_clear_tree(tree);
                
                
                double halfregion = 1;
                for(int i = 0; i<n; i++){
                  halfregion *= chi(-ref.objectives[i], -mu.objectives[i], sigma.objectives[i]);
                }
                return halfregion - hyperv;
    //return hyperv;
              }






//CDD13 algorithm modified WFG
              typedef struct{
                POINT upper;
                POINT lower;
                int sign;
              }BOX;

              typedef struct{
                int nBoxes;
                int ndim;
                BOX *boxesp;
              } SETBOX;

              void boxdecomp(FRONT ps, SETBOX *SetBp);

              void initialbox(SETBOX *SetBp){
    //just add a new box to Setbox and set the sign and space for upper and lower, not even set value for upper and lower bound
                int boxnum = 0;
                boxnum = SetBp->nBoxes;
                SetBp->nBoxes++;
                SetBp->boxesp = realloc(SetBp->boxesp, sizeof(BOX)*SetBp->nBoxes);
    //SetBp->boxesp[boxnum].upper = ps.points[p];
                BOX *f = &SetBp->boxesp[boxnum];
    //f->upper.objectives = ps.points[p].objectives;
                f->upper.objectives = malloc(sizeof(double)*n);
                f->lower.objectives = malloc(sizeof(double)*n);
                SetBp->boxesp[boxnum].sign = 1;
                return;
              }

              void exclbox(FRONT ps, int p, SETBOX *SetBp)
// returns the exclusive hypervolume of ps[p] relative to ps[p+1 ..]
              {
                int boxnum = 0;
                
    //add a box here, the box should include a multiple dimension upper and lower bound
                boxnum = SetBp->nBoxes;
                initialbox(SetBp);
                memcpy(SetBp->boxesp[boxnum].upper.objectives, ps.points[p].objectives, sizeof(double)*n);
                memcpy(SetBp->boxesp[boxnum].lower.objectives, ref.objectives, sizeof(double)*n);
                if (ps.nPoints > p + 1)
                {
                  makeDominatedBit(ps, p);
                  boxdecomp(fs[fr - 1], SetBp);
        //boxset = boxset + hv(fs[fr-1])
                  for (int i = boxnum+1; i<SetBp->nBoxes; i++){
                    SetBp->boxesp[i].sign *= -1;
                  }
                  fr--;
                }
                
                return;
              }


              void boxdecomp(FRONT ps, SETBOX *SetBp){
    // returns the decomposition box set
                qsort(ps.points, ps.nPoints, sizeof(POINT), greater);
                
    for (int i = 0; i < ps.nPoints; i++) //boxset = box + exclboxset
      exclbox(ps, i, SetBp);
    
    return;
  }

  double ehvi_cdd(FRONT ps){
    //ehvi equal to the integral over region = (half region) - (dominated region)
    SETBOX *SetBp = malloc(sizeof(SETBOX));
    SetBp->nBoxes = 0;
    SetBp->boxesp = NULL;
    initialbox(SetBp);
    for (int i =0; i<n ; i++){
      SetBp->boxesp[0].upper.objectives[i] = INFINITY;
    }
    memcpy(SetBp->boxesp[0].lower.objectives, ref.objectives, sizeof(double)*n);
    boxdecomp(ps, SetBp);
    for(int i = 1; i<SetBp->nBoxes; i++){
      SetBp->boxesp[i].sign *= -1;
    }
    /*for(int i = 0; i<SetBp->nBoxes;i++){
        for(int j = 0; j<n; j++){
            printf("%f ", SetBp->boxesp[i].upper.objectives[j]);
        }
        printf("%d\n", SetBp->boxesp[i].sign);
    }*/
    double ehvi = 0;
    for(int i = 0; i<SetBp->nBoxes; i++){
      for(int k = 0; k<SetBp->nBoxes; k++){
        double boxcontribution = 1;
        double *il, *iu, *kl, *ku;
        int isign, ksign;
        
        il = SetBp->boxesp[i].lower.objectives;
        iu = SetBp->boxesp[i].upper.objectives;
        isign = SetBp->boxesp[i].sign;
        kl = SetBp->boxesp[k].lower.objectives;
        ku = SetBp->boxesp[k].upper.objectives;
        ksign = SetBp->boxesp[k].sign;
        for (int j = 0; j<n; j++){
                //in this example, the objective are maximized
          boxcontribution *= g_function(-iu[j], -il[j], -ku[j], -kl[j], -mu.objectives[j], sigma.objectives[j]);
        }
        ehvi += boxcontribution*isign*ksign;
      }
    }
    printf("%d\n",SetBp->nBoxes);
    return ehvi;
  }











  double timeDif(struct rusage before, struct rusage after){
    return ((after.ru_utime.tv_sec * 1000000 + after.ru_utime.tv_usec) -
      (before.ru_utime.tv_sec * 1000000 + before.ru_utime.tv_usec))*1e-6;
  }
double timeDif2(//suit for any platform
 struct timeval *before,
 struct timeval *after)
{
  return ((after->tv_sec * 1000000 + after->tv_usec) -
    (before->tv_sec * 1000000 + before->tv_usec))*1e-6;
}


int main(int argc, char *argv[]) 
// processes each front from the file
{
  FILECONTENTS *f = readFile(argv[1]);
    //FILECONTENTS *p = readFile(argv[2]);

  // find the biggest fronts
  for (int i = 0; i < f->nFronts; i++)
    {if (f->fronts[i].nPoints > maxm) maxm = f->fronts[i].nPoints;
     if (f->fronts[i].n       > maxn) maxn = f->fronts[i].n;
   }

  // allocate memory
   mu.objectives = malloc(sizeof(OBJECTIVE) * maxn);
   sigma.objectives = malloc(sizeof(OBJECTIVE) * maxn);
   for (int i = 0; i < maxn; i++){
    mu.objectives[i] = MU;
    sigma.objectives[i] = SIGMA;
      //mu.objectives[i] = p->fronts[0].points[0].objectives[i];
  	//sigma.objectives[i] = p->fronts[0].points[0].objectives[i+maxn];
  }

  #if opt == 0
  fs = malloc(sizeof(FRONT) * maxm);
  #else

  // slicing (opt > 1) saves a level of recursion
  int maxd = maxn - (opt / 2 + 1); 
  fs = malloc(sizeof(FRONT) * maxd);

  // 3D base (opt = 3) needs space for the sentinels
  int maxp = maxm + 2 * (opt / 3);
  //int maxp = 100000;
  for (int i = 0; i < maxd; i++) 
    {fs[i].points = malloc(sizeof(POINT) * maxp); 
     for (int j = 0; j < maxp; j++) 
     {
       fs[i].points[j].tnode = malloc(sizeof(avl_node_t));
       // slicing (opt > 1) saves one extra objective at each level
       fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (maxn - (i + 1) * (opt / 2)));
     }
   }
  #endif

   tree = avl_alloc_tree ((avl_compare_t) compare_tree_asc,
     (avl_freeitem_t) free);
   
  //printf("opt is %d\n", opt);
  // initialize the reference point
   ref.objectives = malloc(sizeof(OBJECTIVE) * maxn);
   ref.tnode = malloc(sizeof(avl_node_t));
   
   for (int i = 0; i < maxn; i++) ref.objectives[i] = REF;

    for (int i = 0; i < f->nFronts; i++) 
    {      
      struct timeval tv1, tv2;
      
      struct rusage ru_before, ru_after;
      
      
      
      n = f->fronts[i].n;
      getrusage (RUSAGE_SELF, &ru_before);
      gettimeofday(&tv1, NULL);
      printf("hv(%d) = %1.10f\n", i+1, ehvi_wfg(f->fronts[i]));
      
      
      
      
        //printf("%1.10f\n", hv3_AVL(f->fronts[i]));

      gettimeofday(&tv2, NULL);
      getrusage (RUSAGE_SELF, &ru_after);
      printf("WFG Time: %f(s)\n", timeDif(ru_before, ru_after));
       // printf("Elapsed time2: %f seconds\n", timeDif2(&tv1, &tv2));
      
      
      
      FILE *fp;
      fp = fopen("/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/time/time_wfg.txt", "a");
      if (i == 0){
        fprintf(fp, "%d ", f->fronts[0].nPoints);
      }
      fprintf(fp, "%1.6f ", timeDif(ru_before, ru_after));
      if (i == f->nFronts-1){
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
    
    for (int i = 0; i < f->nFronts; i++)
    {
      struct timeval tv1, tv2;
      struct rusage ru_before, ru_after;
      getrusage (RUSAGE_SELF, &ru_before);
      gettimeofday(&tv1, NULL);
      
      n = f->fronts[i].n;
      
      printf("hv(%d) = %1.10f\n", i+1, ehvi_cdd(f->fronts[i]));
      
      gettimeofday(&tv2, NULL);
      getrusage (RUSAGE_SELF, &ru_after);
      printf("cdd time: %f(s)\n", timeDif(ru_before, ru_after));
      
      FILE *fp;
      fp = fopen("/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/time/time_cdd.txt", "a");
      if (i == 0){
        fprintf(fp, "%d ", f->fronts[0].nPoints);
      }
      fprintf(fp, "%1.6f ", timeDif(ru_before, ru_after));
      if (i == f->nFronts-1){
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
    if (n == 3){
      for (int i = 0; i < f->nFronts; i++)
      {
        struct rusage ru_before, ru_after;
        
            for(int j = 0; j < f->fronts[i].nPoints; j++) {//convert the value to minus to suit for avl-tree algorithm with minimizing objectives
                //printf("%f\n",f->fronts[i].points[j].objectives[2]);
              for(int k = 0; k < f->fronts[i].n; k++) {
                f->fronts[i].points[j].objectives[k] = -f->fronts[i].points[j].objectives[k];
              }
            }
            
            getrusage (RUSAGE_SELF, &ru_before);
            
            printf("hv2(%d) = %1.10f\n",i+1, hv3_AVL(f->fronts[i]));
            getrusage (RUSAGE_SELF, &ru_after);
            printf("avl time: %f(s)\n", timeDif(ru_before, ru_after));
            
            FILE *fp;
            fp = fopen("/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/time/time_avl.txt", "a");
            if (i == 0){
              fprintf(fp, "%d ", f->fronts[0].nPoints);
            }
            fprintf(fp, "%1.6f ", timeDif(ru_before, ru_after));
            if (i == f->nFronts-1){
              fprintf(fp, "\n");
            }
            fclose(fp);
          }
        }

        return 0;
      }
