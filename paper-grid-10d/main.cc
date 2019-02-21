//Command-line application for calculating the 3-Dimensional
//Expected Hypervolume Improvement.
//By Iris Hupkens, 2013
#include <deque>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "ehvi_calculations.h"
#include "ehvi_sliceprob.h"
#include "ehvi_montecarlo.h"
#include "string.h"
#include "ehvi_multi.h"
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <sys/resource.h>
#include "ehvi_sliceupdate.h"
#include <sstream>
#include <algorithm>
#include <stdlib.h>

using namespace std;

int dimension;

double timeDif(
    struct rusage before,
    struct rusage after)
{
    return ((after.ru_utime.tv_sec * 1000000 + after.ru_utime.tv_usec) -
        (before.ru_utime.tv_sec * 1000000 + before.ru_utime.tv_usec))*1e-6;
}



//Checks if p dominates P. Removes points dominated by p from P and return the number of points removed.
int checkdominance(deque<individual*> & P, individual* p){
  int nr = 0;
  for (int i=P.size()-1;i>=0;i--){
    if (p->f[0] >= P[i]->f[0] && p->f[1] >= P[i]->f[1] && p->f[2] >= P[i]->f[2]){
      cerr << "Individual " << (i+1) << " is dominated or the same as another point; removing." << endl;
      P.erase(P.begin()+i);
      nr++;
    }
  }
  return nr;
}

typedef struct
{
    int nFronts;
    deque<individual*> fronts;
} FILECONTENTS;


//Loads a testcase from the file with the name filename.
void loadtestcase(char *filename, deque<individual*> & testcase, double r[], vector<mus*> & pdf, vector<int*> &num_array){
//void loadtestcase(char *filename, vector<deque<individual*>*> & casearray, double r[], vector<mus*> & pdf, vector<int*> &num_array){
    string line, token;
    int n = 0, casenum = 0;
    int *m;
    ifstream file2;
    file2.open(filename, ios::in);
    if (file2.is_open()){
        cout<<"Reading from the file"<<endl;
        while(getline(file2, line)){
            if (line == "#"){
                
                    //casearray.push_back(testcase);
                    //testcase = *new deque<individual*>;
                    m = new int;
                    *m = casenum;
                    num_array.push_back(m);
                    //casenum = 0;
                
            }
            else{
                casenum++;
                //cout<<line<<endl;
                stringstream ss(line);
                individual * tempvidual = new individual;
                
                string tok;
                
                /*while(ss.rdbuf()->in_avail()){
                    k = ss.rdbuf()->in_avail();
                    ss>>tempvidual->f[t];
                    cout<<tempvidual->f[t]<<endl;
                    t++;
                }*/
                n = 0;
                while(getline(ss, tok, ' ')){
                    //cout<<tok<<endl;
                    stringstream ss(tok);
                    ss>>tempvidual->f[n];
                    //cout<<tempvidual->f[t]<<endl;
                    n++;
                }
                //checkdominance(testcase,tempvidual);
                testcase.push_back(tempvidual);
            }
            /*getline(file2, line);
            cout<<line<<endl;
            stringstream ss(line);
            while(getline(ss, line, ' ')){
                cout<<line<<endl;
            }*/
        }
    }
    else cout<<"Unable to open file";
    dimension = n;
    for(int i = 0; i<n;i++){
        r[i] = 0;
        pdf[0]->mu[i] = 10;
        pdf[0]->s[i] = 2.5;
    }
    file2.close();
    
    
    /*ifstream file;
    int n, inds = 0;
    
    file.open(filename, ios::in);
    
    file >> n;
    n = 5;
    for (int i=0;i<n;i++){
        individual * tempvidual = new individual;
        file >> tempvidual->f[0] >> tempvidual->f[1] >> tempvidual->f[2];
        checkdominance(testcase,tempvidual);
        testcase.push_back(tempvidual);
    }
    file >> r[0] >> r[1] >> r[2];
    while (!file.eof()){
        if (inds > 0)
            pdf.push_back(new mus);
        file >> pdf[inds]->mu[0] >> pdf[inds]->mu[1] >> pdf[inds]->mu[2];
        file >> pdf[inds]->s[0] >> pdf[inds]->s[1] >> pdf[inds]->s[2];
        if (file.fail()){
            //We discover this while trying to read an individual and will end it here.
            pdf.pop_back();
            file.close();
            return;
        }
        inds++;
    }
    file.close();*/
    

}
/*    for(int i = 0; i<3;i++){
        r[i] = 0;
        pdf[0]->mu[i] = 10;
        pdf[0]->s[i] = 2.5;
    }
    //file.close();
}*/

int main(int argc, char *argv[]){
    
    struct rusage ru_before, ru_after;
  
  
  int n = 5;
    vector<int*> num_array;
  deque<individual*> testcase;
    vector<deque<individual*>*> casearray;
  double r[DIMENSIONS];
  double mu[DIMENSIONS];
  double s[DIMENSIONS];
  vector<mus*> pdf;
    
  mus * tempmus = new mus;
  pdf.push_back(tempmus);
  cout << setprecision(10);
    
    cerr << "Loading testcase from file..." << endl;
    loadtestcase(argv[1],testcase,r,pdf, num_array);

    
    
    
    
    cerr << "No scheme specified, defaulting to slice-update..." << endl;
   
    ofstream test2("/Users/guangzhao/Google-Drive/git/2018_exact_EHVI/time/time_exp.txt", ios::app);
    
    test2<<*num_array[1]<<" ";
    
    
    for (int i = 0; i<num_array.size()-1; i++){
        deque<individual*> testcase2 ;
        for (int j = *num_array[i]; j<*num_array[i+1]; j++){
            testcase2.push_back(testcase[j]);
            //cout<<testcase2[i]->f[0]<<endl;
        }
        
        //the better algorithm
        getrusage (RUSAGE_SELF, &ru_before);
        clock_t begin = clock();
        cout << ehvi3d_sliceprob(testcase2, r, pdf[0]->mu, pdf[0]->s, dimension) << endl;
        //ehvi3d_sliceprob(testcase2, r, pdf[0]->mu, pdf[0]->s);
        getrusage (RUSAGE_SELF, &ru_after);
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Elapsed time2: %f seconds\n", timeDif(ru_before, ru_after));
        //printf("Elapsed time2: %f seconds\n", time_spent);
        
        test2<<timeDif(ru_before, ru_after)<<" ";
    }
    test2<<endl;
    
    
    
	
    return 0;
}
