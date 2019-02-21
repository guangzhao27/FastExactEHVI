   README: 'nlogn_3d V.1.6', MEX code version

 * Compute truncated expected hypervolume improvement in 3-D in time O(n log n). Described in: 
 * Yang, K., Emmerich, M., Deutz, A., & Fonseca, C. M. (2017, March). Computing 3-D Expected  
 * Hypervolume Improvement and Related Integrals in Asymptotically Optimal Time. In International  
 * Conference on Evolutionary Multi-Criterion Optimization (pp. 685-700). Springer, Cham.
 * 
 * The algorithm considers maximization.
 *
 * Implementation: Kaifeng Yang, Michael T. M. Emmerich
 * Algorithm: Kaifeng Yang and Michael T. M. Emmerich, LIACS, Leiden University
 *  
 * 
 * Released under GNU General Public License (GPL) version 2
 * (C) by Kaifeng Yang, Michael T. M. Emmerich, 
 * k.yang@liacs.leidenuniv.nl
 * m.t.m.emmerich@liacs.leidenuniv.nl
 * Date: Nov. 21, 2015
 *
 * 
 *
 * COMPILE:
 * mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" main.cpp avl.cpp boxlist.cpp cmppoint.cpp ehvi_hvol.cpp ehvi_multi.cpp getLRlist.cpp helper.cpp quicksortPZ.cpp
 *
 * Data format:
 * main(pf,ref,[mu std])
 * 	pf: 	current Pareto front
 * 	ref: 	reference points
 * 	mu: 	mean
 * 	std: 	standard deviation
 *  
 * Example:
 * pf = [1 3 4;4 2 3;2 4 2;3 5 1];
 * ref = [0 0 0];
 * [mu std] = [0.5 0.5 4.5 1 1 1;4.5 0.5 4.5 1 1 1];
 * then the results should be around 0.4631 and 5.8292. You can just run the MATLAB script "example.m" in the current folder.
 *
 * Limitations/remarks to this version:
 *    -  maximization is considered. 
 *    -  all the points in the Pareto front should be non dominated
 *
 ********************************************************************
 * The following limitations maybe still exists, NEED TO BE VERIFY. *
 *    -  coordinates should not exceed 100000000000.0		    *
 *    -  duplicate coordinates are not ext. tested		    *	
 *    -  20000 points is maximum				    *
 ********************************************************************
 * 
 * Includes:
 * ANSI C Library for maintainance of AVL Balanced Trees
 *  G. M. Adelson-Velskij & E. M. Landis
 *  Doklady Akad. Nauk SSSR 146 (1962), 263-266
 *  (C) 2000 Daniel Nagy, Budapest University of Technology and Economics
 * Released under GNU General Public License (GPL) version 2
 *

cd Z:\Google-Drive\git\2018_exact_EHVI\KMAC\3D_EHVI_KMAC\
 
