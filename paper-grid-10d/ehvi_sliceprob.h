//Include this to use the O(n^3) slice-update scheme for calculating the EHVI.
#include <deque>

double ehvi3d_sliceprob(deque<individual*> P, double r[], double mu[], double s[], int dimension);
