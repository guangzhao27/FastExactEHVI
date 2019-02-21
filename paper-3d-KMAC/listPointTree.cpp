#include <stdio.h>
#include "avl.h"
#include "listPointTree.h"
#include "hycon3d.h"

void listPointTree(point_avl* a,int m){
	int n=m;
	if(a==0) return;
	if(a->right) listPointTree((point_avl*)a->right,m+1);
	while(n--) printf("   ");
	printf("x:%.3f y:%.3f z:%.3f i:%d\n",a->x, a->y, a->z, a->myindex);
	if(a->left) listPointTree((point_avl*)a->left,m+1);
}


