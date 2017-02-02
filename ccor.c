#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void ccor(double *mat, double *Ans, int *nrow)
{
	int n = *nrow;
	int m = n - 2;
	double mi,mj,si,sj,sum;
	//double *Ans = calloc(n*n, sizeof(double));
	for(int i = 0; i < n; ++i){
	   for(int j = i; j < n; ++j){
	     
	   	mi = 0;
	    for(int k = 0;k < n; ++k){
	      if(k!=i&&k!=j){
	        mi += mat[i+k*n];
	      }  
	    }
	    mj = 0;
	    for(int k = 0;k < n; ++k){
	      if(k!=i&&k!=j){
	        mj += mat[j+k*n];
	      }  
	    }
	    si = 0;
	    for(int k = 0; k < n; ++k){
	      if(k!=i&&k!=j){
	        si += mat[i+k*n]*mat[i+k*n];
	      }
	    }
	    sj = 0;
	    for(int k = 0; k < n; ++k){
	      if(k!=i&&k!=j){
	        sj += mat[j+k*n]*mat[j+k*n];
	      }
	    }
	    sum = 0;
	    for(int k = 0; k < n; ++k) {
	      if(k!=i&&k!=j){
	        sum += mat[i+k*n]*mat[j+k*n];
	      }
	    }
	    Ans[i+j*n] = (sum-mi*mj/m)/sqrt((si-mi*mi/m)*(sj-mj*mj/m));
	    Ans[j+i*n] = Ans[i+j*n];
	   }
    }
}