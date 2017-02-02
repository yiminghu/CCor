#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void mccor(double *mat, double *Ans, int *nrow, double *threshold, double *w, double *w1)
{
	int n = *nrow;
	//int m = n - 2;
	double thd = *threshold;
	double w2 = *w1;
	double mi,mj,si,sj,sum;
	//double *Ans = calloc(n*n, sizeof(double));
	//printf("%f",thd);
	for(int i = 0; i < n; ++i){
	   for(int j = i; j < n; ++j){
	   	double sw = 0;
	    for(int k = 0;k < n; ++k){
	    	if(k!=i&&k!=j){
	    		if(fabs(mat[i+k*n])<=thd&&fabs(mat[j+k*n])<=thd){
	    			w[k] = 1;
	    		}else{
	    			if(fabs(mat[i+k*n])>thd&&fabs(mat[j+k*n])>thd){
	    				w[k] = w2*w2;
	    			}else{
	    				w[k] = w2;
	    			}	    			
	    		}
	    	}else{
	    		w[k] = 0;
	    	}
	    	sw += w[k];
	    }
	    //printf("%f",sw);
	    //for(int k = 0;k < n; ++k)
    	//	printf("%f",w[k]);
	   	mi = 0;
	    for(int k = 0;k < n; ++k){
	      if(k!=i&&k!=j){
	        mi += w[k]*mat[i+k*n];
	      }  
	    }
	    mj = 0;
	    for(int k = 0;k < n; ++k){
	      if(k!=i&&k!=j){
	        mj += w[k]*mat[j+k*n];
	      }  
	    }
	    si = 0;
	    for(int k = 0; k < n; ++k){
	      if(k!=i&&k!=j){
	        si += w[k]*mat[i+k*n]*mat[i+k*n];
	      }
	    }
	    sj = 0;
	    for(int k = 0; k < n; ++k){
	      if(k!=i&&k!=j){
	        sj += w[k]*mat[j+k*n]*mat[j+k*n];
	      }
	    }
	    sum = 0;
	    for(int k = 0; k < n; ++k) {
	      if(k!=i&&k!=j){
	        sum += w[k]*mat[i+k*n]*mat[j+k*n];
	      }
	    }
	    Ans[i+j*n] = (sum-mi*mj/sw)/sqrt((si-mi*mi/sw)*(sj-mj*mj/sw));
	    Ans[j+i*n] = Ans[i+j*n];
	   }
    }
}