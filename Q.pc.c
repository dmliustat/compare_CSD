/**********************************************************

 * program to calcalate Q matrices for the 2-allele case * 

 * L-loci, simple version                                *

 *********************************************************/

/*include NAG routines for inverting a matrix*/

#include <stdio.h>
#include "header.h"
#include "common.h"
#include <math.h>
#define ERROR_MIN (0.000000001)
#define TINY 1.0e-20;


extern long L; /*DL: number of loci*/
extern long M;

extern double P[2][2]; /*DL: mutation transition matrix*/



void Qmat (const double theta, double rho);
void matrix_inverse_new(double **a,long n);

void exp_matrix(double **out, double a);

extern double ****Q; /*Q matrices- Q[][][j][l] is for loci l, rate theta/(j+theta-1) */
/*DL: Q is for single-locus model*/ 
extern double *****Qb;/*joint q matrices allele to allele at loci, in time, no of inds*/
/*DL: Qb is for multilocus model*/ 

void Qmat (const  double theta,double rho)
{
  long i,j,l,k,ii,jj;
  double lambda,con, **Id,**temp;
  double s[]={0.3225476896192312,1.74576110115834658,4.53662029692112798,9.39507091230113313};
  
  
  Id = (double **)calloc(2,sizeof(double *));
  temp = (double **)calloc(2,sizeof(double *));
  for (i=0;i<2;i++){
    temp[i]=dou_vec_init(2);
    Id[i]= dou_vec_init(2);
  }		 

  
  /*generate Identity  matrix*/
  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
      Id[i][j]= (i==j);
    }
  }
  
  /* do first entry- matrix of rows of mu */

  for(i=0;i<2;i++){
    for(j=0;j<2;j++){
      for(l=0;l<L;l++) Q[i][j][l][0]= 0.5;
    }
  }
 
  /* do other M entries */
  
  for(i=1;i<(M+1);i++){
    /*DL: when H contains i chromosomes*/
    for(l=0;l<L;l++){
      /*DL: for each locus l*/
      lambda= theta/((double)i+theta);
      
      /* calculate (I-lambda*P) */
      
      for(j=0;j<2;j++){
	      for(k=0;k<2;k++){
	        temp[k][j]=Id[k][j]-lambda*P[k][j];
	      }
      }
      
      matrix_inverse_new(temp,2); /* invert - to temp */    
      
       /* update Q - (1-lambda)*temp */
      for(j=0;j<2;j++){
	      for(k=0;k<2;k++){
	        Q[k][j][l][i]=(1.0-lambda)*temp[k][j];
	      }
      }
    }
  }
   
  /*calculate Qb matrix*/ /*DL: Qb is defined in Appendix A, Q(t)=exp{theta*t*(P-I)/L}. Note that *(theta+l) is equivalent to theta/L*/
  for(i=0;i<(NODE_MAX);i++){
    for(j=0;j<4;j++){
      for(l=0;l<L;l++){
	      lambda = (i+1+rho); /*should include a rho-which one (all different possibilities?!)*/
        /*DL: lambda should be the number of chromosomes contained in the current configuration H?????????*/
	      /* calculate exp_matrix( - theta s / lambda *(I- P)) */
	      exp_matrix(temp, theta*s[j]/lambda); /*DL: temp=exp{(*(theta+l)*s[j]/lambda)*P}*/
	      con=exp(-1* (theta*s[j]/lambda));
	      for(ii=0;ii<2;ii++){
	        for(jj=0;jj<2;jj++){
	          Qb[ii][jj][l][j][i]= temp[ii][jj]*con; /*DL: Qb is 2*2*L*4*QB_MAX*/ 
	        }
	      }
	    }
    }
    /*(void)fprintf(stderr,"Calc matrix: %d",i+1);*/
  }
  for(i=0;i<2;i++){
    free(Id[i]);
    free(temp[i]);
  }
  free(Id);
  free(temp);

} 

/**/
void exp_matrix(double **out, double a) /*DL: out=exp(a*P)*/
{
  double **temp,**mult,**temp2, **out_temp,err,sum;
  long i,j,l,ii;

  temp =  (double **)calloc(2,sizeof(double *));
  for (i=0;i<2;i++) temp[i]= dou_vec_init(2);
	

  temp2 =  (double **)calloc(2,sizeof(double *));
  for (i=0;i<2;i++) temp2[i]= dou_vec_init(2);
	

  mult =  (double **)calloc(2,sizeof(double *));
  for (i=0;i<2;i++) mult[i]= dou_vec_init(2);
	

  out_temp  = (double **)calloc(2,sizeof(double *));
  for (i=0;i<2;i++) out_temp[i]= dou_vec_init(2);
       
 
  /* initial temp is the identity */
  for (i=0;i<2;i++){
    for (j=0;j<2;j++){
      temp[i][j]= (i==j);
      out_temp[i][j]= temp[i][j];
    }
  }

  /* mult is the matrix in the exponent */

  for (i=0;i<2;i++){
    for (j=0;j<2;j++){
      mult[i][j]=a* P[i][j];
      /*mult[i][j]*= -a;*/
    }
  }

  /*approx exponent*/

  err=1;
  for(l=1;l<25 * (10) && err>(ERROR_MIN);l++){
    /* calculate mult*temp/l =mult^l/l!*/
    for (i=0;i<2;i++){
      for (j=0;j<2;j++){
        temp2[i][j]=0;
        for( ii=0;ii<2;ii++){
	        temp2[i][j]+=mult[i][ii]*temp[ii][j]/((double)l);
        }
      }
    }
    err*= a/((double)l);
    /* set temp=temp2 out +=temp2/(l!) */
    for (i=0;i<2;i++){
      for (j=0;j<2;j++){
	      temp[i][j]=temp2[i][j];
	      out_temp[i][j] +=temp2[i][j];
	    }
    } 
  }

  
  for (i=0;i<2;i++){
    for (j=0;j<2;j++){
      out[i][j] =(double)out_temp[i][j];
    }   
  }
  
  
  
  for(i=0;i<2;i++){
    free(temp[i]);
    free(temp2[i]);
    free(out_temp[i]);
    free(mult[i]);
  }
  free(temp);
  free(temp2);
  free(out_temp);
  free(mult);
}

/*inefficient Gauss-Jordan elimination*/
void matrix_inverse_new(double **a,long n)
{
  long i,j,k,imax;
  double big,temp;
  double **I;

  /*set-up I*/
  I = (double **)calloc(n,sizeof(double *));
  for(i=0;i<n;i++) I[i]=dou_vec_init(n);

  for(i=0;i<n;i++) for(j=0;j<n;j++) I[i][j]=0.0;
  for(i=0;i<n;i++) I[i][i]=1.0;

  /*loop over columns*/

  for(j=0;j<n;j++){
    /*firstly find largest a[i][j], i>=j*/
    big=0.0;
    for(i=j;i<n;i++){
      if(fabs(a[i][j])>big){
		big=fabs(a[i][j]);
		imax=i;
      }
    }
    if(big<=0.0){
      (void)fprintf(stderr,"Error in matrix_inverse routine: singular matrix? \n");
      abort();
    }

    /*swap rows*/
    for(k=0;k<n;k++){
      temp=a[j][k];
      a[j][k]=a[imax][k];
      a[imax][k]=temp;
      temp=I[j][k];
      I[j][k]=I[imax][k];
      I[imax][k]=temp;
    }

    /*divide row by a[j][j]*/
    temp=a[j][j];
     for(k=0;k<n;k++){
       a[j][k]/=temp;
       I[j][k]/=temp;
     }

    /*calculate jth column*/
    for(i=0;i<n;i++){
      if(j!=i){
	      temp=a[i][j];
	      for(k=0;k<n;k++){
	        a[i][k]-=temp*a[j][k];
	        I[i][k]-=temp*I[j][k];
	 }
       }
     }
  }

  /*replace A by I*/
  for(i=0;i<n;i++) for(j=0;j<n;j++) a[i][j]=I[i][j];

  /*free I*/
  for(i=0;i<n;i++) free(I[i]);
  free(I);
}

