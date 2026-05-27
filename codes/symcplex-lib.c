/********** symcplex-lib.c **********/

#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "symlatin.h"
#include "symcplex-lib.h"

int initCPXSymInstance( CPXInstance I, Boolean ***X, int **L, int n ){
  Nall = (n+1)*n/2;
  I->probname = (char*)mallocE(32*sizeof(char));
  strcpy(I->probname,"PSLSE");
  I->numForbidden = getNumOfForbiddenMatches( n, X, L );
  I->numFixed = getNumOfFixedMatches( n, L );
  I->numNonZero = Nall + 3*n*(n*(n-1)/2) + 2*n*n
    + I->numForbidden + I->numFixed;
  I->numCols = Nall + Nall*n;
  I->numRows = Nall + n*n + I->numForbidden + I->numFixed;
  I->objSense = CPX_MAX;
  I->c = (double*)mallocE(I->numCols*sizeof(double));
  I->b = (double*)mallocE(I->numRows*sizeof(double));
  I->sense = (char*)mallocE(I->numRows*sizeof(char));
  I->matbeg = (int*)mallocE(I->numCols*sizeof(int));
  I->matcnt = (int*)mallocE(I->numCols*sizeof(int));
  I->matind = (int*)mallocE(I->numNonZero*sizeof(int));
  I->matval = (double*)mallocE(I->numNonZero*sizeof(double));
  I->LB = (double*)mallocE(I->numCols*sizeof(double));
  I->UB = (double*)mallocE(I->numCols*sizeof(double));
  I->Range = NULL;
  I->ctype = (char*)mallocE(I->numCols*sizeof(char));
  return 0;
}

int setCPXSymInstance( CPXInstance I, Boolean ***X, int **L, int n ){
  int p,q,i,j,k,r,cntForbidden=0,cntFixed=0;

  /*** coefficients of objective ***/
  for(q=0;q<I->numCols;q++)
    if( q < Nall )
      I->c[q] = 1.;
    else
      I->c[q] = 0.;

  /*** righthand side values of constraints ***/
  for(p=0;p<I->numRows;p++)
    if( p < Nall || ( p >= Nall+n*n &&
		      p <  Nall+n*n+I->numForbidden ) )
      I->b[p] = 0.;
    else
      I->b[p] = 1.;

  /*** sign of each constraint: 'L' as <=, 'E' as ==, 'G' as >=
       note: >= was indicated by 'R' in v12.4 and earlier versions ***/
  for(p=0;p<I->numRows;p++)
    if( p < Nall )
      I->sense[p] = 'E';
    else if( p < Nall+n*n )
      I->sense[p] = 'L';
    else
      I->sense[p] = 'E';

  /*** the starting index and the number of non-zero coeffs in each column ***/
  i = 0;
  j = 0;
  k = 0;
  for(q=0;q<I->numCols;q++){
    /* the starting index */
    if( q == 0 )
      I->matbeg[0] = 0;
    else
      I->matbeg[q] = I->matbeg[q-1] + I->matcnt[q-1];
    /* the number of non-zero coeffs */
    if( q < Nall )
      I->matcnt[q] = 1;
    else{
      if( i == j )
	I->matcnt[q] = 2;
      else
	I->matcnt[q] = 3;
      if( L[i][j] == k )
	(I->matcnt[q])++;
      else if( X[i][j][k] == FALSE )
	(I->matcnt[q])++;
      // go on to the next triple
      getNextCell( n, &i, &j, &k );
    }
  }

  /*** row indices and non-zero coeffs in each column ***/
  r = 0;
  i = 0;
  j = 0;
  k = 0;
  for(q=0;q<I->numCols;q++){
    /* y_{i,j} */
    if( q < Nall ){
      I->matind[r] = q;
      I->matval[r] = -1.0;
      r++;
      continue;
    }
    /* x_{i,j,k} */
    // correspondance with y_{i,j}
    I->matind[r] = (i+1)*i/2+j;
    I->matval[r] = 1.0;
    r++;
    // for each (min(i,j),k)...
    I->matind[r] = Nall + k*n + smaller(i,j);
    I->matval[r] = 1.0;
    r++;
    // for each (max(i,j),k)...
    if( i != j ){
      I->matind[r] = Nall + k*n + larger(i,j);
      I->matval[r] = 1.0;
      r++;
    }
    // fixed matches
    if( L[i][j] == k ){
      I->matind[r] = Nall + n*n + I->numForbidden + cntFixed;
      I->matval[r] = 1.0;
      r++;
      cntFixed++;
    }
    // forbidden matches
    else if( X[i][j][k] == FALSE ){
      I->matind[r] = Nall + n*n + cntForbidden;
      I->matval[r] = 1.0;
      r++;
      cntForbidden++;
    }
    // go on to the next triple
    getNextCell( n, &i, &j, &k );
  }

  /*** lower and upper bounds and types of the variables ***/
  for(q=0;q<I->numCols;q++){
    I->LB[q] = 0.0;
    I->UB[q] = 1.0;
  }    
  return 0;
}


void getSolutionFromCPXInstance( CPXENVptr env, CPXLPptr IP, CPXInstance I, int **S, int n ){
  double *x;
  int i=0,j=0,k=0,q;
  x = (double*)malloc(I->numCols*sizeof(double));
  if( CPXgetx( env, IP, x, 0, I->numCols-1 ) ){
    fprintf( stderr, "error: unfailed to call CPXgetx()\n" );
    exit( EXIT_FAILURE );
  }
  for(q=Nall;q<I->numCols;q++){
    if( x[q] > 0.99 ){
      S[i][j] = k;
      S[j][i] = k;
    }
    getNextCell( n, &i, &j, &k );
  }
  free( x );
}


void getNextCell( int n, int *ip, int *jp, int *kp ){
  (*jp)++;
  if( *jp>*ip ){
    *jp = 0;
    (*ip)++;
    if( *ip==n ){
      *jp = 0;
      *ip = 0;
      (*kp)++;
    }
  }
}


void freeCPXSymInstance( CPXInstance I ){
  free( I->probname );
  free( I->c );
  free( I->b );
  free( I->sense );
  free( I->matbeg );
  free( I->matcnt );
  free( I->matind );
  free( I->matval );
  free( I->LB );
  free( I->UB );
}

