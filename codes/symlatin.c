#include <stdio.h>
#include <stdlib.h>

typedef int Boolean;

#include "misc.h"
#include "mt19937ar.h"
#include "symlatin.h"

int generateBooleanCube( Boolean ***X, int **L, Boolean **H, 
			 int n, int r, int seed, int type ){
  int input_ass,i,j,k;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      for(k=0;k<n;k++)
	X[i][j][k] = TRUE;
  input_ass = generatePSLS( L, n, r, seed, type );
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if( L[i][j] != EMPTY )
	for(k=0;k<n;k++){
	  X[i][j][k]       = FALSE;
	  X[i][k][L[i][j]] = FALSE;
	  X[k][j][L[i][j]] = FALSE;
	}
  if( H != NULL )
    reduceByHAT( X, H, n );
  if( type == HOLEHAT || type == HOLERAND )
    makeHole( X, n, r, seed ); 
  return input_ass;
}


int generatePSLS( int **L, int n, int r, int seed, int type ){
  int holes,i;
  switch( type ){
  case QWH:
    generateSLS( L, n, seed );
    holes = removeSLS( L, n, r, seed );
    return n*(n+1)/2-holes;
  case QC:
    return constructSLS( L, n, r, seed );
  case RANDHAT:
  case INPUTHAT:
  case HOLEHAT:
    if( n%2 != 0 ){
      fprintf( stderr, "error: for q=HAT, n should be even.\n" );
      exit( EXIT_FAILURE );
    }
    for(i=0;i<n;i++)
      L[i][i] = n-1;
    return n;
  case HOLERAND:
    return 0;
  default:
    fprintf(stderr,"error: undefined type (%d)\n",type);
    exit(EXIT_FAILURE);
  }
}


void generateSLS( int **L, int n, int seed ){
  int *R,k,l,i,j,tmp;
  polygonMethod( L, n );
  init_genrand( seed );
  // row&col shuffler
  for(k=0;k<n-1;k++){
    l = k+genrand_int31()%(n-k);
    for(i=0;i<n;i++){
      tmp = L[i][k];
      L[i][k] = L[i][l];
      L[i][l] = tmp;
    }
    for(j=0;j<n;j++){
      tmp = L[k][j];
      L[k][j] = L[l][j];
      L[l][j] = tmp;
    }
  }  
  // value shuffler
  R = (int*)malloc(n*sizeof(int));
  for(k=0;k<n;k++)
    R[k] = k;
  shuffle( R, n-1, sizeof(int), genrand_real2 );
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      L[i][j] = R[L[i][j]];
  free( R );
}


void polygonMethod( int **L, int n ){
  int t,k,i,j;
  if( n%2==1 ){
    fprintf( stderr, "error: n should be even.\n" );
    exit( EXIT_FAILURE );
  }
  for(k=0;k<n;k++)
    L[k][k] = n-1;
  for(t=0;t<n-1;t++){
    L[t][n-1] = t;
    L[n-1][t] = t;
    for(k=1;k<n/2;k++){
      i = (t+k) % (n-1);
      j = (t-k+n-1) % (n-1);
      L[i][j] = t;
      L[j][i] = t;
    }
  }
}


int removeSLS( int **L, int n, int r, int seed ){
  int *A,i,j,k=0,maxRemov;
  init_genrand( seed );
  A = (int*)malloc(((n+1)*n/2)*sizeof(int));
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++){
      A[k] = i*n+j;
      k++;
    }
  shuffle( A, (n+1)*n/2, sizeof(int), genrand_real2 );
  maxRemov = ((n+1)*n/2)*(100-r)/100;
  for(k=0;k<maxRemov;k++){
    i = A[k]/n;
    j = A[k]%n;
    if( i == j )
      L[i][i] = EMPTY;
    else{
      L[i][j] = EMPTY;
      L[j][i] = EMPTY;
    }
  }
  return maxRemov;
}

void makeHole( Boolean ***X, int n, int r, int seed ){
  int *A,i,j,k,t,numA=0,maxRemov;
  init_genrand( seed );
  A = (int*)malloc(n*n*n*sizeof(int));
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++)
      for(k=0;k<n;k++)
	if( X[i][j][k] == TRUE ){
	  A[numA] = i*n*n + j*n + k;
	  numA++;
	}
  shuffle( A, numA, sizeof(int), genrand_real2 );
  maxRemov = numA*r/100;
  for(t=0;t<maxRemov;t++){
    i = A[t] / (n*n);
    j = ( A[t] % (n*n) ) / n;
    k = ( A[t] % (n*n) ) % n;
    if( i == j )
      X[i][j][k] = FALSE;
    else{
      X[i][j][k] = FALSE;
      X[j][i][k] = FALSE;
    }
  }
  free( A );
  return;
}



void outputLatinSquare( FILE *out, int **L, int n ){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++)
      if( L[i][j] <= EMPTY )
	fprintf( out, "  -" );
      else
	fprintf( out, " %2d", L[i][j]+1 );
    fprintf( out, "\n" );
  }
}


int isPLS( int **L, int n ){
  int i,j,k,*A;
  A = (int*)malloc(n*sizeof(int));
  for(i=0;i<n;i++){
    for(k=0;k<n;k++)
      A[k] = 0;
    for(j=0;j<n;j++)
      if( L[i][j] != EMPTY ){
	if( A[L[i][j]] )
	  return 0;
	A[L[i][j]] = 1;
      }
  }
  for(j=0;j<n;j++){
    for(k=0;k<n;k++)
      A[k] = 0;
    for(i=0;i<n;i++)
      if( L[i][j] != EMPTY ){
	if( A[L[i][j]] )
	  return FALSE;
	A[L[i][j]] = 1;
      }
  }
  free( A );
  return TRUE;
}


Boolean isSymmetric( int **L, int n ){
  int i,j;
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      if( L[i][j] != L[j][i] )
	return FALSE;
  return TRUE;
}


Boolean isExtension( int **S, int **L, int n ){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if( L[i][j] != EMPTY && L[i][j] != S[i][j] )
	return FALSE;
  return TRUE;
}


Boolean isHATCompatible( int **L, Boolean **H, int n ){
  int i,j;
  if( H == NULL )
    return TRUE;
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++)
      if( i==j && L[i][j] != n-1 )
	return FALSE;
      else if( i!=j && L[i][j] != EMPTY && H[i][L[i][j]] == H[j][L[i][j]] )
	return FALSE;
  return TRUE;
}



int constructSLS( int **L, int n, int r, int seed ){
  Boolean flag;
  int *T,numT=0,i,j,k,l,max,assign,row,col,val;
  init_genrand( seed );
  max = (n*(n+1)/2) * r / 100;
  if( max == 0 )
    return 0;
  T = (int*)malloc((n*n*(n+1)/2)*sizeof(int));
  for(k=0;k<n;k++)
    for(i=0;i<n;i++)
      for(j=0;j<=i;j++){
	T[numT] = numT;
	numT++;
      }

  /*** iteration ***/
  while( 1 ){
    assign = 0;
    // make empty Latin square
    for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	L[i][j] = EMPTY;
    // make shuffled entries
    shuffle( T, numT, sizeof(int), genrand_real2 );
    // assign
    for(k=0;k<numT;k++){
      val = T[k] / (n*(n+1)/2);
      l = -1;
      flag = FALSE;
      for(i=0;i<n;i++){
	for(j=0;j<=i;j++){
	  l++;
	  if( l == T[k] % (n*(n+1)/2) ){
	    row = i;
	    col = j;
	    flag = TRUE;
	    break;
	  }
	}
	if( flag )
	  break;
      }
      if( isAssignable( L, n, row, col, val ) ){
	if( row == col )
	  L[row][col] = val;
	else{
	  L[row][col] = val;
	  L[col][row] = val;
	}
	assign++;
	if( assign==max )
	  break;
      }
    }
    if( assign==max )
      break;
  }
  free( T );
  return max;
}


Boolean isAssignable( int **L, int n, int row, int col, int val ){
  int i,j;
  if( L[row][col] != EMPTY || L[col][row] != EMPTY )
    return FALSE;
  for(i=0;i<n;i++)
    if( i != row && ( L[i][col] == val || L[col][i] == val ) )
      return FALSE;
  for(j=0;j<n;j++)
    if( j != col && ( L[row][j] == val || L[j][row] == val ) )
      return FALSE;
  return TRUE;
}

void makeRandHAT( Boolean **H, int n, int seed ){
  int *P,i,k;
  init_genrand( seed );
  P = (int*)malloc(n*sizeof(int));
  for(i=0;i<n;i++)
    P[i] = i;
  /*** make random HAT ***/
  for(k=0;k<n-1;k++){
    shuffle( P, n, sizeof(int), genrand_real2 );
    /* in slot k,
       P[0],...,P[n/2-1] play at home and
       P[n/2],...,P[n-1] play at away */
    for(i=0;i<n/2;i++)
      H[P[i]][k] = HOME;
    for(i=n/2;i<n;i++)
      H[P[i]][k] = AWAY;
  }
  free( P );

}


void makeHATMBfromArray( int n, int **H, int *A ){
  int i,j;
  for(i=0;i<n/2;i++){
    if( i==0 ){
      for(j=0;j<n-1;j++)
	if( j%2 == 0 )
	  H[i][j] = HOME;
	else
	  H[i][j] = AWAY;
      continue;
    }
    H[i][0] = AWAY;
    for(j=1;j<n-1;j++)
      if( j==A[i-1]+1 )
	H[i][j] = H[i][j-1];
      else{
	if( H[i][j-1] == HOME )
	  H[i][j] = AWAY;
	else
	  H[i][j] = HOME;
      }
	
  }
  for(i=n/2;i<n;i++)
    for(j=0;j<n-1;j++)
      if( H[i-n/2][j] == HOME )
	H[i][j] = AWAY;
      else
	H[i][j] = HOME;

  /**********/
  /*
    int k;
  printf("[");
  for(i=0;i<n/2-1;i++)
    printf(" %d",A[i]);
  printf("]\n\n");
  for(i=0;i<n;i++){
    printf("%d:",i);
    for(k=0;k<n-1;k++)
      printf(" %d",H[i][k]);
    printf("\n");
  }
  printf("---\n");
  */
  /**********/
}


void reduceByHAT( Boolean ***X, Boolean **H, int n ){
  int i,j,k;
  /*** reflect H to X ***/
  for(k=0;k<n-1;k++)
    for(i=0;i<n;i++)
      for(j=i+1;j<n;j++)
	if( H[i][k] == H[j][k] ){
	  X[i][j][k] = FALSE;
	  X[j][i][k] = FALSE;
	}
}

int getNumOfForbiddenMatches( int n, int ***X, int **L ){
  int i,j,k,num=0;
  for(k=0;k<n;k++)
    for(i=0;i<n;i++)
      for(j=0;j<=i;j++)
	if( L[i][j] != k && X[i][j][k] == FALSE )
	  num++;
  return num;
}


int getNumOfFixedMatches( int n, int **L ){
  int i,j,num=0;
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++)
      if( L[i][j] != EMPTY )
	num++;
  return num;
}

