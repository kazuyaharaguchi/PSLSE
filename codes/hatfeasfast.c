/******************************
   hatfeasfast.c

   search for all feasible HATs with minimum breaks & output them 
   feasibility check is done by MIM algorithm in PATAT 2002
******************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "misc.h"
#include "symlatin.h"

#define TMAX 1000000

typedef int Boolean;
typedef struct _Hash Hash;
struct _Hash{
  int *A;
  struct _Hash *next;
};

void searchDist( int *D, int m, int k, int lev, int residue );
void searchVector( int *D, int m, int k, int *V, int d, int begin );
Boolean satisfyHATMBCondition( int n, Boolean **H );
int getAlphaValue( int n, int **S, int i, int k );
Boolean existHash( int *A, int m, int k );
Boolean isEquivHash( int *A, int *B, int m, int k );
int getHashValue( int *A, int m, int k );

Hash T[TMAX];

int main( int argc, char *argv[] ){
  int     *D,d,k,m,n,t;

  /***** process arguments *****/
  if( argc < 2 ){
    fprintf( stderr, "usage: %s (n)\n\n", argv[0] );
    exit( 1 );
  }
  n = atoi( argv[1] );
  m = n-1;
  k = n/2;
  for(t=0;t<TMAX;t++){
    T[t].A = NULL;
    T[t].next = NULL;
  }
  D = (int*)malloc((m-k+1)*sizeof(int));
  for(d=0;d<=m-k;d++)
    D[d] = 0;
  searchDist( D, m, k, m-k, m-k );
  return 0;
}




void searchDist( int *D, int m, int k, int lev, int residue ){
  int *V,d,q,v;
  if( lev == 1 ){
    D[lev] = residue;
    V = (int*)malloc(k*sizeof(int));
    for(v=0;v<k;v++)
      V[v] = 0;
    for(d=m-k;d>0;d--)
      if( D[d] > 0 ){
	V[k-1] = d;
	break;
      }
    D[d]--;
    searchVector( D, m, k, V, d, 0 );
    D[d]++;
    free( V );
    return;
  }
  for(q=residue/lev;q>=0;q--){
    D[lev] = q;
    searchDist( D, m, k, lev-1, residue-q*lev );
  }
}


void searchVector( int *D, int m, int k, int *V, int d, int begin ){
  int l;
  if( d == 0 ){
    Boolean **H=NULL;
    int     *A,i;
    // home-away pattern
    H = (Boolean**)malloc((m+1)*sizeof(Boolean*));
    for(i=0;i<m+1;i++)
      H[i] = (Boolean*)malloc(m*sizeof(Boolean));
    // combination
    A = (int*)malloc((k-1)*sizeof(int));
    A[0] = V[0];
    for(l=1;l<k-1;l++)
      A[l] = A[l-1] + V[l] + 1;
    // make HAT from A
    makeHATMBfromArray( m+1, H, A );
    // examine feasibility
    if( satisfyHATMBCondition( m+1, H ) ){
      if( existHash( A, m, k ) == FALSE )
	for(i=0;i<k-1;i++){
	  printf("%d",A[i]);
	  if(i<k-2)
	    printf(",");
	  else{
	    printf("\n");
	    fflush( stdout );
	  }
	}
    }
    free( A );
    for(i=0;i<m+1;i++)
      free( H[i] );
    free( H );
    return;
  }
  if( D[d] == 0 )
    return searchVector( D, m, k, V, d-1, 0 );
  for(l=begin;l<k-1;l++){
    if( V[l] != 0 )
      continue;
    V[l] = d;
    D[d]--;
    searchVector( D, m, k, V, d, l+1 );
    D[d]++;
    V[l] = 0;
  }
}


Boolean satisfyHATMBCondition( int n, Boolean **H ){
  Boolean flag=TRUE;
  int **S,i,k,j;
  
  /***** construct cumulative sum table *****/
  S = (int**)malloc(n*sizeof(int*));
  for(i=0;i<n;i++)
    S[i] = (int*)malloc((n-1)*sizeof(int));
  for(j=0;j<n-1;j++){
    if( H[0][j] == HOME )
      S[0][j] = 1;
    else
      S[0][j] = 0;
    for(i=1;i<n;i++)
      if( H[i][j] == HOME )
	S[i][j] = S[i-1][j]+1;
      else
	S[i][j] = S[i-1][j];
  }

  /***** evaluate whether the condition is satisfied *****/
  for(i=0;i<n/2;i++){
    for(k=i+1;k<i+n/2;k++)
      if( getAlphaValue( n, S, i, k ) < 0 ){
	flag = FALSE;
	break;
      }
    if( !flag )
      break;
  }
  
  /***** postprocess *****/
  for(i=0;i<n;i++)
    free( S[i] );
  free( S );
  
  return flag;
}


int getAlphaValue( int n, int **S, int i, int k ){
  int alpha,j,m,home,away;
  m = k-i+1;
  alpha = -(m*(m-1))/2;
  for(j=0;j<n-1;j++){
    if( i == 0 )
      home = S[k][j];
    else
      home = S[k][j]-S[i-1][j];
    away = m-home;
    if( home < away )
      alpha += home;
    else
      alpha += away;
  }
  return alpha;
}


Boolean existHash( int *A, int m, int k ){
  Hash *current,*new;
  int h,i;
  h = getHashValue( A, m, k );
  if( T[h].A == NULL ){
    T[h].A = (int*)malloc((k-1)*sizeof(int));
    for(i=0;i<k-1;i++)
      T[h].A[i] = A[i];
    T[h].next = NULL;
    return FALSE;
  }
  current = &(T[h]);
  while( 1 ){
    if( isEquivHash( A, current->A, m, k ) )
      break;
    if( current->next == NULL ){
      new = (Hash*)malloc(sizeof(Hash));
      new->A = (int*)malloc((k-1)*sizeof(int));
      for(i=0;i<k-1;i++)
	new->A[i] = A[i];
      new->next = NULL;
      current->next = new;
      return FALSE;
    }
    current = current->next;    
  }
  return TRUE;
}


int getHashValue( int *A, int m, int k ){
  int *B,i,val;
  B = (int*)malloc(k*sizeof(int));
  B[0] = A[0]+1;
  for(i=1;i<k-1;i++)
    B[i] = A[i]-A[i-1];
  B[k-1] = m-A[k-2];
  qsort( B, k, sizeof(int), cmpNegInt );
  val = B[0]*m*m*m + B[1]*m*m + B[2];
  free( B );
  return val%TMAX;
}


Boolean isEquivHash( int *A, int *B, int m, int k ){
  Boolean Flag=FALSE,flag;
  int *X,*Y,i,r,base;
  X = (int*)malloc(k*sizeof(int));
  Y = (int*)malloc(k*sizeof(int));
  X[0] = 0;
  Y[0] = 0;
  for(i=1;i<k;i++){
    X[i] = A[i-1]+1;
    Y[i] = B[i-1]+1;
  }

  /*
  printf("---\n");
  for(i=0;i<k-1;i++)
    printf(" %d",A[i]);
  printf("\t--->\t");
  for(i=0;i<k;i++)
    printf(" %d",X[i]);
  printf("\n");
  for(i=0;i<k-1;i++)
    printf(" %d",B[i]);
  printf("\t--->\t");
  for(i=0;i<k;i++)
    printf(" %d",Y[i]);
  printf("\n");
  */

  for(r=0;r<k;r++){
    // rotate Y
    if( r>0 ){
      base = Y[r];
      for(i=0;i<k;i++){
	Y[i] -= base;
	if( Y[i]<0 )
	  Y[i] += m;
      }
    }
    // check whether X and Y are the same
    flag = TRUE;
    for(i=0;i<k;i++)
      if( X[i] != Y[(r+i)%k] ){
	flag = FALSE;
	break;
      }
    if( flag ){
      Flag = TRUE;
      break;
    }
  }
  free( X );
  free( Y );

  //printf("%d\n\n",Flag);

  return Flag;
}

