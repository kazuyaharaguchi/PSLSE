/******************************
   hatfeas_ls.c

   search for all feasible HATs with minimum breaks
   & execute local search on them 
******************************/

#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "symlatin.h"
#include "elementary.h"
#include "greedy.h"
#include "localsearch.h"
#include "misc.h"
#include "symcplex-lib.h"

typedef struct _Hash Hash;
struct _Hash{
  int *A;
  int counter;
  Boolean isFeasible;
  struct _Hash *next;
};



void processArgsCPLEX( int argc, char *argv[], int *np, double *timelimitp );
void initILSParamForHatfeas( ILSParam *IP, int n, double timelimit );
Hash *getHash( Hash *T, int Tmax, int *A, int n, int m );
Boolean isEquivHash( int *A, int *B, int n, int m );
int getHashValue( int *A, int n, int m );


int main( int argc, char *argv[] ){
  CPXENVptr     env = NULL;
  CPXLPptr      LP = NULL;
  CPXInstance   I;
  Solution      Sol;
  Hash          *T,*h;
  ILSParam      IParam;
  ILSRec        IRec;
  Boolean       ***X,**H=NULL;
  double        objval,timelimit;
  int           Tmax;
  int           status,i,j,n;
  int           **L,*A,input_ass,numLPsatisfy=0,numSkipped=0,p,t=0;

  /*** process arguments ***/
  processArgsCPLEX( argc, argv, &n, &timelimit );

  /*** memory allocation ***/
  X = (Boolean***)malloc(n*sizeof(Boolean**));
  for(i=0;i<n;i++){
    X[i] = (Boolean**)malloc(n*sizeof(Boolean*));
    for(j=0;j<n;j++){
      X[i][j] = (Boolean*)malloc(n*sizeof(Boolean));
    }
  }
  L = (int**)malloc(n*sizeof(int*));
  for(i=0;i<n;i++){
    L[i] = (int*)malloc(n*sizeof(int));
    for(j=0;j<n;j++)
      L[i][j] = EMPTY;
  }
  H = (Boolean**)malloc(n*sizeof(Boolean*));
  for(i=0;i<n;i++)
    H[i] = (Boolean*)malloc((n-1)*sizeof(Boolean));
  I = (CPXInstance)malloc(sizeof(struct _CPXInstance));
  A = (int*)malloc((n/2-1)*sizeof(int));
  initCombi(n-2,n/2-1,A);
  Nall = (n+1)*n/2;
  initILSParamForHatfeas( &IParam, n, timelimit );
  Tmax = (n-2)*(n-2)*(n-2) + (n-2)*(n-2) + (n-2);
  T = (Hash*)malloc(Tmax*sizeof(Hash));
  for(p=0;p<Tmax;p++){
    T[p].A = NULL;
    T[p].counter = 0;
    T[p].isFeasible = UNDEF;
    T[p].next = NULL;
  }
  
  /*** open CPLEX environment ***/
  env = (CPXENVptr)CPXopenCPLEX( &status );
  if( env == NULL ){
    fprintf( stderr, "error: unable to open the CPLEX environment (errno=%d)\n", status );
    exit( 1 );
  }
  status = CPXsetintparam( env, CPXPARAM_LPMethod, CPX_ALG_AUTOMATIC );
  status = CPXsetintparam( env, CPXPARAM_Read_DataCheck, CPX_ON );
  status = CPXsetdblparam( env, CPX_PARAM_CLOCKTYPE, 1 );
  status = CPXsetdblparam( env, CPX_PARAM_TILIM, timelimit );

  /*** check all combinations ***/
  while( 1 ){
    /*** increase the counter ***/
    t++;
    printf("[%d]",t);

    /*** get hash ***/
    h = getHash( T, Tmax, A, n, n/2-1 );
    if( h->isFeasible != UNDEF ){
      numSkipped++;
      if( !getNextCombi( n-2, n/2-1, A ) )
	break;
      printf("\n");
      continue;
    }

    printf("*\n");

    /*** make HAT from A ***/
    makeHATMBfromArray( n, H, A );

    /***** generate instance *****/
    input_ass = generateBooleanCube( X, L, H, n, 0, 1, INPUTHAT );

    /*** create CPLEX instance ***/
    initCPXSymInstance( I, X, L, n );
    setCPXSymInstance( I, X, L, n );
    LP = (CPXLPptr)CPXcreateprob( env, &status, I->probname );
    status = CPXcopylp( env, LP, I->numCols, I->numRows, I->objSense, I->c, I->b,
			I->sense, I->matbeg, I->matcnt, I->matind, I->matval,
			I->LB, I->UB, NULL );

    //CPXwriteprob( env, LP, "TETE.lp", NULL );

    /*** solve LP ***/
    status = CPXlpopt( env, LP );

    /***** output *****/
    status = CPXgetobjval( env, LP, &objval );
    if( objval > (double)(n+1.0)*n/2.0 - 0.5 ){
      h->isFeasible = TRUE;
      numLPsatisfy++;
      initSolution( n, X, L, input_ass, &Sol );
      greedy( &Sol, IParam.init, IParam.sinit );
      initILSRec( &IRec, &Sol ); 
      printf("ILS starts... %d --->\n",Sol.num[NUM_SOL]+Sol.input_ass);
      fflush( stdout );
      ILS( &IParam, &IRec, &Sol );
      printf("%d\t%d\n",Sol.num[NUM_SOL]+Sol.input_ass,(n+1)*n/2);
      if(Sol.num[NUM_SOL]+Sol.input_ass<(n+1)*n/2)
	printf("give up\n");
      freeILSRec( n, &IRec );
      freeSolution( n, &Sol );
    }
    else
      h->isFeasible = FALSE;

    if( !getNextCombi( n-2, n/2-1, A ) )
      break;
    CPXfreeprob( env, &LP );
    free( LP );
    freeCPXSymInstance( I );
   }

  printf("%d\t%d\n",numLPsatisfy,t);
  
  numLPsatisfy=0;
  for(p=0;p<Tmax;p++){
    if( T[p].A == NULL )
      continue;
    h = &(T[p]);
    while( h != NULL ){
      if( h->isFeasible == TRUE )
	numLPsatisfy += h->counter;
      h = h->next;
    }
  }
  printf("numLPsatisfy: %d\n",numLPsatisfy);
  printf("numSkipped: %d\n",numSkipped);

  freeILSParam( &IParam );
  
  return 0;
}


void processArgsCPLEX( int argc, char *argv[], int *np, double *timelimitp ){
  if( argc < 3 ){
    fprintf( stderr, "usage: %s (n)(timelimit)\n\n",
	     argv[0] );
    exit( 1 );
  }
  *np = atoi( argv[1] );
  *timelimitp = atof( argv[2] );
}


void initILSParamForHatfeas( ILSParam *IP, int n, double timelimit ){
  IP->n = n;
  IP->q = INPUTHAT;
  IP->r = 0;
  IP->alg = ILS2;
  IP->init = G1;
  IP->itrmax = DEFAULT_ITRMAX;
  IP->output = (char*)malloc(LENGTH_OF_FILENAME*sizeof(char));
  strcpy( IP->output, "" );
  IP->s = 1;
  IP->seed = 1;
  IP->sinit = 1;
  IP->timemax = timelimit;
}

Hash *getHash( Hash *T, int Tmax, int *A, int n, int m ){
  Hash *current,*new;
  int h,i;
  h = getHashValue( A, n, m );
  if( T[h].A == NULL ){
    T[h].A = (int*)malloc(m*sizeof(int));
    for(i=0;i<m;i++)
      T[h].A[i] = A[i];
    T[h].counter = 1;
    T[h].isFeasible = UNDEF;
    T[h].next = NULL;
    return &(T[h]);
  }
  current = &(T[h]);
  while( 1 ){
    if( isEquivHash( A, current->A, n, m ) ){
      (current->counter)++;
      break;
    }
    if( current->next == NULL ){
      new = (Hash*)malloc(sizeof(Hash));
      new->A = (int*)malloc(m*sizeof(int));
      for(i=0;i<m;i++)
	new->A[i] = A[i];
      new->counter = 1;
      new->isFeasible = UNDEF;
      new->next = NULL;
      current->next = new;
      current = new;
      break;
    }
    current = current->next;    
  }
  return current;
}


int getHashValue( int *A, int n, int m ){
  int *B,i,val;
  B = (int*)malloc((m+1)*sizeof(int));
  B[0] = A[0]+1;
  for(i=1;i<m;i++)
    B[i] = A[i]-A[i-1];
  B[m] = n-2-A[m-1];
  qsort( B, m+1, sizeof(int), cmpNegInt );
  val = B[0]*m*m + B[1]*m;
  if( m>1 )
    val += B[2];
  free( B );
  return val;
}


Boolean isEquivHash( int *A, int *B, int n, int m ){
  Boolean Flag=FALSE,flag;
  int *X,*Y,i,r,base;
  X = (int*)malloc((m+1)*sizeof(int));
  Y = (int*)malloc((m+1)*sizeof(int));
  X[0] = 0;
  Y[0] = 0;
  for(i=1;i<=m;i++){
    X[i] = A[i-1]+1;
    Y[i] = B[i-1]+1;
  }
  for(r=0;r<=m;r++){
    // rotate Y
    if( r>0 ){
      base = Y[r];
      for(i=0;i<=m;i++){
	Y[i] -= base;
	if( Y[i]<0 )
	  Y[i] += n-1;
      }
    }
    // check whether X and Y are the same
    flag = TRUE;
    for(i=0;i<=m;i++)
      if( X[i] != Y[(r+i)%(m+1)] ){
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
  return Flag;
}
