/******************************
   symCPLEX.c
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

void outputUsage( FILE *out, char *cmd );
void outputErrorMsg( char *msg );
void processArgsCPLEX( int argc, char *argv[], ILSParam *IP );
void MIPstart( CPXENVptr env, CPXLPptr IP, CPXInstance I, int **L,
	       Solution *S );

int main( int argc, char *argv[] ){
  CPXENVptr   env = NULL;
  CPXLPptr    IP = NULL;
  CPXInstance I;
  Solution    Sol;
  ILSParam    IParam;
  Boolean     ***X,**H=NULL;
  double      objval;
  double      timeGrdStart,timeGrdEnd;
  double      start_cplextime,finish_cplextime;
  int         status,solstat,i,j,k,n;
  int         **L,**S,input_ass,init_ass;

  /***** process arguments *****/
  initILSParam( &IParam );
  processArgsCPLEX( argc, argv, &IParam );

  /*** memory allocation ***/
  n = IParam.n;
  X = (Boolean***)malloc(n*sizeof(Boolean**));
  for(i=0;i<n;i++){
    X[i] = (Boolean**)malloc(n*sizeof(Boolean*));
    for(j=0;j<n;j++){
      X[i][j] = (Boolean*)malloc(n*sizeof(Boolean));
      for(k=0;k<n;k++)
	X[i][j][k] = TRUE;
    }
  }
  L = (int**)malloc(n*sizeof(int*));
  for(i=0;i<n;i++){
    L[i] = (int*)malloc(n*sizeof(int));
    for(j=0;j<n;j++)
      L[i][j] = EMPTY;
  }
  if( IParam.q == RANDHAT || IParam.q == INPUTHAT || IParam.q == HOLEHAT ){
    H = (Boolean**)malloc(n*sizeof(Boolean*));
    for(i=0;i<n;i++)
      H[i] = (Boolean*)malloc((n-1)*sizeof(Boolean));   
    if( IParam.q == RANDHAT )
      makeRandHAT( H, n, IParam.s );
    else{
      char **str;
      int *A;
      if( IParam.hatarray == NULL ){
	fprintf( stderr, "error: hatarray should be specified.\n" );
	exit( EXIT_FAILURE );
      }
      str = split( IParam.hatarray, ',' );
      A = (int*)malloc((n/2-1)*sizeof(int));
      for(i=0;i<n/2-1;i++){
	A[i] = atoi( str[i] );
	free( str[i] );
      }
      free( str );
      makeHATMBfromArray( n, H, A );
      free( A );
    }    
  }


  /***** initialize Sol *****/
  input_ass = generateBooleanCube( X, L, H, n, IParam.r, IParam.s, IParam.q );
  printf("Input:\t%d\n",input_ass);
  initSolution( n, X, L, input_ass, &Sol );

  
  /***** construct initial solution *****/
  timeGrdStart = cpu_time();
  if( IParam.init != UNDEF )
    init_ass = greedy( &Sol, IParam.init, IParam.sinit ); 
  else
    init_ass = 0;
  timeGrdEnd = cpu_time();
  printf("Input+Init:\t%d\n",input_ass+init_ass);
  printf("InitCpuTime:\t%g\n",timeGrdEnd-timeGrdStart);
  printf("UpperBound:\t%d\n",input_ass+Sol.UB);


#ifdef DEBUG
  outputLatinSquare( stdout, L, n );
  if( !isPLS( L, n ) ){
    fprintf( stderr, "error: instance is not a PLS\n" );
    return 1;
  }
  if( !isSymmetric( L, n ) ){
    fprintf( stderr, "error: instance is not symmetric\n" );
    return 1;
  }
#endif


  /***** solve PLSE by integer programming model *****/
  /*** open CPLEX environment ***/
  env = (CPXENVptr)CPXopenCPLEX( &status );
  if( env == NULL ){
    fprintf( stderr, "error: unable to open the CPLEX environment (errno=%d)\n", status );
    exit( 1 );
  }
  status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
  if( env == NULL ){
    fprintf( stderr, "error: unable to set CPX_PARAM_DATA_CHECK (errno=%d)\n", status );
    exit( 1 );
  }

  /*** create CPLEX instance from L ***/
  Nall = (n+1)*n/2;
  I = (CPXInstance)malloc(sizeof(struct _CPXInstance));
  initCPXSymInstance( I, X, L, n );
  setCPXSymInstance( I, X, L, n );
  for(j=0;j<I->numCols;j++)
    I->ctype[j] = 'B'; // B:binary, C:continuous
  IP = (CPXLPptr)CPXcreateprob( env, &status, I->probname );
  status = CPXcopylp( env, IP, I->numCols, I->numRows, I->objSense, I->c, I->b,
		      I->sense, I->matbeg, I->matcnt, I->matind, I->matval,
		      I->LB, I->UB, NULL );
  if( status ){
    fprintf( stderr, "error: failed to copy problem data (errno=%d)\n", status );
    exit( EXIT_FAILURE );
  }
  status = CPXcopyctype( env, IP, I->ctype );
  if( status ){
    fprintf( stderr, "error: failed to copy ctype (errno=%d)\n", status );
    exit( EXIT_FAILURE );
  }

  /*** compute optimal solution ***/
  // MIP start (give initial solution)
  if( IParam.init != UNDEF )
    MIPstart( env, IP, I, L, &Sol );

  // set the timelimit (CPX_PARAM_CLOCKTYPE=1 indicates CPU time, not elasped time)
  if( IParam.timemax >= 0.0 ){
    status = CPXsetdblparam( env, CPX_PARAM_CLOCKTYPE, 1 );
    status = CPXsetdblparam( env, CPX_PARAM_TILIM, IParam.timemax );
  }
  
#ifdef DEBUG
  CPXwriteprob( env, IP, "TETE.lp", NULL );
#endif

  // 計測開始
  CPXgettime( env, &start_cplextime );
  // 最適化
  status = CPXmipopt( env, IP );
  // 計算終了
  CPXgettime( env, &finish_cplextime );
  if( status ){
    fprintf (stderr, "Failed to optimize MIP.\n");
    exit(1);
  }

  /***** 計算結果の出力 *****/
  solstat = CPXgetstat( env, IP );
  status = CPXgetobjval( env, IP, &objval );
  if( status )
    printf( "CPXOutput:\tNO_SOLUTION_IS_FOUND\t-\n" );
  else{
    printf( "CPXOutput:\t%g\t", objval );
    if( solstat == CPXMIP_OPTIMAL )
      printf( "+\n" );
    else
      printf( "-\n" );
  }
  printf( "CPXCpuTime:\t%g\n",finish_cplextime - start_cplextime );

  S = (int**)mallocE(n*sizeof(int*));
  for(i=0;i<n;i++){
    S[i] = (int*)mallocE(n*sizeof(int));
    for(j=0;j<n;j++)
      S[i][j] = L[i][j];
  }
  getSolutionFromCPXInstance( env, IP, I, S, n ); 
  if( !isPLS( S, n ) ){
    fprintf( stderr, "error: found solution is not a PLS\n" );
    return 1;
  }
  if( !isSymmetric( S, n ) ){
    fprintf( stderr, "error: found solution is not symmetric\n" );
    return 1;
  }
  if( !isExtension( S, L, n ) ){
    fprintf( stderr, "error: found solution is not an extension of given PLS\n" );
    return 1;
  }
  if( !isHATCompatible( S, H, n ) ){
    fprintf( stderr, "error: found solution is not compatible with HAT\n" );
    return 1;
  }
#ifdef DEBUG
  outputLatinSquare( stdout, S, n );
#endif

  return 0;
}


void outputUsage( FILE *out, char *cmd ){
  fprintf( out, "\nCPLEX IP solver for partial symmetric Latin square extension problem\n\n" );
  fprintf( out, "Usage: %s -n=<INT> [OPTIONS]\n\n", cmd );
  fprintf( out, "  -n=<INT>        Order (i.e., grid length)\n\n" );
  fprintf( out, "Options:\n" );
  fprintf( out, "  -q=<INT>        Problem type: {0,...,6}--->{QWH,QC,RANDHAT,INPUTHAT,HOLERAND,HOLEHAT,HOLEFILE} (default:%d)\n", DEFAULT_Q );
  fprintf( out, "  -r=<INT>        Pre-assigned ratio in percent (default:%d)\n", DEFAULT_R );
  fprintf( out, "  -s=<INT>        Random seed for the randinst generator (default:1)\n\n" );
  fprintf( out, "  -hatarray=<STR>\n" );
  fprintf( out, "  -init=<INT>     Initial solution generator: {-1,0,1,2}--->{-,Rand,G1,G5} (default:%d)\n", DEFAULT_INIT );
  fprintf( out, "  -output=<FILE>  Write the output solution to <FILE> (default:none)\n" );
  fprintf( out, "  -sinit=<INT>    Random seed for initial solution generator (default:1)\n" );
  fprintf( out, "  -timemax=<INT>  Time limit in sec (default:%d)\n", DEFAULT_TIMEMAX );
}


void outputErrorMsg( char *msg ){
  fprintf( stderr, "\nError: %s\n\n", msg );
  exit( EXIT_FAILURE );
}

void processArgsCPLEX( int argc, char *argv[], ILSParam *IP ){
  char **S;
  int k,eqs;
  if( argc < 2 ){
    outputUsage( stdout, argv[0] );
    exit( EXIT_FAILURE );
  }

  /*** read the arguments ***/
  for(k=1;k<argc;k++){
    eqs = getNumChar( argv[k], '=' );
    if( argv[k][0] != '-' || eqs > 1 )
      outputErrorMsg("invalid argument");
    S = split( argv[k]+1, '=' );
    if( k == 1 ){
      if( strcmp( S[0], "n" ) != 0 || eqs != 1 )
	outputErrorMsg("grid length should be specified appropriately.");
      IP->n = atoi( S[1] );
    }
    else{
      if( eqs != 1 || strlen( S[1] ) == 0 )
	outputErrorMsg("option value should be specified appropriately.");
      else if( strcmp( S[0], "q" ) == 0 )
	IP->q = atoi( S[1] );
      else if( strcmp( S[0], "r" ) == 0 )
	IP->r = atoi( S[1] );
      else if( strcmp( S[0], "hatarray" ) == 0 ){
	if( IP->hatarray != NULL )
	  free( IP->hatarray );
	IP->hatarray = (char*)malloc(LENGTH_OF_STRING*sizeof(char));
	strcpy( IP->hatarray, S[1] );
      }
      else if( strcmp( S[0], "init" ) == 0 )
	IP->init = atoi( S[1] );
      else if( strcmp( S[0], "itrmax" ) == 0 )
	IP->itrmax = atoi( S[1] );
      else if( strcmp( S[0], "output" ) == 0 ){
	if( IP->output != NULL )
	  free( IP->output );
	IP->output = (char*)malloc(LENGTH_OF_FILENAME*sizeof(char));
	strcpy( IP->output, S[1] );
      }
      else if( strcmp( S[0], "s" ) == 0 )
	IP->s = atoi( S[1] );
      else if( strcmp( S[0], "sinit" ) == 0 )
	IP->sinit = atoi( S[1] );
      else if( strcmp( S[0], "timemax" ) == 0 )
	IP->timemax = atoi( S[1] );
      else
	outputErrorMsg("illegal parameter name.");
    }
    free( S[0] );
    if( eqs == 1 )
      free( S[1] );
    free( S );
  }
  /*** automatically determined parameters ***/
  IP->N = IP->n + IP->ext; 
  return;
}


void MIPstart( CPXENVptr env, CPXLPptr IP, CPXInstance I, 
	       int **L, Solution *S ){
  int *beg,*varindices,*effortlevel,status;
  double *values;
  int p,q,i,j,k;
  beg = (int*)malloc(sizeof(int));
  varindices = (int*)malloc(I->numCols*sizeof(int));
  values = (double*)malloc(I->numCols*sizeof(double));
  effortlevel = (int*)malloc(sizeof(int));
  beg[0] = 0;
  effortlevel[0] = CPX_MIPSTART_CHECKFEAS;
  for(q=0;q<I->numCols;q++){
    varindices[q] = q;
    values[q] = 0.0;
  }
  for(i=0;i<S->n;i++)
    for(j=0;j<=i;j++){
      if( L[i][j] == EMPTY )
	continue;
      k = L[i][j];
      /*** y_{i,j} ***/
      q = ((i+1)*i)/2 + j;
      values[q] = 1.0;
      /*** x_{i,j,k} ***/
      q = Nall + Nall*k + ((i+1)*i)/2 + j;
      values[q] = 1.0;
    }
  
  for(p=0;p<S->num[NUM_SOL];p++){
    i = S->Pi[p]->pos[ROW];
    j = S->Pi[p]->pos[COL];
    k = S->Pi[p]->pos[VAL];
    /*** y_{i,j} ***/
    q = ((i+1)*i)/2 + j;
    values[q] = 1.0;
    /*** x_{i,j,k} ***/
    q = Nall + Nall*k + ((i+1)*i)/2 + j;
    values[q] = 1.0;
  }

  status = CPXsetintparam( env, CPX_PARAM_REPAIRTRIES, -1 );
  if( CPXaddmipstarts( env, IP, 1, I->numCols, beg, varindices,
		       values, effortlevel, NULL ) ){
    fprintf( stderr, "error: failed to execute CPXaddmipstarts\n" );
    exit( EXIT_FAILURE );
  }
  free( beg );
  free( varindices );
  free( values );
  free( effortlevel );

#ifdef DEBUG
  int rows=0,cols=0;
  CPXrefinemipstartconflict( env, IP, 0, &rows, &cols );
  printf("CPXrefinemipstartconflict ... <%d,%d>\n",rows,cols);

  /*
  if( rows || cols ){
    int confstat,confnumrows,confnumcols;
    int *rowind,*rowbdstat,*colind,*colbdstat;
    rowind = (int*)malloc(1000000*sizeof(int));
    rowbdstat = (int*)malloc(1000000*sizeof(int));
    colind = (int*)malloc(1000000*sizeof(int));
    colbdstat = (int*)malloc(1000000*sizeof(int));

    CPXgetconflict( env, IP, &confstat, rowind, rowbdstat,
		    &confnumrows, colind, colbdstat, &confnumcols );
    CPXclpwrite( env, IP, "TETE.lp" );
    printf("%d\n",confstat);
    printf("%d\t%d\n",rows,confnumrows);
    printf("%d\t%d\n",cols,confnumcols);
    exit(1);
  }
  */
#endif


}
