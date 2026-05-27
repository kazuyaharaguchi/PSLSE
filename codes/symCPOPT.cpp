/********************
   symCPOPT.cpp
********************/

#include <ilcp/cp.h>
#include <ilcp/cpext.h> 
#include <stdlib.h>

typedef int Boolean;

#include "define.h"
#include "symlatin.h"
#include "elementary.h"
#include "greedy.h"
#include "localsearch.h"
#include "misc.h"

#define TIMELIMIT 10.0
#define WORKERS 1
#define VERBOSE IloCP::Quiet
#define DEBUG
#undef DEBUG


typedef IloArray<IloIntVarArray> Assignment;
typedef struct{
  int    workers;
  IloInt verbose;
} CPOPTConfig;


extern double cpu_time( void );
void outputUsage( FILE *out, char *cmd );
void outputErrorMsg( const char *msg );
void processArgsCPOPT( int argc, char *argv[], ILSParam *IP, CPOPTConfig *CP );


int main( int argc, char *argv[] ){
  CPOPTConfig cfg;
  Solution   Sol;
  ILSParam   IParam;
  double     start,finish;
  double     timeGrdStart,timeGrdEnd;
  Boolean    ***X,**H=NULL;
  int        **L,**S,**P,input_ass,init_ass;
  int        i,j,k,n;

  /********** process arguments **********/
  initILSParam( &IParam );
  processArgsCPOPT( argc, argv, &IParam, &cfg );

  /********** memory allocation **********/
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

  
  /********** create instance data **********/  
  IloEnv env;
  IloModel model( env );
  Assignment A( env, n );

  /***** set the domain of each cell *****/
  P = (int**)malloc(n*sizeof(int*));
  for(i=0;i<n;i++)
    P[i] = (int*)malloc(n*sizeof(int));
  polygonMethod( P, n );
  for(i=0;i<n;i++){
    A[i] = IloIntVarArray( env, (IloInt)(i+1) );
    for(j=0;j<=i;j++){
      IloIntArray domain( env, n+1 );
      for(k=0;k<=n;k++)
	if( k==0 )
	  domain[k] = -(P[i][j]+1);
	else
	  domain[k] = k;
      A[i][j] = IloIntVar( env, domain );
    }
  }
  
  /***** set the constraints *****/
  /*** fixed cells & forbidden cells ***/
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++)
      if( L[i][j] != EMPTY )
	model.add( A[i][j] == L[i][j]+1 );
      else{
	for(k=0;k<n;k++)
	  if( X[i][j][k] == FALSE )
	    model.add( A[i][j] != k+1 );
      }

  /*** All-different for lines ***/
  for(i=0;i<n;i++){
    IloIntVarArray Line( env );
    for(j=0;j<n;j++)
      if( i>=j )
	Line.add( A[i][j] );
      else
	Line.add( A[j][i] );
    model.add( IloAllDiff( env, Line ) );
  }

  /***** set the objective *****/
  IloIntExpr filled( env );
  IloIntVar obj( env );
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++)
      filled += (A[i][j]>0);
  model.add( obj == filled );
  model.add( IloMaximize( env, obj ) );

  /***** solve PLSE by constraint optimization model *****/
  IloCP  cp( model );
  IloTimer timer( env );

  /*** MIP start ***/
  if( IParam.init != UNDEF ){
    IloSolution initsol( env );
    int p;
    for(p=0;p<Sol.num[NUM_SOL];p++){
      i = Sol.Pi[p]->pos[ROW];
      j = Sol.Pi[p]->pos[COL];
      k = Sol.Pi[p]->pos[VAL];
      initsol.setValue( A[i][j], k+1 );
    }
    cp.setStartingPoint( initsol );      
  }  

  /*** パラメータの設定 ***/
  /* 推論レベル */
  cp.setParameter( IloCP::DefaultInferenceLevel, IloCP::Extended );
  /* All-different 推論レベル */
  cp.setParameter( IloCP::AllDiffInferenceLevel, IloCP::Extended );
  /* verbosity レベル */
  cp.setParameter( IloCP::LogVerbosity, cfg.verbose );
  /* 時計の種類: CPU時間 */
  cp.setParameter( IloCP::TimeMode, IloCP::CPUTime );
  /* 計算時間 */
  if( IParam.timemax >= 0.0 )
    cp.setParameter( IloCP::TimeLimit, IParam.timemax );
  /* コア数 */
  cp.setParameter( IloCP::Workers, cfg.workers );

  /*** 計算 ***/
  timer.start();
  start = timer.getTime();
  cp.solve();
  printf( "CPXOutput:\t%ld\t", cp.getValue(obj) );
  if( cp.getInfo( IloCP::FailStatus ) == IloCP::SearchHasFailedNormally ||
      (int)(cp.getValue(obj)) == (n+1)*n/2 )
    printf( "+\n" );
  else
    printf( "-\n" );
  finish = timer.getTime();
  printf( "CPXCpuTime:\t%g\n", finish-start );

  /*** 解を得る ***/
  S = (int**)mallocE(n*sizeof(int*));
  for(i=0;i<n;i++){
    S[i] = (int*)mallocE(n*sizeof(int));
    for(j=0;j<=i;j++){
      if( cp.getValue( A[i][j] ) > 0 )
	S[i][j] = cp.getValue( A[i][j] ) - 1;
      else
	S[i][j] = EMPTY;
      S[j][i] = S[i][j];
    }
  }
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
  fprintf( out, "\nCPLEX CP solver for partial symmetric Latin square extension problem\n\n" );
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
  fprintf( out, "  -verbose=<INT>  verbose level: {%d,%d,%d,%d} ---> {Quiet/Terse/Normal/Verbose} (default:%d)\n",
	   IloCP::Quiet, IloCP::Terse, IloCP::Normal, IloCP::Verbose, VERBOSE );
  fprintf( out, "  -workers=<INT>  workers (default:%d)\n", WORKERS );
}


void outputErrorMsg( const char *msg ){
  fprintf( stderr, "\nError: %s\n\n", msg );
  exit( EXIT_FAILURE );
}


void processArgsCPOPT( int argc, char *argv[], ILSParam *IP, CPOPTConfig *CP ){
  char **S;
  int k,eqs;
  if( argc < 2 ){
    outputUsage( stdout, argv[0] );
    exit( EXIT_FAILURE );
  }
  CP->workers = WORKERS;
  CP->verbose = VERBOSE; 

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
	IP->q = ILSinsttype( atoi( S[1] ) );
      else if( strcmp( S[0], "r" ) == 0 )
	IP->r = atoi( S[1] );
      else if( strcmp( S[0], "hatarray" ) == 0 ){
	if( IP->hatarray != NULL )
	  free( IP->hatarray );
	IP->hatarray = (char*)malloc(LENGTH_OF_STRING*sizeof(char));
	strcpy( IP->hatarray, S[1] );
      }
      else if( strcmp( S[0], "init" ) == 0 )
	IP->init = ILSinit( atoi( S[1] ) );
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
      else if( strcmp( S[0], "workers" ) == 0 )
	CP->workers = atoi( S[1] );
      else if( strcmp( S[0], "verbose" ) == 0 )
	CP->verbose = atoi( S[1] );
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
