/******************************
  main.c
******************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "symlatin.h"
#include "elementary.h"
#include "greedy.h"
#include "localsearch.h"
#include "misc.h"
#include "mt19937ar.h"


extern double cpu_time( void );
void solve( int input_ass, Boolean ***X, int **L,
	    ILSParam *IP, ILSRec *IR );
void solveCore( int input_ass, Boolean ***X, int **L,
		ILSParam *IP, ILSRec *IR );
void outputUsage( FILE *out, char *cmd );
void outputErrorMsg( char *msg );
void processArgs( int argc, char *argv[], ILSParam *IP );


int main( int argc, char *argv[] ){
  ILSParam IParam;
  ILSRec   IRec;
  Boolean ***X,**H=NULL;
  int **L,input_ass;
  int i,j,k;

  /*** initialize IParam ***/
  cpu_time();
  initILSParam( &IParam );
  processArgs( argc, argv, &IParam );

  /*** memory allocation ***/
  X = (Boolean***)malloc(IParam.N*sizeof(Boolean**));
  for(i=0;i<IParam.N;i++){
    X[i] = (Boolean**)malloc(IParam.N*sizeof(Boolean*));
    for(j=0;j<IParam.N;j++){
      X[i][j] = (Boolean*)malloc(IParam.N*sizeof(Boolean));
      for(k=0;k<IParam.N;k++)
	X[i][j][k] = TRUE;
    }
  }
  L = (int**)malloc(IParam.N*sizeof(int*));
  for(i=0;i<IParam.N;i++){
    L[i] = (int*)malloc(IParam.N*sizeof(int));
    for(j=0;j<IParam.N;j++)
      L[i][j] = EMPTY;
  }

  if( IParam.q == RANDHAT || IParam.q == INPUTHAT || IParam.q == HOLEHAT ){
    H = (Boolean**)malloc(IParam.n*sizeof(Boolean*));
    for(i=0;i<IParam.n;i++)
      H[i] = (Boolean*)malloc((IParam.n-1)*sizeof(Boolean));   
    if( IParam.q == RANDHAT )
      makeRandHAT( H, IParam.n, IParam.s );
    else{
      char **str;
      int *A;
      if( IParam.hatarray == NULL ){
	fprintf( stderr, "error: hatarray should be specified.\n" );
	exit( EXIT_FAILURE );
      }
      str = split( IParam.hatarray, ',' );
      A = (int*)malloc((IParam.n/2-1)*sizeof(int));
      for(i=0;i<IParam.n/2-1;i++){
	A[i] = atoi( str[i] );
	free( str[i] );
      }
      free( str );
      makeHATMBfromArray( IParam.n, H, A );
      free( A );
    }    
  }

  /*** initialize Sol ***/
  input_ass = generateBooleanCube( X, L, H, IParam.n, IParam.r, IParam.s, IParam.q );
  if( IParam.ext > 0 ){
    for(i=0;i<IParam.n;i++)
      for(j=0;j<IParam.n;j++)
	if( L[i][j] != EMPTY )
	  for(k=IParam.n;k<IParam.N;k++)
	    X[i][j][k] = FALSE;
    for(i=IParam.n;i<IParam.N;i++)
      for(j=0;j<IParam.n;j++)
	for(k=0;k<IParam.n;k++)
	  X[i][j][k] = FALSE;
    for(i=0;i<IParam.n;i++)
      for(j=IParam.n;j<IParam.N;j++)
	for(k=0;k<IParam.n;k++)
	  X[i][j][k] = FALSE;
    for(i=IParam.n;i<IParam.N;i++)
      for(j=IParam.n;j<IParam.N;j++)
	for(k=0;k<IParam.N;k++)
	  X[i][j][k] = FALSE;
  }

  /*** solve ***/
  solve( input_ass, X, L, &IParam, &IRec );
  
  /*** postprocess ***/
  printf("SolutionStatus:\t");
  if( IRec.Best->num[NUM_SOL] == IRec.Best->UB ){
    if( input_ass + IRec.Best->num[NUM_SOL] == IParam.N*IParam.N )
      printf("OPTIMAL_BY_COMPLETE_LATIN_SQUARE\n");
    else
      printf("OPTIMAL_BY_UPPER_BOUND\n");
  }
  else
    printf("NOT_OPTIMAL\n");

  printf("Input+LS:\t%d\n",input_ass+IRec.Best->num[NUM_SOL]);
  printf("LSTotalCpuTime:\t%g\n",IRec.timeIlsEnd-IRec.timeIlsStart);
  printf("LSitrtimes:\t%d\n",IRec.t);
  printf("LSCpuTimePerSingle:\t%g\n",
	 (IRec.timeIlsEnd-IRec.timeIlsStart)/(double)IRec.t);
  printf("LSitrmax:\t%d\n",IParam.itrmax);
  printf("LSnumNodes:\t%d\n",IRec.Best->Nodes);

  for(i=0;i<IRec.Best->num[NUM_SOL];i++){
    L[IRec.Best->Pi[i]->pos[ROW]][IRec.Best->Pi[i]->pos[COL]] = IRec.Best->Pi[i]->pos[VAL];
    L[IRec.Best->Pi[i]->pos[COL]][IRec.Best->Pi[i]->pos[ROW]] = IRec.Best->Pi[i]->pos[VAL];
  }
  if( !isPLS( L, IParam.N ) ){
    fprintf( stderr, "error: solution is not a PLS.\n" );
    exit( 1 );
  }
  if( !isSymmetric( L, IParam.N ) ){
    fprintf( stderr, "error: solution is not symmetric.\n" );
    exit( 1 );
  }
  if( IParam.ext == 0 && !isHATCompatible( L, H, IParam.n ) ){
    fprintf( stderr, "error: found solution is not compatible with HAT\n" );
    return 1;
  }
  if( IParam.output != NULL ){
    FILE *out;
    out = fopen( IParam.output, "w" );
    fprintf( out, "%d\n", input_ass+IRec.Best->num[NUM_SOL] );
    outputLatinSquare( out, L, IParam.N );
  }
  printf("complete.\n");

  /*** free allocated memories ***/
  freeILSRec( IParam.n, &IRec ); 
  freeILSParam( &IParam );
  for(i=0;i<IParam.N;i++){
    for(j=0;j<IParam.N;j++)
      free( X[i][j] );
    free( X[i] );
  }
  free( X );
  for(i=0;i<IParam.N;i++)
    free( L[i] );
  free( L );
  if( H != NULL ){
    for(i=0;i<IParam.N;i++)
      free( H[i] );
    free( H );
  }
    

  return 0;
}


void solve( int input_ass, Boolean ***X, int **L,
	    ILSParam *IP, ILSRec *IR ){
  double timeStart;
  timeStart = cpu_time();
  printf("SOLVE_STARTS:\n");
  
  if( IP->itrdrop == FALSE )
    solveCore( input_ass, X, L, IP, IR );
  else{
    Boolean ***Y;
    int i,j,k;
    // initialize Y
    Y = (Boolean***)malloc(IP->N*sizeof(Boolean**));
    for(i=0;i<IP->N;i++){
      Y[i] = (Boolean**)malloc(IP->N*sizeof(Boolean*));
      for(j=0;j<IP->N;j++){
	Y[i][j] = (Boolean*)malloc(IP->N*sizeof(Boolean));
	for(k=0;k<IP->N;k++)
	  Y[i][j][k] = TRUE;
      }
    }    
    for(i=0;i<IP->N;i++)
      for(j=0;j<=i;j++)
	for(k=0;k<IP->N;k++)
	  if( X[i][j][k] == FALSE && genrand_real2() < 0.5 ){
	    Y[i][j][k] = FALSE;
	    Y[j][i][k] = FALSE;
	  }
    while( 1 ){
      int numIncompatible=0,p;
      solveCore( 0, Y, L, IP, IR );
      for(p=0;p<IR->Best->num[NUM_SOL];p++){
	Node *x;
	x = IR->Best->Pi[p];
	if( X[x->pos[ROW]][x->pos[COL]][x->pos[VAL]] == FALSE ){
	  Y[x->pos[ROW]][x->pos[COL]][x->pos[VAL]] = FALSE;
	  numIncompatible++;
	}
      }
      printf("NUM_NODES:\t%d\n",IR->Best->Nodes);
      printf("NUM_INCOMPATIBLE:\t%d\n",numIncompatible);
      printf("---\n");
      if( numIncompatible == 0 )
	break;
      freeILSRec( IP->n, IR );
    }
    for(i=0;i<IP->N;i++){
      for(j=0;j<IP->N;j++)
	free( Y[i][j] );
      free( Y[i] );
    }
    free( Y );
  }

  printf("SOLVE_ENDS:\t%g\n",cpu_time()-timeStart);
}

void solveCore( int input_ass, Boolean ***X, int **L,
		ILSParam *IP, ILSRec *IR ){
  Solution Sol;
  double timeInitStart,timeInitEnd;
  int init_ass;

  printf("Input:\t%d\n",input_ass);
  initSolution( IP->N, X, L, input_ass, &Sol );
  
  /*** construct initial solution ***/
  timeInitStart = cpu_time();
  init_ass = greedy( &Sol, IP->init, IP->sinit ); 
  timeInitEnd = cpu_time();
  printf("Input+Init:\t%d\n",input_ass+init_ass);
  printf("InitCpuTime:\t%g\n",timeInitEnd-timeInitStart);
  printf("UpperBound:\t%d\n",input_ass+Sol.UB);

  /*** initialize IRec ***/
  initILSRec( IR, &Sol ); 

  /*** iterated local search ***/
  if( IR->Best->num[NUM_SOL] < IR->Best->UB ){
    IR->timeIlsStart = cpu_time();
    ILS( IP, IR, &Sol );
    IR->timeIlsEnd = cpu_time();
  }

  /*** postprocess ***/
  freeSolution( IP->n, &Sol );
}


void outputUsage( FILE *out, char *cmd ){
  fprintf( out, "\nIterated local search (ILS) algorithm\nfor partial symmetric Latin square extension problem\n\n" );
  /********************/
  fprintf( out, "Usage: %s -n=<INT> [OPTIONS]\n\n", cmd );
  fprintf( out, "  -n=<INT>        Order (i.e., grid length)\n\n" );
  fprintf( out, "Options:\n" );
  fprintf( out, "  -q=<INT>        Problem type: {0,...,6}--->{QWH,QC,RANDHAT,INPUTHAT,HOLERAND,HOLEHAT,HOLEFILE} (default:%d)\n", DEFAULT_Q );
  fprintf( out, "  -r=<INT>        Pre-assigned ratio in percent (default:%d)\n", DEFAULT_R );
  fprintf( out, "  -s=<INT>        Random seed for the randinst generator (default:1)\n\n" );
  fprintf( out, "  -alg=<INT>      ILS type: {0,1}--->{1-ILS,2-ILS} (default:%d)\n", DEFAULT_ALG );
  fprintf( out, "  -breakprob=<DOUBLE>     probability parameter in perturb (default:%g)\n", DEFAULT_BREAKPROB );
  fprintf( out, "  -datafile=<STR>\n" );
  fprintf( out, "  -hatarray=<STR>\n" );
  fprintf( out, "  -init=<INT>     Initial solution generator: {0,1,2}--->{Rand,G1,G5} (default:%d)\n", DEFAULT_INIT );
  fprintf( out, "  -itrdrop=<BOOL> Iterative drop scheme (default:%d)\n", DEFAULT_ITRDROP );
  fprintf( out, "  -itrmax=<INT>   Iteration times limit (default:infty)\n" );
  fprintf( out, "  -output=<FILE>  Write the output solution to <FILE> (default:none)\n" );
  fprintf( out, "  -seed=<INT>     Random seed for the ILS algorithm (default:%d)\n", DEFAULT_SEED );
  fprintf( out, "  -sinit=<INT>    Random seed for initial solution generator (default:1)\n" );
  fprintf( out, "  -timemax=<INT>  Time limit in sec (default:%d)\n", DEFAULT_TIMEMAX );
  fprintf( out, "  -ext=<INT>      Number of extended lines (default:%d)\n", DEFAULT_EXT );
  /********************/
}


void outputErrorMsg( char *msg ){
  fprintf( stderr, "\nError: %s\n\n", msg );
  exit( EXIT_FAILURE );
}

/********** process arguments **********/
void processArgs( int argc, char *argv[], ILSParam *IP ){
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
      else if( strcmp( S[0], "alg" ) == 0 )
	IP->alg = atoi( S[1] );
      else if( strcmp( S[0], "breakprob" ) == 0 )
	IP->breakprob = atof( S[1] );
      else if( strcmp( S[0], "datafile" ) == 0 ){
	if( IP->datafile != NULL )
	  free( IP->hatarray );
	IP->datafile = (char*)malloc(LENGTH_OF_FILENAME*sizeof(char));
	strcpy( IP->datafile, S[1] );
      }
      else if( strcmp( S[0], "hatarray" ) == 0 ){
	if( IP->hatarray != NULL )
	  free( IP->hatarray );
	IP->hatarray = (char*)malloc(LENGTH_OF_STRING*sizeof(char));
	strcpy( IP->hatarray, S[1] );
      }
      else if( strcmp( S[0], "init" ) == 0 )
	IP->init = atoi( S[1] );
      else if( strcmp( S[0], "itrdrop" ) == 0 ){
	IP->itrdrop = atoi( S[1] );
	fprintf( stderr, "\n===== CAUTION! =====\nCONSTRUCTION OF -itrdrop SCHEME IS NOT COMPLETE!\n\n\n" );
      }
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
      else if( strcmp( S[0], "seed" ) == 0 )
	IP->seed = atoi( S[1] );
      else if( strcmp( S[0], "sinit" ) == 0 )
	IP->sinit = atoi( S[1] );
      else if( strcmp( S[0], "timemax" ) == 0 )
	IP->timemax = atoi( S[1] );
      else if( strcmp( S[0], "ext" ) == 0 ){
	IP->ext = atoi( S[1] );
	fprintf( stderr, "\n===== CAUTION! =====\nCONSTRUCTION OF -ext SCHEME IS NOT COMPLETE!\n\n\n" );
      }
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
