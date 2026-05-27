/******************************
  symLSSOL.c
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

#define ENCODE_INT  0
#define ENCODE_BOOL 1

extern double cpu_time( void );
void outputUsage( FILE *out, char *cmd );
void outputErrorMsg( char *msg );
void processArgs( int argc, char *argv[], ILSParam *IP );

void outputSymLSSOLInt( FILE *out, Boolean ***X, int **L, ILSParam *IP, Solution *S );
void outputSymLSSOLBool( FILE *out, Boolean ***X, int **L, ILSParam *IP, Solution *S );


int main( int argc, char *argv[] ){
  FILE *fp;
  Solution Sol;
  ILSParam IParam;
  Boolean ***X,**H=NULL;
  double timeGrdStart,timeGrdEnd;
  int **L,input_ass,init_ass;
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
  initSolution( IParam.N, X, L, input_ass, &Sol );

  
  /*** construct initial solution ***/
  timeGrdStart = cpu_time();
  init_ass = greedy( &Sol, IParam.init, IParam.sinit ); 
  timeGrdEnd = cpu_time();
  printf("Input+Init:\t%d\n",input_ass+init_ass);
  printf("InitCpuTime:\t%g\n",timeGrdEnd-timeGrdStart);
  printf("UpperBound:\t%d\n",input_ass+Sol.UB);


  /*** output data in LocalSolver format ***/
  if( IParam.output == NULL )
    fp = stdout;
  else 
    fp = fopen( IParam.output, "w" );
  switch( IParam.alg ){
  case ENCODE_INT:
    outputSymLSSOLInt( fp, X, L, &IParam, &Sol );
    break;
  case ENCODE_BOOL:
    outputSymLSSOLBool( fp, X, L, &IParam, &Sol );
    break;
  }

  /*** free allocated memories ***/
#if 0
  freeSolution( IParam.n, &Sol );
  freeILSParam( &IParam );
#endif 
  for(i=0;i<IParam.n;i++){
    for(j=0;j<IParam.n;j++)
      free( X[i][j] );
    free( X[i] );
  }
  free( X );
  for(i=0;i<IParam.n;i++)
    free( L[i] );
  free( L );

  return 0;
}



void outputUsage( FILE *out, char *cmd ){
  fprintf( out, "\noutput problem data in LocalSolver format\nfor partial symmetric Latin square extension problem\n\n" );
  /********************/
  fprintf( out, "Usage: %s -n=<INT> [OPTIONS]\n\n", cmd );
  fprintf( out, "  -n=<INT>        Order (i.e., grid length)\n\n" );
  fprintf( out, "Options:\n" );
  fprintf( out, "  -q=<INT>        Problem type: {0,...,6}--->{QWH,QC,RANDHAT,INPUTHAT,HOLERAND,HOLEHAT,HOLEFILE} (default:%d)\n", DEFAULT_Q );
  fprintf( out, "  -r=<INT>        Pre-assigned ratio in percent (default:%d)\n", DEFAULT_R );
  fprintf( out, "  -s=<INT>        Random seed for the randinst generator (default:1)\n\n" );
  fprintf( out, "  -alg=<INT>      encoding algorithm: {0,1}--->{integer-based,bool-based} (default:%d)\n", DEFAULT_ALG );
  //fprintf( out, "  -breakprob=<DOUBLE>     probability parameter in perturn (default:%g)\n", DEFAULT_BREAKPROB );
  fprintf( out, "  -datafile=<STR>\n" );
  fprintf( out, "  -hatarray=<STR>\n" );
  fprintf( out, "  -init=<INT>     Initial solution generator: {0,1,2}--->{Rand,G1,G5} (default:%d)\n", DEFAULT_INIT );
  //fprintf( out, "  -itrdrop=<BOOL> Iterative drop scheme (default:%d)\n", DEFAULT_ITRDROP );
  //fprintf( out, "  -itrmax=<INT>   Iteration times limit (default:infty)\n" );
  fprintf( out, "  -output=<FILE>  Write the output solution to <FILE> (default:none)\n" );
  //fprintf( out, "  -seed=<INT>     Random seed for the ILS algorithm (default:%d)\n", DEFAULT_SEED );
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
      /*
      else if( strcmp( S[0], "breakprob" ) == 0 )
	IP->breakprob = atof( S[1] );
      */
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
      /*
      else if( strcmp( S[0], "itrdrop" ) == 0 ){
	IP->itrdrop = atoi( S[1] );
	fprintf( stderr, "\n===== CAUTION! =====\nCONSTRUCTION OF -itrdrop SCHEME IS NOT COMPLETE!\n\n\n" );
      }
      else if( strcmp( S[0], "itrmax" ) == 0 )
	IP->itrmax = atoi( S[1] );
      */
      else if( strcmp( S[0], "output" ) == 0 ){
	if( IP->output != NULL )
	  free( IP->output );
	IP->output = (char*)malloc(LENGTH_OF_FILENAME*sizeof(char));
	strcpy( IP->output, S[1] );
      }
      else if( strcmp( S[0], "s" ) == 0 )
	IP->s = atoi( S[1] );
      /*
      else if( strcmp( S[0], "seed" ) == 0 )
	IP->seed = atoi( S[1] );
      */
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



void outputSymLSSOLInt( FILE *out, Boolean ***X, int **L, ILSParam *IP, Solution *S ){
  Boolean available;
  Node *v;
  int i,j,k,n,ub=0;
  n = IP->n;
  /*** input ***/
  fprintf( out, "function input(){\n" );
  fprintf( out, "\tprintln(\"Input:\\t%d\");\n", S->input_ass );
  fprintf( out, "\tprintln(\"Input+Init:\\t%d\");\n}\n\n", S->input_ass+S->num[NUM_SOL] );

  /*** model ***/
  fprintf( out, "function model(){\n" );
  // declare variables
  fprintf( out, "\tx[i in 1..%d][1..i] <- int(0,%d);\n", n,n );
  // constraint: Boolean cube 
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++){
      available = FALSE;
      for(k=0;k<n;k++)
	if( X[i][j][k] )
	  available = TRUE;
	else 
	  fprintf( out, "\tconstraint x[%d][%d] != %d;\n", i+1, j+1, k+1 );
      if( available )
	ub++;
    }
  // constraint: pairwise-different
  fprintf( out, "\tfor [i in 1..%d]{\n", n );
  fprintf( out, "\t\tfor [j in 1..%d]{\n", n-1 );
  fprintf( out, "\t\t\tif ( i >= j ){\n" );
  fprintf( out, "\t\t\t\ta <- x[i][j];\n\t\t\t}\n" );
  fprintf( out, "\t\t\telse{\n" );
  fprintf( out, "\t\t\t\ta <- x[j][i];\n\t\t\t}\n" );
  fprintf( out, "\t\t\tfor [k in (j+1)..%d]{\n", n );
  fprintf( out, "\t\t\t\tif ( i >= k ){\n" );
  fprintf( out, "\t\t\t\t\tb <- x[i][k];\n\t\t\t\t}\n" );
  fprintf( out, "\t\t\t\telse{\n" );
  fprintf( out, "\t\t\t\t\tb <- x[k][i];\n\t\t\t\t}\n" );
  fprintf( out, "\t\t\t\tconstraint a == 0 || b == 0 || a != b;\n" );
  fprintf( out, "\t\t\t}\n" );
  fprintf( out, "\t\t}\n" );
  fprintf( out, "\t}\n" );
  // objective
  fprintf( out, "\tfilled[i in 1..%d][j in 1..i] <- x[i][j]>0 ? 1:0;\n", n );
  fprintf( out, "\tobj <- sum [i in 1..%d][j in 1..i](filled[i][j]);\n", n );
  fprintf( out, "\tconstraint obj <= %d;\n", ub );
  fprintf( out, "\tmaximize obj;\n}\n\n" );

  /*** param ***/ 
  fprintf( out, "function param(){\n" );
  fprintf( out, "\tlsTimeLimit=%d;\n", IP->timemax );
  fprintf( out, "\tlsVerbosity=2;\n" );
  fprintf( out, "\tsetObjectiveBound(0,%d);\n",ub );
  for(i=0;i<S->num[NUM_SOL];i++){
    v = S->Pi[i];
    fprintf( out, "\tsetValue(x[%d][%d],%d);\n",
	     v->pos[ROW]+1, v->pos[COL]+1, v->pos[VAL]+1 );
  }
  fprintf( out, "}\n\n" );

  /*** display ***/
  fprintf( out, "function display(){\n" );
  fprintf( out, "}\n\n" );

  /*** output ***/
  fprintf( out, "function output(){\n" );
  fprintf( out, "\tprintln(\"LSSOLOutput:\\t\"+(getValue(obj)+%d));\n", S->input_ass );
  fprintf( out, "\tprintln(\"LSSOLCpuTime:\\t\"+lsTimeLimit);\n" );
  fprintf( out, "}\n\n" );
  
}



void outputSymLSSOLBool( FILE *out, Boolean ***X, int **L, ILSParam *IP, Solution *S ){ 
  Boolean available;
  Node *v;
  int i,j,k,n,ub=0;
  n = IP->n;
  /*** input ***/
  fprintf( out, "function input(){\n" );
  fprintf( out, "\tprintln(\"Input:\\t%d\");\n", S->input_ass );
  fprintf( out, "\tprintln(\"Input+Init:\\t%d\");\n}\n\n", S->input_ass+S->num[NUM_SOL] );

  /*** model ***/
  fprintf( out, "function model(){\n" );
  // declare variables
  fprintf( out, "\tx[i in 1..%d][1..i][1..%d] <- bool();\n", n,n );
  // constraint: Boolean cube 
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++){
      available = FALSE;
      for(k=0;k<n;k++)
	if( X[i][j][k] )
	  available = TRUE;
	else 
	  fprintf( out, "\tconstraint x[%d][%d][%d] == false;\n", i+1, j+1, k+1 );
      if( available )
	ub++;
    }
  // constraint: v-edge
  fprintf( out, "\tfor [i in 1..%d]{\n", n );
  fprintf( out, "\t\tfor [j in 1..i]{\n" );
  fprintf( out, "\t\t\ta <- sum[ k in 1..%d ] (x[i][j][k]);\n", n );
  fprintf( out, "\t\t\tconstraint a <= 1;\n" );
  fprintf( out, "\t\t}\n" );
  fprintf( out, "\t}\n" );
  // constraint: h-edge
  fprintf( out, "\tfor [k in 1..%d]{\n", n );
  fprintf( out, "\t\tfor [i in 1..%d]{\n", n );
  fprintf( out, "\t\t\ta <- sum[ j in 1..i ] (x[i][j][k]);\n" );
  fprintf( out, "\t\t\tb <- sum[ j in (i+1)..%d ] (x[j][i][k]);\n", n );
  fprintf( out, "\t\t\tconstraint a+b <= 1;\n" );
  fprintf( out, "\t\t}\n" );
  fprintf( out, "\t}\n" );
  // objective
  fprintf( out, "\tobj <- sum [i in 1..%d][j in 1..i][k in 1..%d] (x[i][j][k]);\n",n,n );
  fprintf( out, "\tconstraint obj <= %d;\n", ub );
  fprintf( out, "\tmaximize obj;\n}\n\n" );

  /*** param ***/ 
  fprintf( out, "function param(){\n" );
  fprintf( out, "\tlsTimeLimit=%d;\n", IP->timemax );
  fprintf( out, "\tlsVerbosity=2;\n" );
  fprintf( out, "\tsetObjectiveBound(0,%d);\n",ub );
  for(i=0;i<S->num[NUM_SOL];i++){
    v = S->Pi[i];
    fprintf( out, "\tsetValue(x[%d][%d][%d],true);\n",
	     v->pos[ROW]+1, v->pos[COL]+1, v->pos[VAL]+1 );
  }
  fprintf( out, "}\n\n" );

  /*** display ***/
  fprintf( out, "function display(){\n" );
  fprintf( out, "}\n\n" );

  /*** output ***/
  fprintf( out, "function output(){\n" );
  fprintf( out, "\tprintln(\"LSSOLOutput:\\t\"+(getValue(obj)+%d));\n", S->input_ass );
  fprintf( out, "\tprintln(\"LSSOLCpuTime:\\t\"+lsTimeLimit);\n" );
  fprintf( out, "}\n\n" ); 
}
