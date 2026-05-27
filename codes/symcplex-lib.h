/********** symcplex-lib.h **********/

typedef int Boolean;

#define TIMELIMIT        10.00
#ifndef LENGTH_OF_STRING
#define LENGTH_OF_STRING 1024
#endif
#define DEBUG
#undef DEBUG

struct _CPXInstance{
  char *probname;
  int numForbidden;
  int numFixed;
  int numNonZero;
  int numCols;
  int numRows;
  int objSense;      // 1:minimize, -1:maximize 
  double *c;         // coefficients of objective
  double *b;         // righthand side values of constraint
  char *sense;       // sense of each constraint: 'L' as <=, 'E' as ==, 'R' as >=
  int *matbeg;    // the starting index of each column in const matrix
  int *matcnt;    // the number of non-zero coefficients of each column in const matrix
  int *matind;    // row indices of non-zero coefficients in const matrix
  double *matval;    // non-zero coefficients in const matrix
  double *LB;        // lower bound of the variables
  double *UB;        // upper bound of the variables
  double *Range;     // bounds for ranged constraints (ignored)
  char   *ctype;     // variable type
};
typedef struct _CPXInstance *CPXInstance;

int Nall;

int initCPXSymInstance( CPXInstance I, Boolean ***X, int **L, int n );
int setCPXSymInstance( CPXInstance I, Boolean ***X, int **L, int n );
void getSolutionFromCPXInstance( CPXENVptr env, CPXLPptr IP, CPXInstance I, int **S, int n );
void getNextCell( int n, int *ip, int *jp, int *kp );
void freeCPXSymInstance( CPXInstance I );
