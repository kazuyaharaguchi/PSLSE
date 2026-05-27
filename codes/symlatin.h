/*** symlatin.h ***/

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef UNDEF
#define UNDEF -1
#endif

#ifndef EMPTY
#define EMPTY -1
#endif

#define QWH       0 
#define QC        1
#define RANDHAT   2
#define INPUTHAT  3
#define HOLERAND  4
#define HOLEHAT   5
#define HOLEFILE  6

#define AWAY 0
#define HOME 1

int generateBooleanCube( int ***X, int **L, int **H,
			 int n, int r, int seed, int type );
int generatePSLS( int **L, int n, int r, int seed, int type );
void generateSLS( int **L, int n, int seed );
void polygonMethod( int **L, int n );
int removeSLS( int **L, int n, int r, int seed );
void makeHole( int ***X, int n, int r, int seed );
void outputLatinSquare( FILE *out, int **L, int n );
int isPLS( int **L, int n );
int isSymmetric( int **L, int n );
int isExtension( int **S, int **L, int n );
int isHATCompatible( int **L, int **H, int n );
int constructSLS( int **L, int n, int r, int seed );
int isAssignable( int **L, int n, int row, int col, int val );
void makeRandHAT( int **H, int n, int seed );
void makeHATMBfromArray( int n, int **H, int *A );
void reduceByHAT( int ***X, int **H, int n );
int getNumOfForbiddenMatches( int n, int ***X, int **L );
int getNumOfFixedMatches( int n, int **L );
