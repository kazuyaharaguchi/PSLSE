/******************************
  elementary.h
******************************/
void initILSParam( ILSParam *IP );
void freeILSParam( ILSParam *IP );
void initSolution( int n, Boolean ***X, int **L, int input_ass, Solution *S );
int initNodeCube( int n, Boolean ***X, int **L, Node ****C, int *UB_ptr );
void initNodeOrdering( int n, Node ****C, Node **Pi ); 
void freeSolution( int n, Solution *S );
void initILSRec( ILSRec *IR, Solution *S );
void freeILSRec( int n, ILSRec *IR );
void insert( Solution *S, int p );
void drop( Solution *S, int p );
void exchange( Node **Pi, int p, int q );
void in1Tight( Solution *S, int p );
void out1Tight( Solution *S, int p );
Node *getNextNeighbor( Node *root, Node *v, Dim *root_dp, Dim *dp, Dir *dirp );
Node *getNextNeighborDim( Node *root, Node *v, Dim root_d, Dim *dp, Dir *dirp );
int get1TightNeighbor( Node *v, Node **N, int start, Dim d );
Boolean isAdjacent( Node *v, Node *w );
Dim getDistance( Node *v, Node *w );
Node *getMostPenalizedNode( Node **N, int numN );
Node *getCommon2TightNeighbor( Node ****C, Node *u, Node *x, Node *y, Dim dx, Dim dy );
Dim getDiffPos( Node *u, Node *v );
Dim getSamePos( Node *u, Node *v );
Dir getOrthogonalDir( Node *u, Node *v );
int cmpByIndex( const void *p, const void *q );
void getSolutionDiff( Solution *S, Solution *U, Node **SolS, Node **SolU,
		      Node **S_U, Node **U_S, int *nS_Uptr, int *nU_Sptr );
void exchange1TightRandomly( Solution *S );
void exchange2TightRandomly( Solution *S );
