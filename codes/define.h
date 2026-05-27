/******************************
  define.h
******************************/

#define DEBUG
#undef DEBUG

/***** macros *****/
#define LENGTH_OF_FILENAME 256
#define LENGTH_OF_STRING   1024

#define DEFAULT_Q       QWH
#define DEFAULT_R       0
#define DEFAULT_ALG     ILS2
#define DEFAULT_INIT    Rand
#define DEFAULT_ITRMAX  -1
#define DEFAULT_SEED    1
#define DEFAULT_TIMEMAX 10
#define DEFAULT_BREAKPROB 0.333333
#define DEFAULT_ITRDROP FALSE
#define DEFAULT_EXT     0

#define DIM 3
#define ROW 0
#define COL 1
#define VAL 2

#define DIR 2
#define HOR 0
#define VER 1

#define LEFT  0
#define RIGHT 1

#define UNDEF -1
#define EMPTY -1
#define TRUE 1
#define FALSE 0

#define SECTION    5
#define NUM_SOL    0
#define NUM_FREE   1
#define NUM_1TIGHT 2
#define NUM_2TIGHT 3
#define NUM_3TIGHT 4

#define FORWARD 0
#define BACKWARD 1

#define BI 2


/***** enum's *****/
enum ILSinsttype  {QWH,QC,RANDHAT,INPUTHAT,HOLERAND,HOLEHAT,HOLEFILE};
enum ILSalg       {ILS1,ILS2};
enum ILSinit      {Rand,G1,G5};

/***** typedefs & structs *****/
typedef int Boolean;
typedef int Dim;
typedef int Dir;

typedef struct{
  int index;    // index in heap (diff from Node.index)
  int deg; 
  Boolean flag;
} GreedyInfo;


typedef struct _Node Node;
struct _Node{
  int index;                    // node index
  int order;                    // ordering
  int pos[DIM];                 // coordinate
  int deg;                      // degree
  struct _Node *next[DIM];      // next neighbor
  struct _Node *prev[DIM];      // prev neighbor
  struct _Node *solneigh[DIM];  // sol neighbor
  Boolean isSol;            // flag for sol node
  int tightness;            // tightness (for non-sol node) 
  int penalty;          // penalty
  GreedyInfo gi;            // info for greedy algorithm 
};


typedef struct{
  int n;                // order (i.e., grid length)
  Boolean ***X;         // given Boolean cube
  int **L;              // given PLS
  int input_ass;        // # of given symbols
  int UB;               // upper bound on num[NUM_SOL]
  Node ****C;           // 3D array of nodes
  int **num1t[DIR];     // # of 1-tight nodes in a grid line
  Node ***first1t[DIR]; // pointer to first 1-tight node among a grid line
  int Nodes;            // # of all nodes
  int num[SECTION];     // # of nodes in each section of ordering
  Node **Pi;            // node ordering
} Solution;


typedef struct{
  int               n;
  enum ILSinsttype  q;
  int               r;
  enum ILSalg       alg;
  char              *datafile;
  char              *hatarray;
  enum ILSinit      init;
  int               itrmax;
  char              *output;
  int               s;
  int               seed;
  int               sinit;
  int               timemax;
  double            breakprob;
  int               itrdrop;
  int               ext;
  int               N;       // N=n+ext
} ILSParam;


typedef struct{
  double timeIlsStart;
  double timeIlsEnd;
  int t;
  Solution *Best;
  unsigned int randptr;
  int lastImproved;
} ILSRec;

