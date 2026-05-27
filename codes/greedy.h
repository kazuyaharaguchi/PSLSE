/******************************

  Iterated local search for
  partial Latin square extension
  by Kazuya Haraguchi (2015)

  greedy.h

  header file for greedy.c 

******************************/

int greedy( Solution *S, enum ILSinit init, int seed );
int greedyRandom( Solution *S, int seed );
int greedyMinDeg( Solution *S, enum ILSinit init, int seed );
void collectMin( Node **H, Node **Min, int heapsize, int *minp, int i );
int insertHeap( Node **H, int heapsize, Node *v );
int deleteHeap( Node **H, int heapsize, int k );
void ascendHeap( Node **H, int k );
void descendHeap( Node **H, int heapsize, int k );
int getEdgeStar( Node ****C, Node *v );
