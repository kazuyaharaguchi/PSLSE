/******************************
  greedy.c
******************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "define.h"
#include "greedy.h"
#include "elementary.h"
#include "misc.h"
#include "mt19937ar.h"

/***** wrapper for constructive algorithms *****/
int greedy( Solution *S, enum ILSinit init, int seed ){
  int newly_ass;
  switch( init ){
  case Rand:
    newly_ass = greedyRandom( S, seed );
    break;
  case G1:
    newly_ass = greedyMinDeg( S, G1, seed );
    break;
  case G5:
    newly_ass = greedyMinDeg( S, G5, seed );
    break;
  default:
    fprintf( stderr, "error: -init=%d is inappropriate.\n",init);
    exit( EXIT_FAILURE );
  }

#ifdef DEBUG
  /*** check whether data is stored appropriately ***/
  /*
  Dim root_d,d;
  Dir dir;
  int i,j,t;
  Node *x,*v;
  for(i=0;i<S->num[NUM_SOL];i++){
    x = S->Pi[i];
    printf("SolNode [%d,%d,%d]:\n",x->pos[ROW],x->pos[COL],x->pos[VAL]);
    root_d = 0;
    d = 0;
    dir = FORWARD;
    v = x;
    while( 1 ){
      v = getNextNeighbor( x, v, &root_d, &d, &dir );
      if( v == NULL )
	break;
      printf("\t[%d,%d,%d] ... tightness=%d",
	     v->pos[ROW],v->pos[COL],v->pos[VAL],v->tightness);
      for(t=0;t<DIM;t++)
	if( v->solneigh[t] != NULL ){
	  Node *y = v->solneigh[t];
	  printf(" %d:[%d,%d,%d]",t,y->pos[ROW],y->pos[COL],y->pos[VAL]);
	}
      printf("\n");
    }
  }
  printf("--- line-val ---\n");
  for(i=0;i<S->n;i++)
    for(j=0;j<S->n;j++)
      if( ( S->num1t[HOR][i][j] == 0 &&
	    S->first1t[HOR][i][j] != NULL ) ||
	  ( S->num1t[HOR][i][j] != 0 &&
	    S->first1t[HOR][i][j] == NULL ) ){
	fprintf( stderr, "num1t or first1t is used inappropriately.\n" );
	exit( EXIT_FAILURE );
      }
      else if( S->num1t[HOR][i][j] == 0 &&
	       S->first1t[HOR][i][j] == NULL )
	continue;
      else{
	v = S->first1t[HOR][i][j];
	printf("  (%d,%d) ... num1t=%d, first1t=[%d,%d,%d] with isSol=%d, t=%d\n",
	       i,j,S->num1t[HOR][i][j],
	       v->pos[ROW],v->pos[COL],v->pos[VAL],v->isSol,v->tightness);
      }
  printf("--- row-col ---\n");
  for(i=0;i<S->n;i++)
    for(j=0;j<S->n;j++)
      if( ( S->num1t[VER][i][j] == 0 &&
	    S->first1t[VER][i][j] != NULL ) ||
	  ( S->num1t[VER][i][j] != 0 &&
	    S->first1t[VER][i][j] == NULL ) ){
	fprintf( stderr, "num1t or first1t is used inappropriately.\n" );
	exit( EXIT_FAILURE );
      }
      else if( S->num1t[VER][i][j] == 0 &&
	       S->first1t[VER][i][j] == NULL )
	continue;
      else{
	v = S->first1t[VER][i][j];
	printf("  (%d,%d) ... num1t=%d, first1t=[%d,%d,%d] with isSol=%d, t=%d\n",
	       i,j,S->num1t[VER][i][j],
	       v->pos[ROW],v->pos[COL],v->pos[VAL],v->isSol,v->tightness);
      }
  */
#endif

  return newly_ass;
}


/***** random insertion *****/
int greedyRandom( Solution *S, int seed ){
  int p,newly_ass=0;
  if( seed > 0 )
    init_genrand( seed );
  while( S->num[NUM_FREE] > 0 ){
    p = (int)( genrand_real2() * S->num[NUM_FREE] );
    insert( S, S->num[NUM_SOL]+p );
    newly_ass++;
  }
  return newly_ass;
}


/***** minimum degree algorithm *****/
int greedyMinDeg( Solution *S, enum ILSinit init, int seed ){
  Node **H,**Min,*v,*w,*u;
  Dim root_d,root_d_,d,d_;
  Dir dir,dir_;
  int **Deg[DIR],tail,heapsize=0,minsize,newly_ass=0;
  int edgestar,maxedgestar;
  int i,j,k,p,q;
  if( seed > 0 )
    init_genrand( seed );
  /*** preprocess ***/
  tail = S->num[NUM_SOL] + S->num[NUM_FREE];
  for(p=S->num[NUM_SOL];p<tail;p++){
    S->Pi[p]->gi.deg = 0;
    S->Pi[p]->gi.flag = FALSE;
  }

  /*** compute the degree ***/
  for(dir=0;dir<DIR;dir++){
    Deg[dir] = (int**)malloc(S->n*sizeof(int*));
    for(i=0;i<S->n;i++){
      Deg[dir][i] = (int*)malloc(S->n*sizeof(int));
      for(j=0;j<S->n;j++)
	Deg[dir][i][j] = 0;
    }
  }
  // horizontal direction
  for(k=0;k<S->n;k++)
    for(p=0;p<S->n;p++){
      v = NULL;
      for(q=0;q<S->n;q++)
	if( q <= p && S->C[p][q][k] != NULL ){
	  v = S->C[p][q][k];
	  break;
	}
	else if( q > p && S->C[q][p][k] != NULL ){
	  v = S->C[q][p][k];
	  break;
	}
      if( v == NULL )
	continue;
      if( p == q )
	root_d = ROW;
      else
	root_d = COL;
      d = root_d;
      dir = FORWARD;
      w = v;
      while( 1 ){
	Deg[HOR][k][p]++;
	w = getNextNeighborDim( v, w, root_d, &d, &dir );
	if( w == NULL )
	  break;
      }
    }
  // vertical direction
  for(p=0;p<S->n;p++)
    for(q=0;q<=p;q++){
      v = NULL;
      for(k=0;k<S->n;k++)
	if( S->C[p][q][k] != NULL ){
	  v = S->C[p][q][k];
	  break;
	}
      if( v == NULL )
	continue;
      root_d = VAL;
      d = root_d;
      dir = FORWARD;
      w = v;
      while( 1 ){
	Deg[VER][p][q]++;
	w = getNextNeighborDim( v, w, root_d, &d, &dir );
	if( w == NULL )
	  break;
      }
    }
  // count the degree
  for(i=0;i<S->n;i++)
    for(j=0;j<=i;j++)
      for(k=0;k<S->n;k++){
	if( S->C[i][j][k] == NULL )
	  continue;
	v = S->C[i][j][k];
	if( i == j )
	  v->gi.deg = Deg[HOR][k][i] + Deg[VER][i][j] - 2;
	else
	  v->gi.deg = Deg[HOR][k][i] + Deg[HOR][k][j] + Deg[VER][i][j] - 3;
      }
  for(dir=0;dir<DIR;dir++){
    for(i=0;i<S->n;i++)
      free( Deg[dir][i] );
    free( Deg[dir] );
  }
  
  /*** insert free nodes into the heap ***/
  H = (Node**)malloc(S->Nodes*sizeof(Node*));
  for(p=S->num[NUM_SOL];p<tail;p++)
    heapsize = insertHeap( H, heapsize, S->Pi[p] );
  Min = (Node**)malloc(heapsize*sizeof(Node*));
  
  /*** iteration ***/
  while( S->num[NUM_FREE] > 0 ){
    /* collect minimum degree nodes */
    minsize = 0;
    collectMin( H, Min, heapsize, &minsize, 0 );
    /* pick up one from Min */
    if( init == G1 )
      p = Min[(int)(genrand_real2()*(double)minsize)]->gi.index;
    else{
      if( minsize == 1 )
	p = Min[0]->gi.index;
      else{
	maxedgestar = -INT_MAX;
	shuffle( Min, minsize, sizeof(Node*), genrand_real2 );
	for(q=0;q<minsize;q++){
	  edgestar = getEdgeStar( S->C, Min[q] );
	  if( edgestar > maxedgestar ){
	    p = Min[q]->gi.index;
	    maxedgestar = edgestar;
	  }
	}
      }
    }
    /* drop the neighbors and
       update the degrees of their neighbors */
    v = H[p];
    root_d = ROW;
    d = root_d;
    dir = FORWARD;
    w = v;
    while( 1 ){
      w = getNextNeighbor( v, w, &root_d, &d, &dir );
      if( w == NULL )
	break;
      if( w->isSol == TRUE || w->tightness != 0 )
	continue;
      // update the degrees of w'-s neighbors
      root_d_ = ROW;
      d_ = root_d_;
      dir_ = FORWARD;
      u = w;
      while( 1 ){
	u = getNextNeighbor( w, u, &root_d_, &d_, &dir_ );
	if( u == NULL )
	  break;
	if( root_d_ == root_d || u->isSol == TRUE || u->tightness != 0 )
	  continue;
	(u->gi.deg)--;
	if( u->gi.index >= 0 && u->gi.index < heapsize ) // この1行たぶん f8-public にも必要
	  ascendHeap( H, u->gi.index );
      }
      // drop w from the heap
      heapsize = deleteHeap( H, heapsize, w->gi.index );
    }
    /* insert the picked node into the solution */
    insert( S, v->order );
    heapsize = deleteHeap( H, heapsize, v->gi.index );
    newly_ass++;
#ifdef DEBUG
    for(i=0;i<heapsize;i++)
      if( i != H[i]->gi.index ){
	fprintf(stderr,"error: heap is broken.\n");
	exit( EXIT_FAILURE );
      }
#endif
  }
  free( H );
  free( Min );
  return newly_ass;
}

void collectMin( Node **H, Node **Min, int heapsize, int *minp, int i ){
  int child[2],c;
  Min[*minp] = H[i];
  (*minp)++;
  child[0] = 2*i+1;
  child[1] = 2*i+2;
  for(c=0;c<2;c++){
    if( child[c] >= heapsize )
      continue;
    if( H[child[c]]->gi.deg > H[0]->gi.deg )
      continue;
    collectMin( H, Min, heapsize, minp, child[c] );
  }
}


int insertHeap( Node **H, int heapsize, Node *v ){
  H[heapsize] = v;
  H[heapsize]->gi.index = heapsize;
  ascendHeap( H, heapsize );
  return heapsize+1;
}


int deleteHeap( Node **H, int heapsize, int k ){
  H[heapsize-1]->gi.index = k;
  swap( H, k, heapsize-1, sizeof(Node*) );
  descendHeap( H, heapsize-1, k );
  return heapsize-1;  
}


void ascendHeap( Node **H, int k ){
  int parent;
  if( k == 0 )
    return;
  parent = (k-1)/2;
  if( H[parent]->gi.deg > H[k]->gi.deg ){
    H[parent]->gi.index = k;
    H[k]->gi.index = parent;
    swap( H, parent, k, sizeof(Node*) );
    ascendHeap( H, parent );
  }
}

void descendHeap( Node **H, int heapsize, int k ){
  int child;
  if( 2*k+1 >= heapsize )
    return;
  else if( 2*k+2 == heapsize )
    child = 2*k+1;
  else{
    if( H[2*k+1]->gi.deg < H[2*k+2]->gi.deg )
      child = 2*k+1;
    else
      child = 2*k+2;
  }
  if( H[k]->gi.deg > H[child]->gi.deg ){
    H[k]->gi.index = child;
    H[child]->gi.index = k;
    swap( H, k, child, sizeof(Node*) );
    descendHeap( H, heapsize, child );
  }
}


int getEdgeStar( Node ****C, Node *v ){
  Dim root_d,root_d_,d,d_;
  Dir dir,dir_;
  Node *w,*u;
  int edgestar=0;
  /*** clear the flag ***/
  root_d = 0;
  d = 0;
  dir = FORWARD;
  w = v;
  while( 1 ){
    w = getNextNeighbor( v, w, &root_d, &d, &dir );
    if( w == NULL )
      break;
    if( w->isSol == TRUE || w->tightness != 0 )
      continue;
    /* search w'-s neighbors */
    root_d_ = 0;
    d_ = 0;
    dir_ = FORWARD;
    u = w;
    while( 1 ){
      u = getNextNeighbor( w, u, &root_d_, &d_, &dir_ );
      if( u == NULL )
	break;
      if( d_ == d )
	continue;
      if( u->isSol == TRUE || u->tightness != 0 )
	continue;
      u->gi.flag = FALSE;
    }
  }
  /*** compute the edgestar ***/
  root_d = 0;
  d = 0;
  dir = FORWARD;
  w = v;
  while( 1 ){
    w = getNextNeighbor( v, w, &root_d, &d, &dir );
    if( w == NULL )
      break;
    if( w->isSol == TRUE || w->tightness != 0 )
      continue;
    /* search w'-s neighbors */
    root_d_ = 0;
    d_ = 0;
    dir_ = FORWARD;
    u = w;
    while( 1 ){
      u = getNextNeighbor( w, u, &root_d_, &d_, &dir_ );
      if( u == NULL )
	break;
      if( d_ == d )
	continue;
      if( u->isSol == TRUE || u->tightness != 0 )
	continue;
      if( u->gi.flag == FALSE ){
	edgestar++;
	u->gi.flag = TRUE;
      }
    }
  }
  return edgestar;
}
