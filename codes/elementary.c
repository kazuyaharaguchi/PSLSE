/******************************
  elementary.c
******************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "define.h"
#include "elementary.h"
#include "misc.h"
#include "mt19937ar.h"

/***** initialize ILSParam: set default parameter values *****/
void initILSParam( ILSParam *IP ){
  IP->q = DEFAULT_Q;
  IP->r = DEFAULT_R;
  IP->alg = DEFAULT_ALG;
  IP->breakprob = DEFAULT_BREAKPROB;
  IP->datafile = NULL;
  IP->hatarray = NULL;
  IP->init = DEFAULT_INIT;
  IP->itrdrop = DEFAULT_ITRDROP;
  IP->itrmax = DEFAULT_ITRMAX;
  IP->output = NULL;
  IP->s = 1;
  IP->seed = 1;
  IP->sinit = 1;
  IP->timemax = DEFAULT_TIMEMAX;
  IP->ext = DEFAULT_EXT;
}


/***** free ILSParam *****/
void freeILSParam( ILSParam *IP ){
  free( IP->output );
}


/***** initialize solution *****/
void initSolution( int n, Boolean ***X, int **L, int input_ass, Solution *S ){
  Dir dir;
  int s,k,l;
  S->n = n;
  S->X = X; 
  S->L = L;
  S->input_ass = input_ass;
  S->UB = n*(n+1)/2-input_ass;
  S->C = (Node****)malloc(n*sizeof(Node***));
  for(dir=0;dir<DIR;dir++){
    S->num1t[dir] = (int**)malloc(n*sizeof(int*));
    S->first1t[dir] = (Node***)malloc(n*sizeof(Node**));
    for(k=0;k<n;k++){
      S->num1t[dir][k] = (int*)malloc(n*sizeof(int));
      S->first1t[dir][k] = (Node**)malloc(n*sizeof(Node*));
      for(l=0;l<n;l++){
	S->num1t[dir][k][l] = 0;
	S->first1t[dir][k][l] = NULL;
      }
    }
  }
  S->Nodes = initNodeCube( n, X, L, S->C, &(S->UB) );
  for(s=0;s<SECTION;s++)
    S->num[s] = 0;
  S->num[NUM_FREE] = S->Nodes;
  S->Pi = (Node**)malloc(S->Nodes*sizeof(Node*));
  initNodeOrdering( n, S->C, S->Pi );
}


/***** initialize 3D array of nodes *****/
int initNodeCube( int n, Boolean ***X, int **L, Node ****C, int *UB_ptr ){
  Boolean flag;
  Node *prev;
  Dim d;
  int i,j,k,p,q,nodes=0;
  
  /*** compute C ***/
  for(i=0;i<n;i++){
    C[i] = (Node***)malloc((i+1)*sizeof(Node**));
    for(j=0;j<=i;j++){
      C[i][j] = (Node**)malloc(n*sizeof(Node*));
      for(k=0;k<n;k++){
	if( X[i][j][k] == FALSE ){
	  C[i][j][k] = NULL;
	  continue;
	}
	C[i][j][k] = (Node*)malloc(sizeof(struct _Node));
	C[i][j][k]->index = nodes;
	C[i][j][k]->order = UNDEF;
	C[i][j][k]->pos[ROW] = i;
	C[i][j][k]->pos[COL] = j;
	C[i][j][k]->pos[VAL] = k;
	C[i][j][k]->deg = 0;
	for(d=0;d<DIM;d++){
	  C[i][j][k]->next[d] = NULL;
	  C[i][j][k]->prev[d] = NULL;
	  C[i][j][k]->solneigh[d] = NULL;
	}
	C[i][j][k]->isSol = FALSE;
	C[i][j][k]->tightness = 0;
	C[i][j][k]->penalty = 0;
	nodes++;
      }
    }
  }

  /*** connect adjacent nodes ***/
  // connection among same-value plane
  for(k=0;k<n;k++)
    for(p=0;p<n;p++){
      prev = NULL;
      for(q=0;q<n;q++){
	if( q<p ){ i=p; j=q; }
	else{ i=q; j=p; }
       	if( C[i][j][k] == NULL )
	  continue;
	C[i][j][k]->prev[ C[i][j][k]->deg ] = prev;
	C[i][j][k]->next[ C[i][j][k]->deg ] = NULL;
	( C[i][j][k]->deg )++;
	if( prev != NULL )
	  prev->next[ prev->deg-1 ] = C[i][j][k];
	prev = C[i][j][k];
      }
    }
  // connection between different planes 
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++){
      prev = NULL;
      for(k=0;k<n;k++){
       	if( C[i][j][k] == NULL )
	  continue;
	C[i][j][k]->prev[VAL] = prev;
	C[i][j][k]->next[VAL] = NULL;
	if( prev != NULL )
	  prev->next[VAL] = C[i][j][k];
	prev = C[i][j][k];
      }
    }

  /*** compute UB if specified ***/
  if( UB_ptr != NULL ){
    *UB_ptr = n*(n+1)/2;
    for(i=0;i<n;i++)
      for(j=0;j<=i;j++){
	flag = TRUE;
	for(k=0;k<n;k++)
	  if( C[i][j][k] != NULL ){
	    flag = FALSE;
	    break;
	  }
	if( flag )
	  (*UB_ptr)--;
      }
  }

  return nodes;
}


/***** initialize node ordering *****/
void initNodeOrdering( int n, Node ****C, Node **Pi ){
  Node *v;
  int i,j,k,p=0;
  for(i=0;i<n;i++)
    for(j=0;j<=i;j++){
      v = NULL;
      for(k=0;k<n;k++)
	if( C[i][j][k] != NULL ){
	  v = C[i][j][k];
	  break;
	}
      while( v != NULL ){
	Pi[p] = v;
	Pi[p]->order = p;
	p++;
	v = v->next[VAL];
      }
    }
}

/***** free Solution *****/
void freeSolution( int n, Solution *S ){
  Dir dir;
  int i,j,k;
  for(i=0;i<n;i++){
    for(j=0;j<=i;j++){
      for(k=0;k<n;k++)
	if( S->C[i][j][k] != NULL )
	  free( S->C[i][j][k] );
      free( S->C[i][j] );
    }
    free( S->C[i] );
  }
  free( S->C );
  for(dir=0;dir<DIR;dir++){
    for(i=0;i<n;i++){
      free( S->num1t[dir][i] );
      free( S->first1t[dir][i] );
    }
    free( S->num1t[dir] );
    free( S->first1t[dir] );
  }
  free( S->Pi );
}


/***** initialize ILSRec *****/
void initILSRec( ILSRec *IR, Solution *S ){
  Node *v;
  Dim d;
  int p,pos[DIM];
  IR->timeIlsStart = 0.0;
  IR->timeIlsEnd = 0.0;
  IR->t = 1; // t=0 corresponds to initial solution step
  IR->Best = (Solution*)malloc(sizeof(Solution));
  initSolution( S->n, S->X, S->L, S->input_ass, IR->Best );
  for(p=0;p<S->num[NUM_SOL];p++){
    for(d=0;d<DIM;d++)
      pos[d] = S->Pi[p]->pos[d];
    v = IR->Best->C[pos[ROW]][pos[COL]][pos[VAL]];
    insert( IR->Best, v->order );
  }
  IR->randptr = genrand_int31() % S->Nodes; 
  IR->lastImproved = 0;
}


/***** free ILSRec *****/
void freeILSRec( int n, ILSRec *IR ){
  freeSolution( n, IR->Best );
  free( IR->Best );
}


/***** insert a node into the solution *****/
void insert( Solution *S, int p ){
  Dim root_d,d;
  Dir dir;
  Node *x,*v;
#ifdef DEBUG
  if( p < S->num[NUM_SOL] || p >= S->num[NUM_SOL]+S->num[NUM_FREE] ){
    fprintf( stderr, "error: you tried to insert non-free node: " );
    fprintf( stderr, "p=%d (at [%d,%d,%d] isSol=%d, tightness=%d), #sol=%d, #free=%d\n",
	     p, S->Pi[p]->pos[ROW], S->Pi[p]->pos[COL], S->Pi[p]->pos[VAL],  
	     S->Pi[p]->isSol, S->Pi[p]->tightness,
	     S->num[NUM_SOL], S->num[NUM_FREE] );
    exit( EXIT_FAILURE );
  }    
#endif
  /*** exchange ***/
  x = S->Pi[p];
  exchange( S->Pi, p, S->num[NUM_SOL] );
  x->isSol = TRUE;
  x->tightness = 0;
  (S->num[NUM_SOL])++;
  (S->num[NUM_FREE])--;
  /*** process neighbors ***/
  root_d = ROW;
  d = root_d;
  dir = FORWARD;
  v = x;
  while( 1 ){
    v = getNextNeighbor( x, v, &root_d, &d, &dir );
    if( v == NULL )
      break;
    v->solneigh[d] = x;
    (v->tightness)++;
    if( v->tightness == 1 ){
      in1Tight( S, v->order );
      exchange( S->Pi, v->order, S->num[NUM_SOL]+S->num[NUM_FREE]-1 );
      (S->num[NUM_FREE])--;
      (S->num[NUM_1TIGHT])++;
    }
    else if( v->tightness == 2 ){
      out1Tight( S, v->order );
      exchange( S->Pi, v->order, S->num[NUM_SOL]+S->num[NUM_FREE]+S->num[NUM_1TIGHT]-1 );
      (S->num[NUM_1TIGHT])--;
      (S->num[NUM_2TIGHT])++;
    }
    else{
      exchange( S->Pi, v->order,
		S->num[NUM_SOL]+S->num[NUM_FREE]+S->num[NUM_1TIGHT]+S->num[NUM_2TIGHT]-1 );
      (S->num[NUM_2TIGHT])--;
      (S->num[NUM_3TIGHT])++;
    }
  }
}

/***** drop a node from the solution *****/
void drop( Solution *S, int p ){
  Dim root_d,d;
  Dir dir;
  Node *v,*w;
#ifdef DEBUG
  if( p < 0 || p >= S->num[NUM_SOL] ){
    fprintf( stderr, "error: you tried to drop a non-solution node\n" );
    exit( EXIT_FAILURE );
  }    
#endif
  /*** exchange ***/
  v = S->Pi[p];
  exchange( S->Pi, p, S->num[NUM_SOL]-1 );
  v->isSol = FALSE;
  v->tightness = 0;
  (S->num[NUM_SOL])--;
  (S->num[NUM_FREE])++;
  /*** process neighbors ***/
  root_d = 0;
  d = 0;
  dir = FORWARD;
  w = v;
  while( 1 ){
    w = getNextNeighbor( v, w, &root_d, &d, &dir );
    if( w == NULL )
      break;
    w->solneigh[d] = NULL;
    (w->tightness)--;
    if( w->tightness == 0 ){
      out1Tight( S, w->order );
      exchange( S->Pi, w->order, S->num[NUM_SOL]+S->num[NUM_FREE] );
      (S->num[NUM_FREE])++;
      (S->num[NUM_1TIGHT])--;
    }
    else if( w->tightness == 1 ){
      in1Tight( S, w->order );
      exchange( S->Pi, w->order, S->num[NUM_SOL]+S->num[NUM_FREE]+S->num[NUM_1TIGHT] );
      (S->num[NUM_1TIGHT])++;
      (S->num[NUM_2TIGHT])--;
    }
    else{
      exchange( S->Pi, w->order, 
		S->num[NUM_SOL]+S->num[NUM_FREE]+S->num[NUM_1TIGHT]+S->num[NUM_2TIGHT] );
      (S->num[NUM_2TIGHT])++;
      (S->num[NUM_3TIGHT])--;
    }
  }  
}


/***** exchange nodes in the ordering *****/
void exchange( Node **Pi, int p, int q ){
  swap( Pi, p, q, sizeof(Node*) );
  Pi[p]->order = p;
  Pi[q]->order = q;
}


/***** process num1t & first1t when S->Pi[p] becomes 1-tight *****/
void in1Tight( Solution *S, int p ){
  Node *v,*first1t;
  Dir dir;
  int t,x,y;
  v = S->Pi[p];
  for(t=0;t<3;t++){
    switch( t ){
    case 0: dir=HOR; x=v->pos[COL]; y=v->pos[VAL]; break;
    case 1: dir=HOR; x=v->pos[ROW]; y=v->pos[VAL]; break;
    case 2: dir=VER; x=v->pos[ROW]; y=v->pos[COL]; break;
    }
    if( t==1 && v->pos[ROW]==v->pos[COL] )
      continue;
    S->num1t[dir][x][y]++;
    first1t = S->first1t[dir][x][y];
    if( first1t == NULL )
      S->first1t[dir][x][y] = v;
    else if( t<=1 && 
	     ( v->pos[ROW] < first1t->pos[ROW] ||
	       v->pos[COL] < first1t->pos[COL] ) )
      S->first1t[dir][x][y] = v;
    else if( t==2 && v->pos[VAL] < first1t->pos[VAL] )
      S->first1t[dir][x][y] = v;
  }
}


/***** process num1t & first1t when S->Pi[p] becomes non-1-tight *****/
void out1Tight( Solution *S, int p ){
  Node *v,*w;
  Dim t,s;
  Dir dir,e;
  int x,y;
  v = S->Pi[p];
  for(t=0;t<3;t++){
    switch( t ){
    case 0: dir=HOR; x=v->pos[COL]; y=v->pos[VAL]; break;
    case 1: dir=HOR; x=v->pos[ROW]; y=v->pos[VAL]; break;
    case 2: dir=VER; x=v->pos[ROW]; y=v->pos[COL]; break;
    }
    if( t==1 && v->pos[ROW]==v->pos[COL] )
      continue;
    S->num1t[dir][x][y]--;
    if( S->num1t[dir][x][y] == 0 )
      S->first1t[dir][x][y] = NULL;
    else if( S->first1t[dir][x][y] == v ){
      S->first1t[dir][x][y] = NULL;
      w = v;
      s = t;
      e = FORWARD;
      while( 1 ){
	w = getNextNeighborDim( v, w, t, &s, &e );
	if( w->tightness == 1 ){
	  S->first1t[dir][x][y] = w;
	  break;
	}
	if( e == BACKWARD ){
	  fprintf( stderr, "error: something unexpected occurred at out1Tight.\n" );
	  exit( EXIT_FAILURE );
	}
      }
    }
  }
}


/***** get the next neighbor of v *****/
Node *getNextNeighbor( Node *root, Node *v, Dim *root_dp, Dim *dp, Dir *dirp ){
  Node *w;
  w = getNextNeighborDim( root, v, *root_dp, dp, dirp );
  if( w != NULL )
    return w;
  while( 1 ){
    (*root_dp)++;
    if( *root_dp == DIM )
      return NULL;
    *dp = *root_dp;
    *dirp = FORWARD;
    w = getNextNeighborDim( root, root, *root_dp, dp, dirp );
    if( w != NULL )
      break;
  }
  return w;
}
 
/***** get the next neighbor of v for given direction *****/
Node *getNextNeighborDim( Node *root, Node *v, Dim root_d, Dim *dp, Dir *dirp ){
  Node *w=NULL;    
  if( *dirp == FORWARD ){
    w = v->next[*dp];
    if( w != NULL ){
      if( root_d == COL && w->pos[ROW] >= root->pos[ROW] && w->pos[COL] == root->pos[ROW] )
	*dp = ROW;
      return w;
    }
    *dp = root_d;
    *dirp = BACKWARD;
    w = root->prev[root_d];
  }
  else
    w = v->prev[*dp];
  if( w != NULL )
    if( root_d == ROW && w->pos[ROW] == root->pos[COL] && w->pos[COL] < root->pos[COL] )
      *dp = COL;
  return w;
}
 
 
int get1TightNeighbor( Node *v, Node **N, int start, Dim d ){
  Node *w;
  Dim e;
  Dir dir = FORWARD;
  int nN = start;
  w = v;
  e = d;
  while( 1 ){
    w = getNextNeighborDim( v, w, d, &e, &dir );
    if( w == NULL )
      break;
    if( w->tightness == 1 ){
      N[nN] = w;
      nN++;
    }
  }
  return nN;
}


Boolean isAdjacent( Node *v, Node *w ){
  if( getDistance( v, w ) <= 1 )
    return TRUE;
  return FALSE;
}

Dim getDistance( Node *v, Node *w ){
  Dim dist=2;
  if( v->pos[ROW] == w->pos[ROW] && v->pos[COL] == w->pos[COL] )
    dist = 0;
  else if( v->pos[ROW] == w->pos[ROW] || v->pos[COL] == w->pos[COL] ||
      v->pos[ROW] == w->pos[COL] || v->pos[COL] == w->pos[ROW] )
    dist = 1;
  if( v->pos[VAL] != w->pos[VAL] )
    dist++;

  return dist;
}


Node *getMostPenalizedNode( Node **N, int numN ){
  Node **F,*v;
  int numF=0,i;
  if( numN == 0 )
    return NULL;
  F = (Node**)malloc(numN*sizeof(Node*));
  F[0] = N[0];
  for(i=1;i<numN;i++)
    if( N[i]->penalty < F[0]->penalty ){
      F[0] = N[i];
      numF = 1;
    }
    else if( N[i]->penalty == F[0]->penalty ){
      F[numF] = N[i];
      numF++;
    }
  v = F[(int)(genrand_real2()*(double)numF)];
  free( F );
  return v;
}


/***** find the 2-tight node that has x & y as its solution neighbors
       and that is not adjacent to u:
       + dx ... direction from u to x
       + dy ... direction from u to y
       + it is assumed dx<dy *****/
Node *getCommon2TightNeighbor( Node ****C, Node *u, Node *x, Node *y, Dim dx, Dim dy ){
  Node *v;
  int px=UNDEF,py=UNDEF;
  if( dy == VAL )
    v = C[ x->pos[ROW] ][ x->pos[COL] ][ y->pos[VAL] ];
  else{
    if( x->pos[ROW] == x->pos[COL] || y->pos[ROW] == y->pos[COL] )
      v = NULL;
    else{
      if( u->pos[ROW] == x->pos[ROW] || u->pos[COL] == x->pos[ROW] )
	px = x->pos[COL];
      else
	px = x->pos[ROW];
      if( u->pos[ROW] == y->pos[ROW] || u->pos[COL] == y->pos[ROW] )
	py = y->pos[COL];
      else
	py = y->pos[ROW];
      if( px < py )
	v = C[ py ][ px ][ x->pos[VAL] ];
      else
	v = C[ px ][ py ][ x->pos[VAL] ];
    }
  }
  if( v == NULL )
    return NULL;
  if( v->tightness != 2 )
    return NULL;
  return v;  
}


Dim getDiffPos( Node *u, Node *v ){
#ifdef DEBUG
  if( u->pos[VAL] != v->pos[VAL] ){
    fprintf( stderr, "error: invalid use of getDiffPos\n" );
    exit( EXIT_FAILURE );
  }
#endif
  if( u->pos[ROW] == u->pos[COL] )
    return UNDEF;
  if( u->pos[ROW] == v->pos[ROW] || u->pos[ROW] == v->pos[COL] )
    return COL;
  if( u->pos[COL] == v->pos[ROW] || u->pos[COL] == v->pos[COL] )
    return ROW;
  fprintf( stderr, "error: u&v are not adjacent in getDiffPos\n" );
  exit( EXIT_FAILURE );    
}


Dim getSamePos( Node *u, Node *v ){
#ifdef DEBUG
  if( u->pos[VAL] != v->pos[VAL] ){
    fprintf( stderr, "error: invalid use of getSamePos\n" );
    exit( EXIT_FAILURE );
  }
#endif
  if( u->pos[ROW] == v->pos[ROW] || u->pos[ROW] == v->pos[COL] )
    return ROW;
  if( u->pos[COL] == v->pos[ROW] || u->pos[COL] == v->pos[COL] )
    return COL;
  fprintf( stderr, "error: u&v are not adjacent in getSamePos\n" );
  exit( EXIT_FAILURE );    
}


Dir getOrthogonalDir( Node *u, Node *v ){
#ifdef DEBUG
  if( u->pos[VAL] != v->pos[VAL] ){
    fprintf( stderr, "error: invalid use of getOrthogonalDir\n" );
    exit( EXIT_FAILURE );
  }
#endif
  if( u->pos[ROW] == u->pos[COL] )
    return UNDEF;
  if( u->pos[ROW] == v->pos[ROW] || u->pos[ROW] == v->pos[COL] )
    return ROW;
  if( u->pos[COL] == v->pos[ROW] || u->pos[COL] == v->pos[COL] )
    return COL;
  fprintf( stderr, "error: u&v are not adjacent in getOrthogonalDir\n" );
  exit( EXIT_FAILURE );    
}


int cmpByIndex( const void *p, const void *q ){
  Node *u,*v;
  u = (Node*)p;
  v = (Node*)q;
  if( u->index < v->index )
    return -1;
  else if( u->index > v->index )
    return 1;
  return 0;
}


void getSolutionDiff( Solution *S, Solution *U,
		      Node **SolS, Node **SolU,
		      Node **S_U, Node **U_S,
		      int *nS_Uptr, int *nU_Sptr ){
  int nSolS,nSolU,i,j,*pos;
  nSolS = S->num[NUM_SOL];
  nSolU = U->num[NUM_SOL];

  /* compute set difference between S and U */
  for(i=0;i<nSolS;i++){
    pos = S->Pi[i]->pos;
    SolS[i] = U->C[pos[ROW]][pos[COL]][pos[VAL]];
  }
  for(i=0;i<nSolU;i++)
    SolU[i] = U->Pi[i];
  qsort( SolS, nSolS, sizeof(Node*), cmpByIndex );
  qsort( SolU, nSolU, sizeof(Node*), cmpByIndex );

  /* elementary algorithm that computes set difference */
  i = 0;
  j = 0;
  *nS_Uptr = 0;
  *nU_Sptr = 0;
  while( i < nSolS || j < nSolU ){
    if( i == nSolS ){
      U_S[(*nU_Sptr)] = SolU[j];
      j++;
      (*nU_Sptr)++;
      continue;
    }
    else if( j == nSolU ){
      S_U[(*nS_Uptr)] = SolS[i];
      i++;
      (*nS_Uptr)++;
      continue;
    }
    if( SolS[i]->index == SolU[j]->index ){
      i++;
      j++;
    }
    else if( SolS[i]->index < SolU[j]->index ){
      S_U[(*nS_Uptr)] = SolS[i];
      i++;
      (*nS_Uptr)++;
    }
    else{
      U_S[(*nU_Sptr)] = SolU[j];
      j++;
      (*nU_Sptr)++;
    }
  }
}

void exchange1TightRandomly( Solution *S ){
  Dim d;
  int i,imax,p;
  Node *u;
  imax = S->num[NUM_1TIGHT];
  for(i=0;i<imax;i++){
    if( genrand_real2() < 0.333333 )
      break;
    p = S->num[NUM_SOL] + S->num[NUM_FREE] 
      + (int)( genrand_real2() * (double)S->num[NUM_1TIGHT] );
    for(d=0;d<DIM;d++)
      if( S->Pi[p]->solneigh[d] != NULL ){
	u = S->Pi[p];
	drop( S, u->solneigh[d]->order );
	insert( S, u->order );
	break;
      }
  }
}


void exchange2TightRandomly( Solution *S ){
  Dim d,dx,dy;
  int i,imax,p;
  Node *u,*v,*x,*y;
  imax = S->num[NUM_2TIGHT];
  for(i=0;i<imax;i++){
    if( genrand_real2() < 0.333333 )
      break;
    p = S->num[NUM_SOL] + S->num[NUM_FREE] + S->num[NUM_1TIGHT] 
      + (int)( genrand_real2() * (double)S->num[NUM_2TIGHT] );
    u = S->Pi[p];
    dx = UNDEF;
    dy = UNDEF;
    x = NULL;
    y = NULL;
    for(d=0;d<DIM;d++)
      if( u->solneigh[d] == NULL || ( d == 1 && u->pos[ROW] == u->pos[COL] ) )
	continue;
      else if( x == NULL ){
	x = u->solneigh[d];
	dx = d;
      }
      else if( y == NULL ){
	y = u->solneigh[d];
	dy = d;
      }
    v = getCommon2TightNeighbor( S->C, u, x, y, dx, dy );
    if( v != NULL ){
      drop( S, x->order );
      drop( S, y->order );
      insert( S, u->order );
      insert( S, v->order );
    }
  }
}
