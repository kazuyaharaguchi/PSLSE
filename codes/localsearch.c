/******************************
  localsearch.c
******************************/

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "elementary.h"
#include "misc.h"
#include "mt19937ar.h"
#include "localsearch.h"

extern double cpu_time( void );

/***** iterated local search *****/
void ILS( ILSParam *IP, ILSRec *IR, Solution *S ){
  IR->timeIlsStart = cpu_time();
  init_genrand( IP->seed );
  switch( IP->alg ){
  case ILS1:
    execILS1( IP, IR, S );
    break;
  case ILS2:
    execILS2( IP, IR, S );
    break;
  default:
    fprintf( stderr, "error: -alg=%d is inappropriate.\n", IP->alg );
    exit( EXIT_FAILURE );
  }
  IR->timeIlsEnd = cpu_time();
}


/***** 1-ILS *****/
void execILS1( ILSParam *IP, ILSRec *IR, Solution *S ){
  while( 1 ){
    LS1( IP, IR, S );
    if( updateBest( IP, IR, S ) )
      break;
    perturb( IP, IR, S );
  } 
}


/***** 2-ILS *****/
void execILS2( ILSParam *IP, ILSRec *IR, Solution *S ){
  while( 1 ){
    if( LS1( IP, IR, S ) )
      continue;
    if( LS2( IP, IR, S ) )
      continue;
    if( updateBest( IP, IR, S ) )
      break;
    perturb( IP, IR, S );
  } 
}


/***** update the best solution &
       decide whether the termination condition is satisfied *****/
Boolean updateBest( ILSParam *IP, ILSRec *IR, Solution *S ){
  Node **SolS,**SolBest,**S_B,**B_S,*v;
  Dim d;
  int nSolS,nSolBest,nS_B=0,nB_S=0,i,j,pos[DIM];

  /*** preprocess ***/
  nSolS = S->num[NUM_SOL];
  nSolBest = IR->Best->num[NUM_SOL];
  SolS    = (Node**)malloc(nSolS*sizeof(Node*));
  SolBest = (Node**)malloc(nSolBest*sizeof(Node*));
  S_B     = (Node**)malloc(nSolS*sizeof(Node*));
  B_S     = (Node**)malloc(nSolBest*sizeof(Node*));

  /*** compute set difference ***/
  getSolutionDiff( S, IR->Best, SolS, SolBest, S_B, B_S, &nS_B, &nB_S );
  shuffle( S_B, nS_B, sizeof(Node*), genrand_real2 );
  shuffle( B_S, nB_S, sizeof(Node*), genrand_real2 );

  /*** update the incumbent solution ***/
  if( nSolS > nSolBest || ( nSolS == nSolBest && ( nS_B || nB_S ) ) ){
    for(j=0;j<nB_S;j++)
      drop( IR->Best, B_S[j]->order );  
    for(i=0;i<nS_B;i++)
      insert( IR->Best, S_B[i]->order );

    // best has been found
    if( nSolS > nSolBest ){
      IR->lastImproved = IR->t;
      printf("Improved:\t%d\t%d\t%g\n",S->input_ass+S->num[NUM_SOL],
	     IR->t+1,cpu_time()-IR->timeIlsStart);
    }
  }
  else{
    exchange1TightRandomly( IR->Best );
    exchange2TightRandomly( IR->Best );
    //printf("%d,%d ",IR->Best->num[NUM_1TIGHT],IR->Best->num[NUM_2TIGHT]); fflush(stdout);
    getSolutionDiff( S, IR->Best, SolS, SolBest, S_B, B_S, &nS_B, &nB_S );
    shuffle( S_B, nS_B, sizeof(Node*), genrand_real2 );
    shuffle( B_S, nB_S, sizeof(Node*), genrand_real2 );
    for(i=0;i<nS_B;i++){
      for(d=0;d<DIM;d++)
	pos[d] = S_B[i]->pos[d];
      v = S->C[pos[ROW]][pos[COL]][pos[VAL]];
      drop( S, v->order );
    }
    for(j=0;j<nB_S;j++){
      for(d=0;d<DIM;d++)
	pos[d] = B_S[j]->pos[d];
      v = S->C[pos[ROW]][pos[COL]][pos[VAL]];
      insert( S, v->order );
    }
  }

  free( SolS );
  free( SolBest );
  free( S_B );
  free( B_S );

  /*** increase the counter ***/
  (IR->t)++;
  
  /*** decide whether the termination condition is satisfied ***/
  /* upper bound */
  if( IR->Best->num[NUM_SOL] == S->UB )
    return TRUE;
  /* computation time */
  if( IP->timemax >= 0.0 && IP->timemax < cpu_time()-IR->timeIlsStart )
    return TRUE;
  /* iteration time */
  if( IP->itrmax >= 0 && IP->itrmax < IR->t )
    return TRUE;
  return FALSE;
}


void perturb( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean firstFlag=TRUE;
  double breakprob=IP->breakprob;
  Node **N,*u,*v,*x;
  Dim d,d_,root_d;
  Dir dir;
  int i,nS0,n123,numN;
 
  /*** kick ***/
  while( 1 ){
    nS0 = S->num[NUM_SOL] + S->num[NUM_FREE];
    n123 = S->Nodes - nS0;
    if( firstFlag ){
      firstFlag = FALSE;
      if( S->num[NUM_1TIGHT] == 0 )
	u = getMostPenalizedNode( S->Pi+nS0, n123 );
      else{
	N = (Node**)malloc(3*S->n*S->num[NUM_1TIGHT]*sizeof(Node*));
	numN = 0;
	for(i=0;i<S->num[NUM_1TIGHT];i++){
	  u = S->Pi[nS0+i];
	  for(d=0;d<DIM;d++)
	    if( u->solneigh[d] != NULL )
	      break;
	  x = u->solneigh[d];
	  root_d = ROW;
	  d_ = root_d;
	  dir = FORWARD;
	  v = x;
	  while( 1 ){
	    v = getNextNeighbor( x, v, &root_d, &d_, &dir );
	    if( v == NULL )
	      break;
	    N[numN] = v;
	    numN++;
	  }
	}
	u = getMostPenalizedNode( N, numN );
	free( N );
      }
    }
    else
      u = S->Pi[ nS0 + (int)( genrand_real2() * (double)n123 ) ];
    
    /* drop solution neighbors of u */
    for(d=0;d<DIM;d++)
      if( u->solneigh[d] != NULL )
	drop( S, u->solneigh[d]->order );
    /* insert u */
    insert( S, u->order );
    /* if free nodes appear, insert them */
    while( S->num[NUM_FREE] > 0 )
      insert( S, S->Pi[ S->num[NUM_SOL] + (int)( genrand_real2() * (double)S->num[NUM_FREE]) ]->order );
    /* break */
    if( genrand_real2() < breakprob || S->num[NUM_SOL] == S->UB )
      break;
  }
  /*** update the penalties ***/
  for(i=0;i<S->num[NUM_SOL];i++)
    S->Pi[i]->penalty = IR->t;
}



void lineSwitch( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean *cov;
  Node *u;
  int i,j,h,k,ub;
  int *Uncov,numUncov=0,row,col,L[2],*V[2],l,tmp;

  printf("LINE_SWITCH!\t%d\t(%d) ---> ",S->num[NUM_SOL]+S->input_ass,IR->Best->num[NUM_SOL]+S->input_ass);

  // preprocess
  ub = S->n * (S->n+1) / 2;
  cov = (Boolean*)malloc(ub*sizeof(Boolean));
  for(h=0;h<ub;h++)
    cov[h] = FALSE;
  for(i=0;i<S->num[NUM_SOL];i++){
    h = S->Pi[i]->pos[ROW] * (S->Pi[i]->pos[ROW]+1) / 2 + S->Pi[i]->pos[COL];
    cov[h] = TRUE;
  }
  Uncov=(int*)malloc(ub*sizeof(int));
  for(h=0;h<S->num[NUM_SOL];h++)
    if( !cov[h] ){
      Uncov[numUncov] = h;
      numUncov++;
    }
  h = Uncov[ (int)(genrand_real2() * (double)numUncov) ];
  row = S->n-1;
  for(i=0;i<S->n;i++)
    if( i*(i+1)/2 > h ){
      row = i-1;
      break;
    }
  col = h - row*(row+1)/2;
  if( genrand_real2() < 0.5 )
    L[0] = row;
  else
    L[0] = col;
  while( 1 ){
    L[1] = (int)( genrand_real2() * (double)S->n );
    if( L[1] != L[0] )
      break;
  }
  if( L[0] > L[1] ){
    tmp = L[0];
    L[0] = L[1];
    L[1] = tmp;
  }

#if 0
  {/********************/
    int **Latin;
    Latin = (int**)malloc(S->n*sizeof(int*));
    for(i=0;i<S->n;i++){
      Latin[i] = (int*)malloc(S->n*sizeof(int));
      for(j=0;j<S->n;j++)
	Latin[i][j] = UNDEF;
    }
    for(k=0;k<S->num[NUM_SOL];k++){
      i = S->Pi[k]->pos[ROW];
      j = S->Pi[k]->pos[COL];
      Latin[i][j] = S->Pi[k]->pos[VAL];
    }
    outputLatinSquare( stdout, Latin, S->n );
    free( Latin );
    /********************/}
#endif

  for(l=0;l<2;l++){
    V[l] = (int*)malloc(S->n*sizeof(int));
    for(j=0;j<S->n;j++){
      V[l][j] = UNDEF;
      for(k=0;k<S->n;k++){
	if( j < L[l] )
	  u = S->C[L[l]][j][k];
	else
	  u = S->C[j][L[l]][k];
	if( u != NULL && u->isSol ){
	  V[l][j] = k;
	  break;
	}
      }
    }	    
  }
  for(l=0;l<2;l++)
    for(j=0;j<S->n;j++)
      for(k=0;k<S->n;k++){
	if( j < L[l] )
	  u = S->C[L[l]][j][k];
	else
	  u = S->C[j][L[l]][k];
	if( u != NULL && u->isSol ){
	  drop( S, u->order );
	  break;
	}
      }


#if 0
  {/********************/
    printf("------------------------------\n");
    int **Latin;
    Latin = (int**)malloc(S->n*sizeof(int*));
    for(i=0;i<S->n;i++){
      Latin[i] = (int*)malloc(S->n*sizeof(int));
      for(j=0;j<S->n;j++)
	Latin[i][j] = UNDEF;
    }
    for(k=0;k<S->num[NUM_SOL];k++){
      i = S->Pi[k]->pos[ROW];
      j = S->Pi[k]->pos[COL];
      Latin[i][j] = S->Pi[k]->pos[VAL];
    }
    outputLatinSquare( stdout, Latin, S->n );
    free( Latin );
    /********************/}
#endif

  for(j=0;j<S->n;j++)
    if( j == L[0] ){
      if( V[1][j] != UNDEF ){
	u = S->C[L[1]][j][V[1][j]];
	if( u != NULL )
	  insert( S, u->order );
      }
      if( V[1][L[1]] != UNDEF ){
	u = S->C[L[0]][L[0]][V[1][L[1]]];
	if( u != NULL )
	  insert( S, u->order );
      }
    }
    else if( j == L[1] ){
      if( V[0][L[0]] != UNDEF ){
	u = S->C[L[1]][L[1]][V[0][L[0]]];
	if( u != NULL )
	  insert( S, u->order );
      }
    }
    else{
      for(l=0;l<2;l++){
	if( V[1-l][j] == UNDEF )
	  continue;
	if( j < L[l] )
	  u = S->C[L[l]][j][V[1-l][j]];
	else if( j > L[l] )
	  u = S->C[j][L[l]][V[1-l][j]];
	if( u != NULL )
	  insert( S, u->order );
      }
    }

#if 0
  {/********************/
    printf("------------------------------\n");
    int **Latin;
    Latin = (int**)malloc(S->n*sizeof(int*));
    for(i=0;i<S->n;i++){
      Latin[i] = (int*)malloc(S->n*sizeof(int));
      for(j=0;j<S->n;j++)
	Latin[i][j] = UNDEF;
    }
    for(k=0;k<S->num[NUM_SOL];k++){
      i = S->Pi[k]->pos[ROW];
      j = S->Pi[k]->pos[COL];
      Latin[i][j] = S->Pi[k]->pos[VAL];
    }
    outputLatinSquare( stdout, Latin, S->n );
    free( Latin );
    /********************/}
#endif

  printf("%d\n",S->num[NUM_SOL]);
	    
  free( cov );
  free( Uncov );
  free( V[0] );
  free( V[1] );
}


void symbolSwitch( ILSParam *IP, ILSRec *IR, Solution *S ){
  Node **V,*u;
  int l,j,k,numV=0;
  
  printf("SYMB_SWITCH!\t%d\t(%d) ---> ",S->num[NUM_SOL]+S->input_ass,IR->Best->num[NUM_SOL]+S->input_ass);

  l = (int)( genrand_real2() * S->n );
  V = (Node**)malloc(S->n*sizeof(Node*));
  for(j=0;j<S->n;j++){
    if( j == l )
      continue;
    for(k=0;k<S->n;k++){
      if( j < l )
	u = S->C[l][j][k];
      else
	u = S->C[j][l][k];
      if( u != NULL && u->isSol ){
	V[numV] = u;
	numV++;
	break;
      }
    }
  }
  if( numV <= 1 ){
    printf("%d\n",S->num[NUM_SOL]);
    free( V );
    return;
  }
  shuffle( V, numV, sizeof(Node*), genrand_real2 );
  symbolSwitchCore( IP, IR, S, l, V );
  free( V );

  printf("%d\n",S->num[NUM_SOL]);

}


void symbolSwitchCore( ILSParam *IP, ILSRec *IR, Solution *S, int l, Node **V ){
  Boolean flag;
  Node **A,*u,*v;
  int numA=0,p,q,r,t,h,k;
#if 0
    {/********************/
      int **Latin,i,j;
    Latin = (int**)malloc(S->n*sizeof(int*));
    for(i=0;i<S->n;i++){
      Latin[i] = (int*)malloc(S->n*sizeof(int));
      for(j=0;j<S->n;j++)
	Latin[i][j] = UNDEF;
    }
    for(k=0;k<S->num[NUM_SOL];k++){
      i = S->Pi[k]->pos[ROW];
      j = S->Pi[k]->pos[COL];
      Latin[i][j] = S->Pi[k]->pos[VAL];
    }
    outputLatinSquare( stdout, Latin, S->n );
    free( Latin );
    printf("--------------------\n");
    /********************/}
#endif

  A = (Node**)malloc(2*S->n*sizeof(Node*));
  for(t=0;t<2;t++){
    u = V[t];
    p = 1;
    q = l;
    while( 1 ){
      flag = FALSE;
      h = V[(t+p)%2]->pos[VAL];
      if( u->pos[ROW] == q )
	r = u->pos[COL];
      else
	r = u->pos[ROW];
      for(k=0;k<S->n;k++){
	if( k < r )
	  v = S->C[r][k][h];
	else
	  v = S->C[k][r][h];
	if( v != NULL && v->isSol ){
	  if( v->pos[ROW] != v->pos[COL] && k != l ){
	    A[numA] = v;
	    numA++;
	    flag = TRUE;
	  }
	  break;
	}
      }
      if( !flag )
	break;
      u = v;
      p++;
      q = r;
    }
    if( k == l )
      break;
  }

  /*** drop ***/
  for(t=0;t<2;t++)
    drop( S, V[t]->order );
  for(k=0;k<numA;k++)
    drop( S, A[k]->order );

  /*** insert ***/
  for(t=0;t<2;t++){
    v = S->C[V[t]->pos[ROW]][V[t]->pos[COL]][V[1-t]->pos[VAL]];
    if( v != NULL )
      insert( S, v->order );
  }
  for(k=0;k<numA;k++){
    if( A[k]->pos[VAL] == V[0]->pos[VAL] )
      v = S->C[A[k]->pos[ROW]][A[k]->pos[COL]][V[1]->pos[VAL]];
    else
      v = S->C[A[k]->pos[ROW]][A[k]->pos[COL]][V[0]->pos[VAL]];
    if( v != NULL )
      insert( S, v->order );
  }
  free( A );

#if 0
    {/********************/
      int **Latin,i,j;
    Latin = (int**)malloc(S->n*sizeof(int*));
    for(i=0;i<S->n;i++){
      Latin[i] = (int*)malloc(S->n*sizeof(int));
      for(j=0;j<S->n;j++)
	Latin[i][j] = UNDEF;
    }
    for(k=0;k<S->num[NUM_SOL];k++){
      i = S->Pi[k]->pos[ROW];
      j = S->Pi[k]->pos[COL];
      Latin[i][j] = S->Pi[k]->pos[VAL];
    }
    outputLatinSquare( stdout, Latin, S->n );
    free( Latin );
    printf("--------------------\n");
    /********************/}
  exit( 1 );
#endif
  
}


/***** find (1,n^2)-maximal solution *****/
Boolean LS1( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean flagOverall=FALSE,flagEachTime;
  Node *x=NULL,*u=NULL;
  Dim d;
  int p,head,tail;
  /*** if there is at most one 1-tight node,
       then S is already 1-maximal ***/
  if( S->num[NUM_1TIGHT] <= 1 )
    return FALSE;
  /*** iteration ***/
  while( 1 ){
    flagEachTime = FALSE;
    head = S->num[NUM_SOL] + S->num[NUM_FREE];
    tail = head + S->num[NUM_1TIGHT];
    for(p=head;p<tail;p++){
      /* detect the only 1-tight neighbor of u */
      u = S->Pi[(IR->randptr++)%(tail-head)+head];
      for(d=0;d<DIM;d++)
	if( u->solneigh[d] != NULL ){
	  x = u->solneigh[d];
	  break;
	}
      if( LS1Core( IP, IR, S, x ) ){
	flagOverall = TRUE;
	flagEachTime = TRUE;
	break;
      }
    }
    if( !flagEachTime )
      break;
  }
  return flagOverall;
}


Boolean LS1Core( ILSParam *IP, ILSRec *IR, Solution *S, Node *x ){
  Node *N[DIM],*v;
  Dir dir;
  int nMIS=0,a,b,i,j,t,tmp,numN=0;

  /*** count the MIS size among free nodes from S-x ***/
  if( x->pos[ROW] == x->pos[COL] ){
    if( S->num1t[HOR][ x->pos[ROW] ][ x->pos[VAL] ] > 0 )
      nMIS++;
  }
  else{
    a = S->num1t[HOR][ x->pos[ROW] ][ x->pos[VAL] ];
    b = S->num1t[HOR][ x->pos[COL] ][ x->pos[VAL] ];
    if( a<b ){ tmp=a; a=b; b=tmp; }
    if( a>1 ){
      nMIS++;
      if( b>0 )
	nMIS++;
    }
    else if( a==1 ){
      nMIS++;
      if( b==1 &&
	  !isAdjacent( S->first1t[HOR][x->pos[ROW]][x->pos[VAL]],
		       S->first1t[HOR][x->pos[COL]][x->pos[VAL]] ) )
	nMIS++;
    }
  }
  if( S->num1t[VER][ x->pos[ROW] ][ x->pos[COL] ] > 0 )
    nMIS++;
  if( nMIS <= 1 )
    return FALSE;

  /**********/
  /*
    printf("\n--------------------\n");
    printf("[%d,%d,%d]\tnMIS=%d\n",x->pos[0],x->pos[1],x->pos[2],nMIS);
    for(t=0;t<3;t++){
    switch( t ){
    case 0: dir=HOR; i=x->pos[COL]; j=x->pos[VAL]; break;
    case 1: dir=HOR; i=x->pos[ROW]; j=x->pos[VAL]; break;
    case 2: dir=VER; i=x->pos[ROW]; j=x->pos[COL]; break;
    }
    printf("\t(t=%d) num1t=%d",t,S->num1t[dir][i][j]);
    if( S->num1t[dir][i][j]>0 ){
    v = S->first1t[dir][i][j];
    printf(" [%d,%d,%d]",v->pos[0],v->pos[1],v->pos[2]);
    }
    printf("\n");
    }
  */
  /**********/

  /*** collect the only 1-tight node from all hyperedges ***/
  for(t=0;t<3;t++){
    switch( t ){
    case 0: dir=HOR; i=x->pos[COL]; j=x->pos[VAL]; break;
    case 1: dir=HOR; i=x->pos[ROW]; j=x->pos[VAL]; break;
    case 2: dir=VER; i=x->pos[ROW]; j=x->pos[COL]; break;
    }
    if( S->num1t[dir][i][j] != 1 )
      continue;
    if( t==1 && S->num1t[HOR][x->pos[COL]][x->pos[VAL]] == 1 &&
	S->num1t[HOR][x->pos[ROW]][x->pos[VAL]] == 1 &&
	isAdjacent( S->first1t[HOR][x->pos[COL]][x->pos[VAL]],
		    S->first1t[HOR][x->pos[ROW]][x->pos[VAL]] ) ) 
      continue;
    
    /********************/
    /*
      v = S->first1t[dir][i][j];
      printf(" [%d,%d,%d] is to be inserted\n",v->pos[0],v->pos[1],v->pos[2]);
    */
    /********************/

    N[numN] = S->first1t[dir][i][j];
    numN++;
  }
  drop( S, x->order );
  for(t=0;t<numN;t++)
    insert( S, N[t]->order );
  while( S->num[NUM_FREE] > 0 ){

    /********************/
    /*
      int i;
      for(i=0;i<S->num[NUM_FREE];i++){
      v = S->Pi[S->num[NUM_SOL]+i];
      printf(" [%d,%d,%d]",v->pos[0],v->pos[1],v->pos[2]);
      }
    */
    /********************/

    if( S->Pi[S->num[NUM_SOL]] == x )
      v = getMostPenalizedNode( &(S->Pi[S->num[NUM_SOL]+1]),
				S->num[NUM_FREE]-1 );
    else
      v = getMostPenalizedNode( &(S->Pi[S->num[NUM_SOL]]),
				S->num[NUM_FREE] );

    /********************/
    //printf(" ---> [%d,%d,%d] is inserted\n",v->pos[0],v->pos[1],v->pos[2]);
    /********************/

    insert( S, v->order );
    numN++;
  }

#ifdef DEBUG
  if( numN != nMIS ){
    fprintf( stderr, "error: theoretical MIS size (%d) != number of nodes that are added (%d)\n",
	     nMIS, numN );
    exit( EXIT_FAILURE );
  }
#endif

  return TRUE;
}


/***** find (2,n^2)-maximal solution *****/
Boolean LS2( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean flag=FALSE;
  Node *u;
  int p,head,tail;
  head = S->num[NUM_SOL] + S->num[NUM_FREE] + S->num[NUM_1TIGHT];
  tail = head + S->num[NUM_2TIGHT];
  for(p=head;p<tail;p++){
    u = S->Pi[(IR->randptr++)%(tail-head)+head];
    if( u->tightness == 2 ){
      if( LS2Core( IP, IR, S, u ) ){
	flag = TRUE;
	break;
      }
    }
  }
  return flag;
}


Boolean LS2Core( ILSParam *IP, ILSRec *IR, Solution *S, Node *u ){
  Boolean is2TightPossible=FALSE,isConflict=FALSE;
  Node **U,**O,*x=NULL,*y=NULL,*v=NULL,*w=NULL;
  Dim d,dx=UNDEF,dy=UNDEF;
  Dir dir;
  int ox,oy,px,py,mx,my,oz,ow,pz,pw,mz,mw,
    nMIS=0,numU=0,numO=0,k,l,a,b;
#ifdef DEBUG
  int prevsize;
#endif

  /*** detect solution neighbors of u and the directions from u ***/
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

  /********************/
  /*
  printf("\n\tI am [%d,%d,%d]. dx=%d, dy=%d. ",
	 u->pos[0],u->pos[1],u->pos[2],dx,dy);
  printf("Solution neighbors are: [%d,%d,%d] (dx=%d) and [%d,%d,%d] (dy=%d)\n",
	 x->pos[0],x->pos[1],x->pos[2],dx,y->pos[0],y->pos[1],y->pos[2],dy);
  */
  /********************/


  /*** count the number of nodes that can be added ***/
  /* 1-tight nodes */
  if( dx == ROW && dy == COL ){
    // those on the same row-col plane
    ox = getDiffPos( x, u );
    if( ox == UNDEF ){
      px = UNDEF;
      mx = 0;
    }
    else{
      px = x->pos[ox];
      mx = S->num1t[HOR][px][u->pos[VAL]];
    }
    oy = getDiffPos( y, u );
    if( oy == UNDEF ){
      py = UNDEF;
      my = 0;
    }
    else{
      py = y->pos[oy];
      my = S->num1t[HOR][py][u->pos[VAL]];
    }
    if( mx >= 1 && my >= 1 ){
      if( mx == 1 && my == 1 && 
	  isAdjacent( S->first1t[HOR][px][u->pos[VAL]],
		      S->first1t[HOR][py][u->pos[VAL]] ) )
	nMIS = 1;
      else
	nMIS = 2;
    }
    else if( mx+my >= 1 )
      nMIS = 1;
    else
      is2TightPossible = TRUE;
    // those on perpendicular grids
    if( S->num1t[VER][x->pos[ROW]][x->pos[COL]] ) nMIS++;
    if( S->num1t[VER][y->pos[ROW]][y->pos[COL]] ) nMIS++;


    /**********/
    /*
      printf("\tpx=%d,py=%d\n",px,py);
      if( px == UNDEF )
      printf("\t+num1t[HOR][%d][%d] = 0\n", -1,u->pos[VAL]);
      else
      printf("\t+num1t[HOR][%d][%d] = %d\n", px,u->pos[VAL],S->num1t[HOR][px][u->pos[VAL]]);
      if( py == UNDEF )
      printf("\t+num1t[HOR][%d][%d] = 0\n", -1,u->pos[VAL]);
      else
      printf("\t+num1t[HOR][%d][%d] = %d\n", py,u->pos[VAL],S->num1t[HOR][py][u->pos[VAL]]);
      printf("\t+num1t[VER][%d][%d] = %d\n", x->pos[ROW],x->pos[COL],
      S->num1t[VER][x->pos[ROW]][x->pos[COL]]);
      printf("\t+num1t[VER][%d][%d] = %d\n", y->pos[ROW],y->pos[COL],
      S->num1t[VER][y->pos[ROW]][y->pos[COL]]);
    */
    /**********/

  }
  else{
    is2TightPossible = TRUE;
    ox = getDiffPos( x, u );
    if( ox == UNDEF )
      px = UNDEF;
    else{
      px = x->pos[ox];
      if( S->num1t[HOR][px][x->pos[VAL]] > 0 )
	if( S->num1t[HOR][px][x->pos[VAL]] > 1 ||
	    !isAdjacent( u, S->first1t[HOR][px][x->pos[VAL]] ) ){
	  nMIS++;
	  /*
	  Node *g = S->first1t[HOR][px][x->pos[VAL]];
	  printf("g:[%d,%d,%d]\n",g->pos[0],g->pos[1],g->pos[2]);
	  */
	}
    }
    if( S->num1t[VER][x->pos[ROW]][x->pos[COL]] ){
      is2TightPossible = FALSE;
      nMIS++;
      /*
      Node *g = S->first1t[VER][x->pos[ROW]][x->pos[COL]];
      printf("g:[%d,%d,%d]\n",g->pos[0],g->pos[1],g->pos[2]);
      */
    }

    /**/
    //printf("<nMIS=%d>\n",nMIS);
    //printf("<S->num1t[VER][%d][%d]=%d>\n",x->pos[ROW],x->pos[COL],S->num1t[VER][x->pos[ROW]][x->pos[COL]]);
    /**/


    if( y->pos[ROW]==y->pos[COL] ){
      oz=ROW; ow=ROW;
    }
    else if( y->pos[ROW]==x->pos[ROW] || y->pos[ROW]==x->pos[COL] ){
      oz=ROW; ow=COL;
    }
    else{
      oz=COL; ow=ROW;
    }
    pz = y->pos[oz];
    pw = y->pos[ow];
    mz = S->num1t[HOR][pz][y->pos[VAL]];
    mw = S->num1t[HOR][pw][y->pos[VAL]];
    
    //printf("oz,ow,pz,pw,mz,mw=%d,%d,%d,%d,%d,%d\n",oz,ow,pz,pw,mz,mw);
    
    // decide whether 2-tight node should be taken into account
    w = S->first1t[HOR][pw][y->pos[VAL]];

    /*
    if( w != NULL )
      printf("w:[%d,%d,%d]\n",w->pos[0],w->pos[1],w->pos[2]);
    */

    if( mw == 1 && ( w->pos[ROW] == x->pos[ROW] ||
		     w->pos[ROW] == x->pos[COL] ||
		     w->pos[COL] == x->pos[ROW] ||
		     w->pos[COL] == x->pos[COL] ) )
      is2TightPossible = FALSE;
    else{
      if( mz>1 )
	is2TightPossible = FALSE;
      else if( mz==1 ){
	if( mw==1 ){
	  if( !isAdjacent( S->first1t[HOR][pz][y->pos[VAL]],
			   S->first1t[HOR][pw][y->pos[VAL]] ) )
	    is2TightPossible = FALSE;
	}
	else
	  is2TightPossible = FALSE;
      }
    }

    // decide whether 1-tight nodes in lines pz & pw cannot be inserted simultaneously
    if( y->pos[ROW]==y->pos[COL] )
      isConflict = TRUE;
    else if( mz == 1 && mw == 1 &&
	     isAdjacent( S->first1t[HOR][pz][y->pos[VAL]], S->first1t[HOR][pw][y->pos[VAL]] ) )
      isConflict = TRUE;
    
    // count nMIS
    if( mz ){
      if( mw ){
	if( isConflict )
	  nMIS++;
       	else
	  nMIS+=2;
      }
      else
	nMIS++;
    }
    else if( mw )
      nMIS++;

    /**********/
    /*
      printf("\tox=%d, px=%d\n",ox,px);
      if( px == UNDEF )
      printf("\t+num1t[HOR][-1][%d] = 0\n", x->pos[VAL]);
      else{
      printf("\t+num1t[HOR][%d][%d] = %d\n", px,x->pos[VAL],S->num1t[HOR][px][x->pos[VAL]]);
      Node *z = S->first1t[HOR][px][x->pos[VAL]];
      if( z != NULL ){
      if( isAdjacent( z, u ) )
      printf("\t\t---> [%d,%d,%d] is adjacent to u\n",z->pos[0],z->pos[1],z->pos[2]);
      else
      printf("\t\t---> [%d,%d,%d] is OK\n",z->pos[0],z->pos[1],z->pos[2]);
      }
      }
      printf("\t+num1t[VER][%d][%d] = %d\n", x->pos[ROW],x->pos[COL],
      S->num1t[VER][x->pos[ROW]][x->pos[COL]]);
      printf("\t+num1t[HOR][%d][%d] = %d\n", y->pos[ROW],y->pos[VAL],
      S->num1t[HOR][y->pos[ROW]][y->pos[VAL]]);
      printf("\t+num1t[HOR][%d][%d] = %d\n", y->pos[COL],y->pos[VAL],
      S->num1t[HOR][y->pos[COL]][y->pos[VAL]]);
    */    
    /**********/
  }

  /* 2-tight node: only one is possible; 
     let u=({row,col},val), and suppose x=({row,a},val) and y=({col,b},val)
     then ({row,b},val),({col,a},val),({a,b},val) are candidates,
     but among these, the first two are adjacent to u. 
  */
  if( is2TightPossible ){
    v = getCommon2TightNeighbor( S->C, u, x, y, dx, dy );
    if( v != NULL )
      nMIS++;
  }

  /**********/
  /*
    if( v != NULL )
    printf("\tv=[%d,%d,%d] (tightness=%d)\n",v->pos[0],v->pos[1],v->pos[2],
    v->tightness);
    printf("\tnMIS=%d\n",nMIS);
    fflush( stdout );
    */
  /**********/

  /*** if no improvement is possible, return FALSE ***/
  if( nMIS <= 1 )
    return FALSE;

  /*** collect the 1-tight nodes that should be inserted ***/
  U = (Node**)malloc(4*sizeof(Node*));
  O = (Node**)malloc(4*S->n*sizeof(Node*));
  if( dx == ROW && dy == COL )
    for(l=0;l<4;l++){
      switch( l ){
      case 0: dir=HOR; a=px; b=u->pos[VAL]; w=x; d=getOrthogonalDir(x,u); break; 
      case 1: dir=HOR; a=py; b=u->pos[VAL]; w=y; d=getOrthogonalDir(y,u); break; 
      case 2: dir=VER; a=x->pos[ROW]; b=x->pos[COL]; w=x; d=VAL; break; 
      case 3: dir=VER; a=y->pos[ROW]; b=y->pos[COL]; w=y; d=VAL; break; 
      }
      if( a == UNDEF )
	continue;
      if( S->num1t[dir][a][b] > 1 )
	numO = get1TightNeighbor( w, O, numO, d );
      else if( S->num1t[dir][a][b] == 1 ){
	U[numU] = S->first1t[dir][a][b];
	numU++;
      }
    }  
  else
    for(l=0;l<4;l++){
      switch( l ){
      case 0:
	dir = HOR;
	a   = px;
	b   = x->pos[VAL];
	w   = x;
	d   = getOrthogonalDir(x,u);
	break; 
      case 1:
	dir = VER;
	a   = x->pos[ROW];
	b   = x->pos[COL];
	w   = x;
	d   = VAL;
	break; 
      case 2: 
	dir = HOR;
	a   = y->pos[ROW];
	b   = y->pos[VAL];
	w   = y;
	d   = COL;
	break; 
      case 3:
	dir = HOR;
	a   = y->pos[COL];
	b   = y->pos[VAL];
	w   = y;
	d   = ROW;
	break; 
      }
      if( a == UNDEF )
	continue;
      //if( l>=2 && isConflict && d==ow )
      //continue;
      if( S->num1t[dir][a][b] > 1 )
	numO = get1TightNeighbor( w, O, numO, d );
      else if( S->num1t[dir][a][b] == 1 ){
	U[numU] = S->first1t[dir][a][b];
	numU++;
      }
    }
  if( numO )
    shuffle( O, numO, sizeof(Node*), genrand_real2 );
  if( numU )
    shuffle( U, numU, sizeof(Node*), genrand_real2 );

  /**********/
  /*
  printf("---------- improvement by LS2 (nMIS=%d) ----------\n",nMIS);
  printf("\t-ow=%d, isConflict=%d\n",ow,isConflict);
    printf("\t-x=[%d,%d,%d]\n",x->pos[0],x->pos[1],x->pos[2]);
    printf("\t-y=[%d,%d,%d]\n",y->pos[0],y->pos[1],y->pos[2]);
    printf("\t+u=[%d,%d,%d]\n",u->pos[0],u->pos[1],u->pos[2]);
    if( v != NULL )
    printf("\t+v=[%d,%d,%d]\n",v->pos[0],v->pos[1],v->pos[2]);

    printf("\t+U:");
    for(k=0;k<numU;k++)
    printf(" [%d,%d,%d](%d)",U[k]->pos[0],U[k]->pos[1],U[k]->pos[2],U[k]->tightness);
    printf("\n");
    printf("\t+O:");
    for(k=0;k<numO;k++)
    printf(" [%d,%d,%d](%d)",O[k]->pos[0],O[k]->pos[1],O[k]->pos[2],O[k]->tightness);
    printf("\n");
    fflush( stdout );
    */
  /**********/

  /*** drop x & y ***/
  drop( S, x->order );
  drop( S, y->order );

  /*** add u ***/
  insert( S, u->order );

#ifdef DEBUG
  prevsize = S->num[NUM_SOL];
#endif

  /*** add 2-tight node ***/
  if( v != NULL )
    insert( S, v->order );

  /** add nodes in U ***/
  for(k=0;k<numU;k++)
    if( U[k]->isSol == FALSE && U[k]->tightness == 0 )
      insert( S, U[k]->order );

  /** add nodes in O ***/
  for(k=0;k<numO;k++)
    if( O[k]->isSol == FALSE && O[k]->tightness == 0 )
      insert( S, O[k]->order );
  
  /********************/
  /*
    printf("\n(Free:");
    for(k=0;k<S->num[NUM_FREE];k++){
    v = S->Pi[k+S->num[NUM_SOL]];
    printf(" [%d,%d,%d](%d)",v->pos[0],v->pos[1],v->pos[2],v->tightness);
    }
    printf(")\n");
    if( S->num[NUM_FREE] ){
    fflush( stdout );
    exit( 1 );
    }
  */
  /********************/

#ifdef DEBUG
  if( S->num[NUM_SOL] - prevsize != nMIS ){
    fprintf( stderr, "error@LS2Core: theoretical MIS size (%d) != number of nodes that are added (%d)\n", nMIS, S->num[NUM_SOL]-prevsize );
    exit( EXIT_FAILURE );
  }
#endif

  free( U );
  free( O );
  return TRUE;
}



#if 0

void perturb_ver2( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean firstFlag=TRUE,*cov;
  double breakprob=IP->breakprob;
  Node **N,*u;
  Dim d;
  int i,j,h,k,K,nS0,n123,ub,*Prev,numPrev=0;
 
  /*** if no improvement is made for a long time, perform switching ***/
  /*
  if( ( IR->t - IR->lastImproved ) % (S->Nodes) == 0 ){
    if( genrand_real2() < 0.5 )
      lineSwitch( IP, IR, S );
    else
      symbolSwitch( IP, IR, S );
      }
  */

  /*** kick ***/
  ub = S->n * (S->n+1) / 2;
  cov = (Boolean*)malloc(ub*sizeof(Boolean));
  Prev = (int*)malloc(ub*sizeof(int));
  while( 1 ){
    nS0 = S->num[NUM_SOL] + S->num[NUM_FREE];
    n123 = S->Nodes - nS0;
    for(h=0;h<ub;h++)
      cov[h] = FALSE;
    for(i=0;i<S->num[NUM_SOL];i++){
      h = S->Pi[i]->pos[ROW] * (S->Pi[i]->pos[ROW]+1) / 2 + S->Pi[i]->pos[COL];
      cov[h] = TRUE;
    }
    for(h=0;h<numPrev;h++)
      cov[Prev[h]] = TRUE;
    N = (Node**)malloc(n123*sizeof(Node*));
    K = 0;
    i = 0;
    j = -1;
    for(h=0;h<ub;h++){
      j++;
      if( j>i ){
	i++;
	j=0;
      }
      if( cov[h] )
	continue;
      for(k=0;k<S->n;k++)
	if( S->C[i][j][k] != NULL ){
	  N[K] = S->C[i][j][k];
	  K++;
	}
    }
    if( K == 0 )
      break;
    if( firstFlag ){
      firstFlag = FALSE;
      u = getMostPenalizedNode( N, K );
    }
    else
      u = N[ (int)( genrand_real2() * (double)K ) ];
    free( N );	 
    /* drop solution neighbors of u */
    numPrev = 0;
    for(d=0;d<DIM;d++)
      if( u->solneigh[d] != NULL ){
	i = u->solneigh[d]->pos[ROW];
	j = u->solneigh[d]->pos[COL];
	Prev[numPrev] = i*(i+1)/2 + j;
	numPrev++;
	drop( S, u->solneigh[d]->order );
      }
    /* insert u */
    insert( S, u->order );
    /* if free nodes appear, insert them */
    while( S->num[NUM_FREE] > 0 )
      insert( S, S->Pi[ S->num[NUM_SOL] + (int)( genrand_real2() * (double)S->num[NUM_FREE]) ]->order );
    /* break */
    if( genrand_real2() < breakprob || S->num[NUM_SOL] == S->UB )
      break;
  }
  free( cov );
  free( Prev );
  /*** update the penalties ***/
  for(i=0;i<S->num[NUM_SOL];i++)
    S->Pi[i]->penalty = IR->t;
}

void perturb_ver1( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean firstFlag=TRUE,*cov;
  double breakprob=DEFAULT_BREAKPROB;
  Node **N,*u;
  Dim d;
  int i,j,h,k,K,nS0,n123,ub,Prev[DIM],numPrev=0;

  ub = S->n * (S->n+1) / 2;
  cov = (Boolean*)malloc(ub*sizeof(Boolean));

  while( 1 ){
    nS0 = S->num[NUM_SOL] + S->num[NUM_FREE];
    n123 = S->Nodes - nS0;
    for(h=0;h<ub;h++)
      cov[h] = FALSE;
    for(i=0;i<S->num[NUM_SOL];i++){
      h = S->Pi[i]->pos[ROW] * (S->Pi[i]->pos[ROW]+1) / 2 + S->Pi[i]->pos[COL];
      cov[h] = TRUE;
    }
    for(h=0;h<numPrev;h++)
      cov[Prev[h]] = TRUE;
    N = (Node**)malloc(n123*sizeof(Node*));
    K = 0;
    i = 0;
    j = -1;
    for(h=0;h<ub;h++){
      j++;
      if( j>i ){
	i++;
	j=0;
      }
      if( cov[h] )
	continue;
      for(k=0;k<S->n;k++)
	if( S->C[i][j][k] != NULL ){
	  N[K] = S->C[i][j][k];
	  K++;
	}
    }
    if( K == 0 ){
      printf("example\n");
      break;
    }
    if( firstFlag ){
      firstFlag = FALSE;
      u = getMostPenalizedNode( N, K );
    }
    else
      u = N[ (int)( genrand_real2() * (double)K ) ];
    free( N );	 

    /* drop solution neighbors of u */
    numPrev = 0;
    for(d=0;d<DIM;d++)
      if( u->solneigh[d] != NULL ){
	i = u->solneigh[d]->pos[ROW];
	j = u->solneigh[d]->pos[COL];
	Prev[numPrev] = i*(i+1)/2 + j;
	numPrev++;
	drop( S, u->solneigh[d]->order );
      }
    /* insert u */
    insert( S, u->order );
    /* if free nodes appear, insert them */
    while( S->num[NUM_FREE] > 0 )
      insert( S, S->Pi[ S->num[NUM_SOL] + (int)( genrand_real2() * (double)S->num[NUM_FREE]) ]->order );
    /* break */
    if( genrand_real2() < breakprob || S->num[NUM_SOL] == S->UB )
      break;
  }

  free( cov );

  /*** update the penalties ***/
  for(i=0;i<S->num[NUM_SOL];i++)
    S->Pi[i]->penalty = IR->t;
}


void perturb_ver0( ILSParam *IP, ILSRec *IR, Solution *S ){
  Boolean firstFlag=TRUE,*cov;
  double breakprob=DEFAULT_BREAKPROB;
  Node **N,*Prev[DIM],*x,*y,*u,*v;
  Dim root_d,d,d_,dx,dy;
  Dir dir;
  int i,j,imax,h,k,p,K,nS0,n123,numPrev,ub;

  /*** Step 1:
       exchange 1-tight node and its only solution neighbor
       this is done at most imax times        
       NOTE: free node does not appear due to 1-maximality ***/
  imax = S->num[NUM_1TIGHT];
  for(i=0;i<imax;i++){
    if( genrand_real2() < breakprob )
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

  /*** Step 2: exchange 2-tight nodes v,w and their common solution neighbors x,y ***/
  imax = S->num[NUM_2TIGHT];
  for(i=0;i<imax;i++){
    if( genrand_real2() < breakprob )
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
 
  /*** Step 3: insert non-free nodes forcibly ***/
  while( 1 ){
    nS0 = S->num[NUM_SOL] + S->num[NUM_FREE];
    n123 = S->Nodes - nS0;
    if( firstFlag ){
      firstFlag = FALSE;
      ub = S->n * (S->n+1) / 2;
      cov = (Boolean*)malloc(ub*sizeof(Boolean));
      for(h=0;h<ub;h++)
	cov[h] = FALSE;
      for(i=0;i<S->num[NUM_SOL];i++){
	h = S->Pi[i]->pos[ROW] * (S->Pi[i]->pos[ROW]+1) / 2 + S->Pi[i]->pos[COL];
	cov[h] = TRUE;
      }
      N = (Node**)malloc(n123*sizeof(Node*));
      K = 0;
      i = 0;
      j = -1;
      for(h=0;h<ub;h++){
	j++;
	if( j>i ){
	  i++;
	  j=0;
	}
	if( cov[h] )
	  continue;
	for(k=0;k<S->n;k++)
	  if( S->C[i][j][k] != NULL ){
	    N[K] = S->C[i][j][k];
	    K++;
	  }
      }
      u = getMostPenalizedNode( N, K );
      free( cov );
      free( N );	

      /*
	if( S->num[NUM_1TIGHT] == 0 )
	u = getMostPenalizedNode( S->Pi+nS0, n123 );
	else{
	N = (Node**)malloc(n123*sizeof(Node*));
	cov = (Boolean*)malloc(S->Nodes*sizeof(Boolean));
	for(p=nS0;p<S->Nodes;p++)
	cov[p] = FALSE;
	for(i=0;i<S->num[NUM_1TIGHT];i++){
	u = S->Pi[nS0+i];
	for(d=0;d<DIM;d++)
	if( u->solneigh[d] != NULL )
	break;
	x = u->solneigh[d];
	root_d = ROW;
	d_ = root_d;
	dir = FORWARD;
	v = x;
	while( 1 ){
	v = getNextNeighbor( x, v, &root_d, &d_, &dir );
	if( v == NULL )
	break;
	cov[v->order] = TRUE;		
	}
	}
	K = 0;
	for(p=nS0;p<S->Nodes;p++)
	if( cov[p] ){
	N[K] = S->Pi[p];
	K++;
	}
	u = getMostPenalizedNode( N, K );
	free( N );
	free( cov );
	}
      */
    }
    else{
      K = 0;
      if( numPrev ){
	x = Prev[ (int)( genrand_real2() * (double)numPrev ) ];	
	N = (Node**)malloc(S->n*sizeof(Node*));
	d_ = VAL;
	dir = FORWARD;
	v = x;
	while( 1 ){
	  v = getNextNeighborDim( x, v, VAL, &d_, &dir );
	  if( v == NULL )
	    break;
	  if( v->isSol )
	    continue;
	  N[K] = v;
	  K++;	
	}
      }
      if( K == 0 )
	u = S->Pi[ nS0 + (int)( genrand_real2() * (double)n123 ) ];
      else
	u = N[ (int)( genrand_real2() * (double)K ) ];
      if( numPrev )
	free( N );
    }

    /* drop solution neighbors of u */
    numPrev = 0;
    for(d=0;d<DIM;d++)
      if( u->solneigh[d] != NULL ){
	if( d != VAL ){
	  Prev[numPrev] = u->solneigh[d];
	  numPrev++;
	}
	drop( S, u->solneigh[d]->order );
      }
    /* insert u */
    insert( S, u->order );
    /* if free nodes appear, insert them */
    while( S->num[NUM_FREE] > 0 )
      insert( S, S->Pi[ S->num[NUM_SOL] + (int)( genrand_real2() * (double)S->num[NUM_FREE]) ]->order );
    /* break */
    if( genrand_real2() < breakprob )
      break;
  }

  /*** update the penalties ***/
  for(i=0;i<S->num[NUM_SOL];i++)
    S->Pi[i]->penalty = IR->t;
}
#endif
