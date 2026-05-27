#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "misc.h"



/***** split a string into strings *****/
char **split( char *s, char c ){
  typedef struct list{
    int i;
    struct list *p;
  } List;
  char **str;
  List head,*tail,*ip;
  int countnum=0,counts=0;
  int i=0,j=0,i0=0;
  head.i=-1;
  head.p=NULL;
  tail=&head;
  while(s[i]!='\0'){
    if(s[i]==c){
      ip=(List*)mallocE(sizeof(List));
      ip->i=countnum;
      ip->p=NULL;
      tail->p=ip;
      tail=ip;
      countnum=-1;
      counts++;
    }
    countnum++;
    i++;
  }
 ip=(List*)mallocE(sizeof(List));
  ip->i=countnum;
  ip->p=NULL;
  tail->p=ip;
  tail=ip;
  str=(char**)mallocE((counts+1)*sizeof(char*));
  ip=&head;
  do{
    int k=0;
    ip=ip->p;
    str[j]=(char*)mallocE(((ip->i)+1)*sizeof(char));
    for(i=i0;i<i0+ip->i;i++){
      str[j][k]=s[i];
      k++;
    }
    str[j][k]='\0';
    i0+=ip->i+1;
    j++;
  }while(ip->p!=NULL);
  // free memory
  tail = head.p;
  while( tail != NULL ){
    ip = tail->p;
    free( tail );
    tail = ip;
  }
  return str;

}

/***** number of char c in string s *****/
int getNumChar( char *s, char c ){
  int numChar = 0;
  register int k = 0;
  while( s[k] != '\0' ){
    if( s[k] == c )
      numChar++;
    k++;
  }
  return numChar;
}

/***** quicksort with indices *****/
int *qsortWi( void *array, size_t num, size_t size, 
	      int cmp( const void *, const void * ) ){
  int *I;
  register int i;
  // prepare index set
  I = (int*)malloc(num*sizeof(int));
  for(i=0;i<(int)num;i++)
    I[i] = i;
  // call subroutine
  qsortWiSub( array, I, 0, num-1, size, cmp ); 	
  return I;
}


/***** sorting routine of qsortWi *****/
void qsortWiSub( void *array, int *I, int head, int tail, size_t size,
		 int cmp( const void *, const void * ) ){
  int uniformity = 1;
  register int i,j;
  // determine pivot
  for(i=head+1;i<=tail;i++)
    if( cmp( (char*)array+head*size, (char*)array+i*size ) < 0 ){
      uniformity = 0;
      break;
    }
    else if( cmp( (char*)array+head*size, (char*)array+i*size ) > 0 ){
      swap( array, head, i, size );
      swap( I, head, i, sizeof(int) );
      uniformity = 0;
      break;
    }
  if( uniformity )
    return;
  // exchange values
  i = head+1;
  j = tail;
  while(1){
    while( cmp( (char*)array+i*size, (char*)array+head*size ) <= 0 && i<j )
      i++;
    while( cmp( (char*)array+j*size, (char*)array+head*size ) > 0  && i<j )
      j--;
    if( i >= j )
      break;
    swap( array, i, j, size );
    swap( I, i, j, sizeof(int) );
  }
  // recursive call
  qsortWiSub( array, I, head, i-1, size, cmp ); 	
  qsortWiSub( array, I, i, tail, size, cmp ); 	
}


/***** shuffle array (yourrand has to return real in [0,1)) *****/
void shuffle( void *array, size_t num, size_t size, double yourrand() ){
  register int i,j;
  for(i=0;i<(int)num;i++){
    j = (int)((double)i + yourrand()*(double)(num-i));
    swap( array, i, j, size );
  }	
}


/***** swap values in array *****/
void swap( void *array, int i, int j, size_t size ){
  size_t r;
  char tmp;
  for(r=0;r<size;r++){
    tmp = ((char*)array)[i*size+r];
    ((char*)array)[i*(int)size+r] = ((char*)array)[j*(int)size+r];
    ((char*)array)[j*(int)size+r] = tmp;
  }
}


/***** functions for memory allocation *****/
void *mallocE( size_t size ){
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc: not enough memory.\n" );
    exit( EXIT_FAILURE );
  }
  return s;
}


/***** compare integers *****/
int cmpInt( const void *p, const void *q ){
  if( *(int*)p < *(int*)q )
    return -1;
  if( *(int*)p > *(int*)q )
    return 1;
  return 0;
}

int cmpNegInt( const void *p, const void *q ){
  return -cmpInt( p, q );
}

/***** compare floating numbers *****/
int cmpDouble( const void *p, const void *q ){
  if( *(double*)p < *(double*)q )
    return -1;
  else if( *(double*)p > *(double*)q )
    return 1;
  return 0;
}


int smaller( int a, int b ){
  if( a<b )
    return a;
  return b;
}

int larger( int a, int b ){
  if( a>b )
    return a;
  return b;
}

void initCombi( int n, int k, int *A ){
  int i;
  for(i=0;i<k;i++)
    A[i] = i;
}

int getNextCombi( int n, int k, int *A ){
  int flag=0;
  int i,j;
  A[k-1]++;
  if( A[k-1] <= n-1 )
    return 1;
  for(i=k-2;i>=0;i--)
    if( A[i] < n-(k-i) ){
      A[i]++;
      for(j=i+1;j<k;j++)
	A[j] = A[j-1]+1;
      flag = 1;
      break;
    }
  return flag;
}


int isCanonicalForm( int n, int k, int *A ){
  int flag = 1;
  int i,max;
  max = A[0]+1;
  for(i=1;i<k;i++)
    if( A[i]-A[i-1] > max ){
      flag = 0;
      break;
    }
  return flag;
}
