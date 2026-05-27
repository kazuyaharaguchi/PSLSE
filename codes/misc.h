/*** misc.h ***/

char **split( char *s, char c );
int getNumChar( char *s, char c );
int *qsortWi( void *array, size_t num, size_t size, 
	      int cmp( const void *, const void *) );
void qsortWiSub( void *array, int *I, int head, int tail, size_t size,
		 int cmp( const void *, const void *) );
void shuffle( void *array, size_t num, size_t size, double yourrand() );
void swap( void *array, int i, int j, size_t size );
void *mallocE( size_t size );
int cmpInt( const void *p, const void *q );
int cmpNegInt( const void *p, const void *q );
int cmpDouble( const void *p, const void *q );
int smaller( int a, int b );
int larger( int a, int b );
void initCombi( int n, int k, int *A );
int getNextCombi( int n, int k, int *A );
int isCanonicalForm( int n, int k, int *A );
