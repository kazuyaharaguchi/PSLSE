/******************************
  localsearch.h
******************************/

/********** ILS-oriented functions **********/
void ILS( ILSParam *IP, ILSRec *IR, Solution *S );
void execILS1( ILSParam *IP, ILSRec *IR, Solution *S );
void execILS2( ILSParam *IP, ILSRec *IR, Solution *S );
Boolean updateBest( ILSParam *IP, ILSRec *IR, Solution *S );
void perturb( ILSParam *IP, ILSRec *IR, Solution *S );
void lineSwitch( ILSParam *IP, ILSRec *IR, Solution *S );
void symbolSwitch( ILSParam *IP, ILSRec *IR, Solution *S );
void symbolSwitchCore( ILSParam *IP, ILSRec *IR, Solution *S, int l, Node **V );

/********** LS-oriented functions **********/
Boolean LS1( ILSParam *IP, ILSRec *IR, Solution *S );
Boolean LS1Core( ILSParam *IP, ILSRec *IR, Solution *S, Node *x );
Boolean LS2( ILSParam *IP, ILSRec *IR, Solution *S );
Boolean LS2Core( ILSParam *IP, ILSRec *IR, Solution *S, Node *u );
