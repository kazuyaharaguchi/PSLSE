########## general macros ##########
CC = gcc
CCC = g++
CFLAGS = -O2 -Wall
RM = rm -f
SRC = symlatin.c misc.c mt19937ar.c
OBJ = $(SRC:%.c=%.o)
LSC = greedy.c elementary.c localsearch.c cpu_time.c
LSO = $(LSC:%.c=%.o)
CPPOBJ = $(SRC:%.c=%.cpp.o) $(LSC:%.c=%.cpp.o)
ALL = SYMILS SYMCPLEX SYMCPOPT SYMLSSOL HATFEAS HATFEAS_LS HATFEASFAST
CLEAN = *.o *~ \\\#*


########## CPLEX macros ##########
LONG_BIT = $(shell getconf LONG_BIT)
ILOGDIR   = /usr/ilog/CPLEX_Studio126
LIBFORMAT  = static_pic
CPLEXDIR      = $(ILOGDIR)/cplex
CONCERTDIR    = $(ILOGDIR)/concert
CPOPTDIR      = $(ILOGDIR)/cpoptimizer
ifeq ($(LONG_BIT),32)
	SYSTEM = x86_linux
else
	SYSTEM = x86-64_linux
endif
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINCDIR   = $(CPLEXDIR)/include
CPLEX_COPT    = -m$(LONG_BIT) -fPIC -fno-strict-aliasing 
CPLEX_CFLAGS  = $(CPLEX_COPT) -I$(CPLEXINCDIR) -Wall
CPLEX_CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -m$(LONG_BIT) -lm -pthread
CPOPT_CFLAGS = -DIL_STD -O -DNDEBUG -I $(CPOPTDIR)/include -I$(CONCERTDIR)/include -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -Wno-long-long -ffloat-store -m$(LONG_BIT) -DILOUSEMT -D_REENTRANT -DILM_REENTRANT # -march=i686 is required for 32bit? 
LDFLAGS = -L$(CPOPTDIR)/lib/$(SYSTEM)/static_pic -lcp -L$(CPLEXDIR)/lib/$(SYSTEM)/static_pic -lcplex -L$(CONCERTDIR)/lib/$(SYSTEM)/static_pic -lconcert  -lpthread 

########## generating executable files ##########
clean:
	$(RM) $(CLEAN)

distclean:
	$(RM) $(ALL) $(CLEAN)

SYMILS:    main.o define.h $(LSO) $(OBJ)
	$(CC) $(CFLAGS) -o $@ main.o $(LSO) $(OBJ)

SYMCPLEX: symCPLEX.o symcplex-lib.o $(LSO) $(OBJ)
	$(CC) -o $@ symCPLEX.o symcplex-lib.o $(LSO) $(OBJ) $(CPLEX_CFLAGS) $(CPLEX_CLNFLAGS)

SYMCPOPT: symCPOPT.cpp.o $(CPPOBJ)
	$(CCC) -o $@ symCPOPT.cpp.o $(CPPOBJ) $(LDFLAGS)

SYMLSSOL:   symLSSOL.o define.h $(LSO) $(OBJ)
	$(CC) $(CFLAGS) -o $@ symLSSOL.o $(LSO) $(OBJ)

HATFEAS: hatfeas.o symcplex-lib.o $(OBJ)
	$(CC) -o $@ hatfeas.o symcplex-lib.o $(OBJ) $(CPLEX_CFLAGS) $(CPLEX_CLNFLAGS)

HATFEAS_LS: hatfeas_ls.o symcplex-lib.o $(LSO) $(OBJ)
	$(CC) -o $@ hatfeas_ls.o symcplex-lib.o $(LSO) $(OBJ) $(CPLEX_CFLAGS) $(CPLEX_CLNFLAGS)

HATFEASFAST: hatfeasfast.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ hatfeasfast.o $(OBJ) -lm

DEBUG: hatfeas-debugger.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ hatfeas-debugger.o $(OBJ) -lm

ALL:	$(ALL)

all:	$(ALL)

########## generating object files ##########
hatfeas.o: hatfeas.c symcplex-lib.h
	$(CC) -o $@ -c $< $(CPLEX_CFLAGS)

hatfeas_ls.o: hatfeas_ls.c symcplex-lib.h
	$(CC) -o $@ -c $< $(CPLEX_CFLAGS)

symCPLEX.o: symCPLEX.c symcplex-lib.h
	$(CC) -o $@ -c $< $(CPLEX_CFLAGS)

symcplex-lib.o: symcplex-lib.c symcplex-lib.h
	$(CC) -o $@ -c $< $(CPLEX_CFLAGS)


symCPOPT.cpp.o: symCPOPT.cpp
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)


symlatin.cpp.o: symlatin.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)

misc.cpp.o: misc.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)

mt19937ar.cpp.o: mt19937ar.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)

greedy.cpp.o: greedy.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)

elementary.cpp.o: elementary.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)

localsearch.cpp.o: localsearch.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)

cpu_time.cpp.o: cpu_time.c
	$(CCC) $(CPOPT_CFLAGS) -c $< -o $@ $(LDFLAGS)


.c.o:
	$(CC) $(CFLAGS) -c $<
