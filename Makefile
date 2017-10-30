SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

HOSTNAME := $(shell hostname)

ifeq ("$(HOSTNAME)","phdstulin2.isye.gatech.edu")
	HPCOPTPATH    = /hpcopt
	IBMPATH       = ibm
	CPLEX_VERSION = CPLEX_Studio126
endif

ifeq ("$(HOSTNAME)","brownflaming-VirtualBox")
	HPCOPTPATH    = /opt
	IBMPATH       = ibm
	CPLEX_VERSION = CPLEX_Studio1263
endif

ifeq ("$(HOSTNAME)","Jikais-MacBook-Pro.local")
	SYSTEM = x86-64_osx
	HPCOPTPATH    = /Users/brownflaming/Applications
	IBMPATH       = IBM
	CPLEX_VERSION = CPLEX_Studio1263
endif

ifeq ("$(HOSTNAME)","tesla2.isye.gatech.edu")
	CPLEXDIR      = /home/dzink3/my_cplex/cplex
	CONCERTDIR    = /home/dzink3/my_cplex/concert
else
	CPLEXDIR      = $(HPCOPTPATH)/$(IBMPATH)/ILOG/$(CPLEX_VERSION)/cplex
	CONCERTDIR    = $(HPCOPTPATH)/$(IBMPATH)/ILOG/$(CPLEX_VERSION)/concert
endif

# --------------------------------------------------------------------
# Compiler selection and option
# --------------------------------------------------------------------

ifeq ("$(HOSTNAME)","Jikais-MacBook-Pro.local")
	CCC = clang++ -O0
	CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -stdlib=libstdc++
else
	CCC = g++ -std=c++0x
	CCOPT = -m64 -O2 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
endif

# -------------------------------------------------------------------
# Link options and libraries
# -------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)

ifeq ("$(HOSTNAME)","Jikais-MacBook-Pro.local")
	CCLNFLAGS = -lconcert -lilocplex -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit -stdlib=libstdc++
else
	CCLNFLAGS = -lilocplex -lcplex -lconcert -lm -pthread -g -fopenmp
endif
#-g -fopenmp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -g -O0 -fopenmp
#-g -O0 -fopenmp

OBJ_DIR = obj
OUT_RELEASE = bin/msuc

OBJ_FILE = $(OBJ_DIR)/mt64.o $(OBJ_DIR)/functions.o $(OBJ_DIR)/sddip.o

print-%: ; @echo $*=$($*)

all: sddip

debug: CCLNFLAGS += -DDEBUG -g
debug: CCFLAGS += -DDEBUG -g
debug: sddip

sddip: $(OBJ_FILE)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(OUT_RELEASE) $(OBJ_FILE) $(CCLNFLAGS)

$(OBJ_DIR)/sddip.o: sddip.cpp
	$(CCC) -c $(CCFLAGS) sddip.cpp -o $(OBJ_DIR)/sddip.o

$(OBJ_DIR)/mt64.o: mt64.cpp
	$(CCC) -c $(CCFLAGS) mt64.cpp -o $(OBJ_DIR)/mt64.o

$(OBJ_DIR)/functions.o: functions.cpp
	$(CCC) -c $(CCFLAGS) functions.cpp -o $(OBJ_DIR)/functions.o

clean:
	/bin/rm -rf $(OBJ_DIR)/*.o *~ *.class
	/bin/rm -rf $(OUT_RELEASE)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp
