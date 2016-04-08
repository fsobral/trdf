# Paths

TRDF_HOME = $(CURDIR)
SRC = $(TRDF_HOME)/src
LIB = $(TRDF_HOME)/lib
TES = $(TRDF_HOME)/tests
SOL = $(SRC)/solvers
BIN = $(TRDF_HOME)/bin
OBJ = $(TRDF_HOME)/objects

ifndef USERMAIN
   USERMAIN = trdf_main
endif

# Compiler options
FC  = gfortran
FCC = "-xf77-cpp-input"

# Solver configuration parameters

SOLVERLIB = /opt/tango/algencan-3.0.0/lib
SOLVER_INTERFACE = algencan_solver
SLOPTS = -lalgencan

# Linking options

LOPTS = $(SLOPTS)

# Executable options

export

# User-defined executable
trdf:
	$(MAKE) -C $(SRC) all install
	$(MAKE) -C $(SOL) install
	$(MAKE) -C $(TES) $(USERMAIN)
	gfortran -L$(SOLVERLIB) -L$(LIB) $^ $(LOPTS) -ltrdf $(OBJ)/solver.o $(OBJ)/prob.o -o $(BIN)/$@

# User-defined C executable
c_trdf: trdf_main.c libtrdf.a
	gcc $^ -L$(SOLVERLIB) $(LOPTS) -lgfortran -lm -o $@

# Hock-Schittkowski test set executable
hstests: $(SOLVER_INTERFACE).o hstests.o TRDF.o functions.o libhs.a
	gfortran -L$(SOLVERLIB) $^  $(LOPTS) -o $@ 

clean:
	rm -vf $(LIB)/* $(BIN)/* $(OBJ)/*
	$(foreach i,$(SRC) $(TES),$(MAKE) -C $(i) clean;)
