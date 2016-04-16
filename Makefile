# Paths

TRDF_HOME = $(CURDIR)
SRC = $(TRDF_HOME)/src
LIB = $(TRDF_HOME)/lib
TES = $(TRDF_HOME)/tests
SOL = $(SRC)/solvers
BIN = $(TRDF_HOME)/bin
OBJ = $(TRDF_HOME)/objects

ifndef PROBLEM
   PROBLEM = tests/examples/trdf_main.f
endif

# Compiler options

CC  = gcc
FC  = gfortran
FCC = "-xf77-cpp-input"

# Solver configuration parameters

SOLVERLIB = /opt/tango/algencan-3.0.0/lib
SOLVER_INTERFACE = algencan_solver
SLOPTS = -lalgencan -ltrdf

# Linking options

LOPTS = $(SLOPTS)

export

all: lib solver

# Generate the main TRDF library
lib:
	$(MAKE) -C $(SRC) all install

# Generate the solver interface object
solver:
	$(MAKE) -C $(SOL) install

# User-defined executable
trdf: all
	$(FC) -L$(SOLVERLIB) -L$(LIB) $(OBJ)/solver.o \
	$(FCC) -I$(SRC) $(PROBLEM) $(LOPTS) -o $(BIN)/$@

# User-defined C executable
c_trdf: all
	$(CC) -L$(SOLVERLIB) -L$(LIB) $(OBJ)/solver.o \
	$(PROBLEM) $(LOPTS) -lgfortran -lm -o $(BIN)/$@

# Hock-Schittkowski test set executable
hstests: $(SOLVER_INTERFACE).o hstests.o TRDF.o functions.o libhs.a
	gfortran -L$(SOLVERLIB) $^  $(LOPTS) -o $@ 

clean:
	rm -vf $(LIB)/* $(BIN)/* $(OBJ)/*
	$(foreach i,$(SRC) $(TES),$(MAKE) -C $(i) clean;)

.PHONY: lib all clean
