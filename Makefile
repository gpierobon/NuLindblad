CC          = g++
OPTIMIZE    = -fopenmp -O3 -g -Wall -Wno-unknown-pragmas -march=haswell

INCL        = -I./src
EIGEN_INCL  = -I./externals

EXEC        = run 
EXEC_X      = runx 

CFLAGS      = $(OPTIMIZE) $(EIGEN_INCL) $(INCL)

SRC         = src/main.cxx
SRC_X       = src/main_jx.cxx
OBJ         = $(SRC:.cxx=.o)
OBJ_X       = $(SRC_X:.cxx=.o)

all: $(EXEC) $(EXEC_X)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC)

$(EXEC_X): $(OBJ_X)
	$(CC) $(CFLAGS) $(OBJ_X) -o $(EXEC_X)

%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(OBJ_X) $(EXEC) $(EXEC_X)

.PHONY: all clean

