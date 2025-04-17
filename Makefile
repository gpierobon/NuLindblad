CC          = g++
OPTIMIZE    = -fopenmp -O3 -g -Wall -Wno-unknown-pragmas -march=haswell

INCL        = -I./src
EIGEN_INCL  = -I./externals

EXEC        = run 

CFLAGS      = $(OPTIMIZE) $(EIGEN_INCL) $(INCL)

SRC         = src/main.cxx
OBJ         = $(SRC:.cxx=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC)

%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)

.PHONY: all clean

