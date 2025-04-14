CC          = g++
OPTIMIZE    = -fopenmp -O3 -g -Wall -Wno-unknown-pragmas -march=native

INCL        = -I./
EIGEN_INCL  = -I/home/g/work/OQS/largeN/eigen

EXEC        = run 

CFLAGS      = $(OPTIMIZE) $(EIGEN_INCL) $(INCL)

SRC         = main.cxx
OBJ         = $(SRC:.cxx=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC)

%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)

.PHONY: all clean

