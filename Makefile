INC=-Iexternal
MPI_CC  = mpicxx
CFLAGS  = -O3 -lm -std=c++98

all: task2

task2: task2.cpp *.h
	$(MPI_CC) $(CFLAGS) task2.cpp -o task2

clean: 
	rm -f task2