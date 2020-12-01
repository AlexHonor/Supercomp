INC=-Iexternal
MPI_CC  = mpicxx
LDFLAGS = -lm
CFLAGS  = -Ofast

all: task2

task2: task2.cpp *.h
	$(MPI_CC) $(LDFLAGS) task2.cpp -o task2 -std=c++98 -pedantic

clean: 
	rm -f task2