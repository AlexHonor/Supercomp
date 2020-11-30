INC=-Iexternal

all: task2

task2: task2.cpp *.h
	mpic++ task2.cpp -o task2 -std=c++98 -pedantic

clean: 
	rm -f task2