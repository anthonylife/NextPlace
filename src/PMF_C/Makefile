CC = g++
CFLAGS = -Wall -O3 -fopenmp -Wcpp 

all: train

train: run.cpp ../utils.hpp model.hpp
	$(CC) $(CFLAGS) -o run  run.cpp 

clean:
	rm run
