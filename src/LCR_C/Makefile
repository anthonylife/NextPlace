CC = g++
CFLAGS = -Wall -O3 -fopenmp -Wcpp 

all: train

train: run.cpp ../utils.hpp model.hpp
	$(CC) $(CFLAGS) -o run  run.cpp 

train1: run.cpp ../utils.hpp model1.hpp
	$(CC) $(CFLAGS) -o run1  run.cpp 


clean:
	rm run
	rm run1
