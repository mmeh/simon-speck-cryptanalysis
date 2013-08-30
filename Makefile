CC=g++
CFLAGS=-std=c++11 -fopenmp -Wall -Wno-format -g -O3

all:
	$(CC) $(CFLAGS) -o simon simon.cpp

clean:
	rm -rf *.o simon