CXX=g++
COMPFLAGS=-lm -Wall -Wpedantic -Winline -Wno-long-long -O3 -std=c++17 -lz

all: count_min_sketch
	$(CXX) count_min_sketch.o CMA.cpp  $(COMPFLAGS) -o cma

count_min_sketch:
	$(CXX) -c count_min_sketch.c -o count_min_sketch.o $(COMPFLAGS)

clean: 
	rm cma;
	rm *.o;