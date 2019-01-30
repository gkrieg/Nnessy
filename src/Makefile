include CompilerOptions.in

all: knn.c dispersionTree.o metricspaces.o word.o r32_point.o partition.o binaryHeapTieBreaker.o compoundHeap.o heap.o nnFinder.o lpgenerator.o lpgenerator

#CC = gcc -g -Wall -O3
CC = gcc -g -gstabs -Wall 
LIBS = -lm

module: knn.c dispersionTree.o metricspaces.o word.o r32_point.o partition.o binaryHeapTieBreaker.o compoundHeap.o heap.o

	$(CC) -o knn knn.c dispersionTree.o metricspaces.o word.o r32_point.o partition.o binaryHeapTieBreaker.o compoundHeap.o heap.o $(LIBS)

dispersionTree.o: dispersionTree.h dispersionTree.c
	$(CC) -c dispersionTree.c

metricspaces.o: metricspaces.h metricspaces.c
	$(CC) -c metricspaces.c

partition.o: partition.h partition.c
	$(CC) -c partition.c

binaryHeapTieBreaker.o: binaryHeapTieBreaker.h  binaryHeapTieBreaker.c
	$(CC) -c binaryHeapTieBreaker.c

nnFinder.o:  nnFinder.cpp
	g++ -g -gstabs -Wall -c  nnFinder.cpp

compoundHeap.o: compoundHeap.h compoundHeap.c
	$(CC) -c compoundHeap.c

heap.o: heap.h heap.c
	$(CC) -c heap.c	

# two example distance functions
word.o: word.h word.c
	$(CC) -c word.c

r32_point.o: r32_point.c r32_point.h
	$(CC) -c r32_point.c

lpgenerator:  main.o
	g++ -g -gstabs -Wall -std=c++0x -pthread -o lpgenerator lpgenerator2.o nnFinder.o word.o metricspaces.o binaryHeapTieBreaker.o compoundHeap.o heap.o r32_point.o partition.o dispersionTree.o main.o
	
lpgenerator.o: main.cpp lpgenerator2.cpp
	g++ -g -gstabs -Wall -std=c++0x -pthread -c main.cpp lpgenerator2.cpp 
 
clean:
	rm *.o lpgenerator
