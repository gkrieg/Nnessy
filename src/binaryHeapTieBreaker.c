/*
    This file is part of the dispersion tree code.

    The dispersion tree code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Foobar is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the code for a dispersion tree.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include "binaryHeapTieBreaker.h"
//#include "binaryHeap.h"


#define LEFT(n) (2*n+1)
#define RIGHT(n) (2*n+2)
#define PARENT(n) ( (n-1)/2)

// because of the tolerance I don't bother making a is greater-than-or-equal-to function (GE)

/*
#define TOLERANCE 0.0001
#define ISLT(H1, H2)   ( (H1.dPrim < H2.dPrim - TOLERANCE) || ( abs(H1.dPrim - H2.dPrim) <= TOLERANCE && H1.dSecond < H2.dSecond)  ) 
#define ISGT(H1, H2)   ( (H1.dPrim > H2.dPrim + TOLERANCE) || ( abs(H1.dPrim - H2.dPrim) <= TOLERANCE && H1.dSecond > H2.dSecond)  ) 

#define ISLT(H1, H2) ( abs(H1.dPrim - H2.dPrim) <= TOLERANCE ? H1.dSecond < H2.dSecond :  H1.dPrim <  H2.dPrim)
#define ISGT(H1, H2) ( abs(H1.dPrim - H2.dPrim) <= TOLERANCE ? H1.dSecond > H2.dSecond :  H1.dPrim >  H2.dPrim)
*/

#define ISLT(H1, H2) ( H1.dPrim ==  H2.dPrim ? H1.dSecond < H2.dSecond :  H1.dPrim <  H2.dPrim)
#define ISGT(H1, H2) ( H1.dPrim ==  H2.dPrim ? H1.dSecond > H2.dSecond :  H1.dPrim >  H2.dPrim)


BinaryHeapTiebreaker * 
initHeapTB(unsigned initialSize, char type) {

    BinaryHeapTiebreaker *init =  (BinaryHeapTiebreaker*) malloc(sizeof(BinaryHeapTiebreaker));
    void * memory = malloc(sizeof(Heapvaltb_distanceitem)*initialSize);
    if (memory==NULL || init==NULL) {
	fprintf(stderr, "not enough memory to init the heap!\n");
	exit(1);
    }

    init->currSize = 0;
    init->currMalloc = initialSize;
    init->heap = (Heapvaltb_distanceitem*) memory;
    if (type)
	init->type = MAXHEAP;
    else
	init->type = MINHEAP;
    
    return init;
}

void
destroyHeapTB(BinaryHeapTiebreaker *h) {
    free(h->heap);
    free(h);
}

void
maxHeapifyTB(BinaryHeapTiebreaker *h, unsigned index) {
    unsigned left = LEFT(index);
    unsigned right = RIGHT(index);
    unsigned largest;

//    if (left < h->currSize && h->heap[left].dPrim > h->heap[index].dPrim) 
    if (left < h->currSize && ISGT(h->heap[left], h->heap[index])  )
	largest = left;
    else
	largest = index;

//    if (right < h->currSize && h->heap[right].dPrim > h->heap[largest].dPrim)
    if (right < h->currSize && ISGT(h->heap[right], h->heap[largest]) )
	largest = right;

    if (largest != index) {
	Heapvaltb_distanceitem tmp;
	tmp = h->heap[index];
	h->heap[index] = h->heap[largest];
	h->heap[largest]=tmp;
	maxHeapifyTB(h, largest);
    }
}

void
minHeapifyTB(BinaryHeapTiebreaker *h, unsigned index) {
    unsigned left = LEFT(index);
    unsigned right = RIGHT(index);
    unsigned smallest;
//    if (left < h->currSize && h->heap[left].dPrim < h->heap[index].dPrim) 
    if (left < h->currSize && ISLT(h->heap[left], h->heap[index]) )
	smallest = left;
    else
	smallest = index;

//    if (right < h->currSize && h->heap[right].dPrim < h->heap[smallest].dPrim)
    if (right < h->currSize && ISLT(h->heap[right], h->heap[smallest]) )
	smallest = right;

    if (smallest != index) {
	Heapvaltb_distanceitem tmp;
	tmp = h->heap[index];
	h->heap[index] = h->heap[smallest];
	h->heap[smallest]=tmp;
	minHeapifyTB(h, smallest);
    }
}


// used to get the values of the best item in the heap
distance
getBestDistanceTB(BinaryHeapTiebreaker *heap) {
    if (heap->currSize < 1) {
	fprintf(stderr, "Heap underflow...\n");
	exit(1);
    } 
    return heap->heap[0].dPrim;
}



distance
getBestDistanceSecondaryTB(BinaryHeapTiebreaker *heap) {
    if (heap->currSize < 1) {
	fprintf(stderr, "Heap underflow...\n");
	exit(1);
    } 
    return heap->heap[0].dSecond;
}



item
getBestItemTB(BinaryHeapTiebreaker *heap) {
    if (heap->currSize < 1) {
	fprintf(stderr, "Heap underflow...\n");
	exit(1);
    } 
    
    return heap->heap[0].i;
}

// used to delete the 0th item in the heap

void
heapExtractTB(BinaryHeapTiebreaker *heap) {
    if (heap->currSize < 1) {
	fprintf(stderr, "Heap underflow... extract\n");
	exit(1);
    }

    --(heap->currSize);
    heap->heap[0] = heap->heap[ heap->currSize];
    if (heap->type)
	maxHeapifyTB(heap, 0);
    else
	minHeapifyTB(heap, 0);
}



// a generic implementation of heap-increase-key (a la cormen)
void
heapSetKeyTB(BinaryHeapTiebreaker *heap, int index, distance primval, distance secondVal, item i ) {
    heap->heap[index].dPrim = primval;
    heap->heap[index].dSecond = secondVal;
    heap->heap[index].i = i;
    Heapvaltb_distanceitem tmp;
    if (heap->type) { // maxheap
//	while (index >= 0 && heap->heap[ PARENT(index)].dPrim < heap->heap[index].dPrim) {
	while (index >= 0 && ISLT(heap->heap[ PARENT(index)], heap->heap[index] ) ) {
	    tmp = heap->heap[ PARENT(index)];
	    heap->heap[ PARENT(index)] = heap->heap[index];
	    heap->heap[index]  = tmp;
	    index = PARENT(index);
	}
    } else {
//	while (index >= 0 && heap->heap[ PARENT(index)].dPrim > heap->heap[index].dPrim) {
	while (index >= 0 && ISGT(heap->heap[ PARENT(index)], heap->heap[index])) {
	    tmp = heap->heap[ PARENT(index)];
	    heap->heap[ PARENT(index)] = heap->heap[index];
	    heap->heap[index]  = tmp;
	    index = PARENT(index);
	}
    }
}

void
heapInsertTB(BinaryHeapTiebreaker *h, distance dPrimary, distance dSecondary, item i) {
    if (h->currSize == h->currMalloc) { // grow the heap
	h->currMalloc = h->currMalloc*2;
	h->heap = (Heapvaltb_distanceitem*)realloc(h->heap, h->currMalloc * sizeof(Heapvaltb_distanceitem));
	if (h->heap == NULL) {
	    fprintf(stderr, "Not enough memory to grow the heap!\n");
	    exit(1);
	}
    }
    heapSetKeyTB(h, h->currSize, dPrimary, dSecondary, i);
    ++(h->currSize);    
}

/*
int
main (int argc, char **argv) {

    Heapvaltb_distanceitem d;
    Intpair x = {2,2};
    float n=6.0;
    d.i = x;
    d.dPrim = n;
    
    BinaryHeapTiebreaker *h = initHeapTB(100, MAXHEAP);

    int i;
    for (i=0; i < 100; ++i) {
	n = rand() % 20;
	heapInsertTB(h, n, rand() % 20, x);
//	fprintf(stderr, "%f %f\n",  getBestDistanceTB(h), n );
    }
    for (i=0; i < 100; ++i) {
	fprintf(stderr, "%f %f\n",  getBestDistanceTB(h), getBestDistanceSecondaryTB(h) );
	heapExtractTB(h);
    }

    return 0;

}

*/
