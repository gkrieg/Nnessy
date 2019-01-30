#ifndef BINARY_HEAP_TIEBREAKER
#define BINARY_HEAP_TIEBREAKER

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

#include "portable.h"

// a twist on: 
// http://stackoverflow.com/questions/7069371/which-way-is-better-for-creating-type-agnostic-structures-in-c
// uses the C-prepocessor to let you modify the heap data types
// I am assuming a compound heap, where D is the *result* of a distance computation
// and I is the item involved


#define HEAP_TRIO(D,I) typedef struct Heapvaltb_##D##I { D dPrim; D dSecond;  I i; } Heapvaltb_##D##I

// change these types if you want to change the types involved with the heap.
#ifndef distance

typedef struct {
    unsigned x;
    unsigned y;
} Intpair;

#define distance distancevalue
#define item Intpair
#endif

HEAP_TRIO(distance, item);

// technically, I treat minheap as 0, and maxheap as everything else
#ifndef MINHEAP
#define MINHEAP 0
#define MAXHEAP 1
#endif

typedef struct {
    Heapvaltb_distanceitem *heap;
    unsigned currSize;
    unsigned currMalloc;
    char type;
} BinaryHeapTiebreaker;


BinaryHeapTiebreaker * 
initHeapTB(unsigned initialSize, char type);

void
destroyHeapTB(BinaryHeapTiebreaker *h);

distance
getBestDistanceTB(BinaryHeapTiebreaker *heap);

distance
getBestDistanceSecondaryTB(BinaryHeapTiebreaker *heap);

item
getBestItemTB(BinaryHeapTiebreaker *heap);

void
heapExtractTB(BinaryHeapTiebreaker *heap);

void
heapInsertTB(BinaryHeapTiebreaker *h, distance dPrimary, distance dSecondary, item i);
// TODO: batch create

#endif

