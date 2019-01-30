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

/*
 * heap.h -- Pairing heap definitions
 */

/*
 * Copyright 1989, 1996 by John Kececioglu
 */


#ifndef CompoundHeapInclude
#define CompoundHeapInclude


#include "portable.h"
//#include "heap.h"

/*
 * Nil pointer
 */
#ifndef Nil
#define Nil 0
#endif

/*
 * Heap key and value
 */



#ifndef HeapBlockSize
#define HeapBlockSize 16 /* Number of heaps allocated per memory request */
#define NodeBlockSize 32 /* Number of nodes allocated per memory request */
#endif

/*
 * Heap node
 */
typedef struct CompoundHeapStruct {
    CompoundKey   Key;
    Pointer Value;
    struct CompoundHeapStruct *Up, *Left, *Right;
       /*
        * `Up' is reused for the pool of free nodes
        */
} CompoundHeapNode;

/*
 * Pairing heap
 */
typedef struct {
    CompoundHeapNode *Root;
} CompoundHeap;


/*
 * Creation and destruction
 */
extern CompoundHeap *CreateCompoundHeap  Proto(());
extern Void  DestroyCompoundHeap Proto(( CompoundHeap *H ));

/*
 * Insertion and deletion
 */
extern CompoundHeapNode *CompoundHeapInsert Proto(( CompoundKey K, Pointer V, CompoundHeap *H ));
extern Void      CompoundHeapDelete Proto(( CompoundHeapNode *N, CompoundHeap *H ));

/*
 * Arbitrary access
 */

extern CompoundKey   CompoundHeapNodeKey   Proto(( CompoundHeapNode *N ));
extern Pointer CompoundHeapNodeValue Proto(( CompoundHeapNode *N ));

/*
 * Finding the minimum
 */
extern CompoundKey   CompoundHeapMinKey    Proto(( CompoundHeap *H ));
extern Pointer CompoundHeapMinValue  Proto(( CompoundHeap *H ));
extern Pointer CompoundHeapDeleteMin Proto(( CompoundHeap *H ));

/*
 * Miscellaneous
 */
extern CompoundHeap *CompoundHeapUnion       Proto(( CompoundHeap *A, CompoundHeap *B ));
extern Void  CompoundHeapDecreaseKey Proto(( CompoundHeapNode *N, CompoundKey K, CompoundHeap *H ));
extern short CompoundHeapIsEmpty     Proto(( CompoundHeap *H ));


#endif /* HeapInclude */
