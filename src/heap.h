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


#ifndef HeapInclude
#define HeapInclude


#include "portable.h"


/*
 * Nil pointer
 */
#ifndef Nil
#define Nil 0
#endif

/*
 * Heap key and value
 */
typedef distancevalue HeapKey; /* Changed by AW. Was typedef Pointer */
typedef Pointer HeapValue;


#define HeapBlockSize 16 /* Number of heaps allocated per memory request */
#define NodeBlockSize 32 /* Number of nodes allocated per memory request */

/*
 * Heap node
 */
typedef struct HeapStruct {
    HeapKey   Key;
    HeapValue Value;
    struct HeapStruct *Up, *Left, *Right;
       /*
        * `Up' is reused for the pool of free nodes
        */
} HeapNode;

/*
 * Pairing heap
 */
typedef struct {
   int (*Compare) Proto(( HeapKey, HeapKey )); /* Key comparison-function */
   HeapNode *Root;
      /*
       * `Root' is reused for the pool of free heaps
       */
} Heap;


/*
 * Creation and destruction
 */
extern Heap *CreateHeap  Proto(( int (*Compare) Proto((HeapKey, HeapKey)) ));
extern Void  DestroyHeap Proto(( Heap *H ));

/*
 * Insertion and deletion
 */
extern HeapNode *HeapInsert Proto(( HeapKey K, HeapValue V, Heap *H ));
extern Void      HeapDelete Proto(( HeapNode *N, Heap *H ));

/*
 * Arbitrary access
 */
extern HeapKey   HeapNodeKey   Proto(( HeapNode *N ));
extern HeapValue HeapNodeValue Proto(( HeapNode *N ));

/*
 * Finding the minimum
 */
extern HeapKey   HeapMinKey    Proto(( Heap *H ));
extern HeapValue HeapMinValue  Proto(( Heap *H ));
extern HeapValue HeapDeleteMin Proto(( Heap *H ));

/*
 * Miscellaneous
 */
extern Heap *HeapUnion       Proto(( Heap *A, Heap *B ));
extern Void  HeapDecreaseKey Proto(( HeapNode *N, HeapKey K, Heap *H ));
extern short HeapIsEmpty     Proto(( Heap *H ));


#endif /* HeapInclude */
