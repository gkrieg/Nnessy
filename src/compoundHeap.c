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
 * heap.c -- Two-pass pairing heaps
 */

/*
 * Copyright 1989, 1996 by John Kececioglu
 */


/*
 * Synopsis
 *
 * O(n) insert, meld, or decrease operations, and O(m) extract or delete
 * operations, are conjectured to take O(n + m log n) time worst-case.
 * A meld operation forms the union of two disjoint heaps, a decrease
 * operation decreases the key associated with an element, and an extract
 * operation deletes the element with the minimum key.
 *
 * References
 *
 * Fredman, Sedgewick, Sleator, and Tarjan, ``The pairing heap:
 * a new form of self-adjusting heap,'' Algorithmica 1, 111-129, 1986.
 *
 * Stasko and Vitter, ``Pairing heaps:  experiments and analysis,''
 * Communications of the ACM 30:3, 234-249, March 1987.
 *
 */

/*
 * Author
 *
 * John Kececioglu
 * kece@cs.uga.edu
 *
 * Department of Computer Science
 * The University of Georgia
 * Athens, GA 30602
 *
 */

/*
 * History
 *
 * 11 September 1998 JDK
 * Fixed a bug in the nonrecursive DestroyHeap that was pointed out by
 * Thomas Christof.  The code, on deleting nodes from a heap, was not ensuring
 * that the up-pointer of the root node is nil (as previously it did not need
 * to); now the root node of a heap always has a nil up-pointer.
 *
 * 10 September 1998 JDK
 * Wrote a nonrecursive version of DestroyHeap, on the urging of Thomas
 * Christof.  The previous recursive version could exceed the machine stack
 * limit when the height of the heap was large.
 *
 * 28 February 1996 JDK
 * Changed the interface to conform with library conventions.  Added block
 * allocation of nodes and heaps.
 *
 * 30 May 1986 JDK
 * Wrote the initial implementation.
 *
 */


#include <stdio.h>
#include "compoundHeap.h"


int
CompoundCompare(CompoundKey k1, CompoundKey k2) {

    if (k1.lowerboundDistance < k2.lowerboundDistance)
	return -1;
    else if (k1.lowerboundDistance > k2.lowerboundDistance)
	return 1;
    return 0;
/*
    if (k1.rank == k2.rank) {
	if (k1.expectedDistance < k2.expectedDistance)
	    return -1;
	else if (k1.expectedDistance > k2.expectedDistance)
	    return 1;
	return 0;
    } 
    return k1.rank - k2.rank;
*/
}



static CompoundHeapNode *CreateCompoundHeapNode
Proto(( CompoundKey K, Pointer V ));

static CompoundHeapNode *LinkCompoundHeap
Proto(( CompoundHeapNode *A, CompoundHeapNode *B));
         
static CompoundHeapNode *DeleteCompoundHeap
Proto(( CompoundHeapNode *N));


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Creation and destruction
 *
 */


static CompoundHeap *HeapPool = Nil; /* Pool of free heaps */
static CompoundHeapNode *NodePool = Nil; /* Pool of free nodes */

/*
 * Heap pool maintenance
 */
#define FreeCompoundHeap(H) (((H)->Root = (CompoundHeapNode *) HeapPool), HeapPool = (H))
#define NewCompoundHeap(H)  (((H) = HeapPool), HeapPool = (CompoundHeap *) HeapPool->Root)

/* 
 * Node pool maintenance
 */
#define FreeCompoundHeapNode(N) (((N)->Up = NodePool), NodePool = (N))
#define NewCompoundHeapNode(N)  (((N) = NodePool), NodePool = NodePool->Up)


int
HeapPairComporator(CompoundKey p1, CompoundKey p2) {
    if (p1.lowerboundDistance < p2.lowerboundDistance)
	return -1;
    else if (p1.lowerboundDistance > p2.lowerboundDistance)
	return 1;
    return 0;
/*
    if (p1.rank == p2.rank) {
	if (p1.expectedDistance < p2.expectedDistance)
	    return -1;
	else if (p1.expectedDistance > p2.expectedDistance)
	    return 1;
	return 0;
    } else
	return p1.rank - p2.rank;
*/

}



/*
 * CreateHeap -- Create an empty heap
 *
 */
CompoundHeap *
CreateCompoundHeap()
{
   register CompoundHeap *H;
   
   if (HeapPool == Nil)
   {
      register CompoundHeap *Block, *P;
      
      /*
       * Allocate a block of heaps
       */
      Block = (CompoundHeap *) Allocate(sizeof(CompoundHeap) * HeapBlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(CreateCompoundHeap) Memory allocation failed.\n");
         Halt();
      }

      /*
       * Place the heaps in the block into the pool
       */
      for (P = Block; P - Block < HeapBlockSize; P++)
         FreeCompoundHeap(P);
   }
   
   NewCompoundHeap(H);
   H->Root = Nil;
   
   return H;
}


/*
 * DestroyHeap -- Destroy a heap
 *
 */
Void DestroyCompoundHeap

#ifdef Ansi
   (CompoundHeap *H)
#else
   (H) CompoundHeap *H;
#endif

{
   register CompoundHeapNode *N, *P, *R, *L;

   N = H->Root;
   while (N != Nil)
   {
      P = N->Up;
      
      R = N->Right;
      if (R != Nil)
      {
         R->Up = P;
         P = R;
      }
      
      L = N->Left;
      if (L != Nil)
      {
         L->Up = P;
         P = L;
      }
      
      FreeCompoundHeapNode(N);
      N = P;
   }
   
   FreeCompoundHeap(H);
}


/*
 * CreateNode -- Create a new heap node
 *
 */
static CompoundHeapNode *CreateCompoundHeapNode

#ifdef Ansi
   (CompoundKey K, Pointer V)
#else
   (K, V) CompoundKey K; Pointer V;
#endif

{
   register CompoundHeapNode *N;

   if (NodePool == Nil)
   {
      register CompoundHeapNode *Block, *P;

      /*
       * Allocate a block of nodes
       */
      Block = (CompoundHeapNode *) Allocate(sizeof(CompoundHeapNode) * NodeBlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(CreateCompoundHeapNode) Memory allocation failed.\n");
         Halt();
      }

      /*
       * Place the nodes in the block into the pool
       */
      for (P = Block; P - Block < NodeBlockSize; P++)
         FreeCompoundHeapNode(P);
   }
   
   NewCompoundHeapNode(N);
   N->Key = K;
   N->Value = V;
   N->Right = N->Left = N->Up = Nil;
   
   return N;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Insertion and deletion
 *
 */


#define IsLeftChild(N) ((N)->Up->Left == (N))
#define IsLessThan(A, B) (CompoundCompare(A, B) < 0)


/*
 * HeapInsert -- Insert a key-value pair into a heap
 *
 */
CompoundHeapNode *CompoundHeapInsert

#ifdef Ansi
   (CompoundKey K, Pointer V, register CompoundHeap *H)
#else
      (K, V, H) CompoundKey K; Pointer V; register CompoundHeap *H;
#endif

{
   register CompoundHeapNode *N;
   
   N = CreateCompoundHeapNode(K, V);
   H->Root = LinkCompoundHeap(H->Root, N);
   return N;
}


/*
 * HeapDelete -- Delete a node from a heap
 */
Void CompoundHeapDelete

#ifdef Ansi
   (register CompoundHeapNode *N, register CompoundHeap *H)
#else
   (N, H) register CompoundHeapNode *N; register CompoundHeap *H;
#endif

{
   if (N != H->Root)
   {
      if (IsLeftChild(N))
         N->Up->Left = N->Right;
      else
         N->Up->Right = N->Right;
      if (N->Right != Nil)
         N->Right->Up = N->Up;
      H->Root = LinkCompoundHeap(H->Root, DeleteCompoundHeap(N));
   }
   else
      H->Root = DeleteCompoundHeap(H->Root);
}


/*
 * Link -- Link together two heap-ordered trees
 *
 */
static CompoundHeapNode *LinkCompoundHeap

#ifdef Ansi
   (register CompoundHeapNode *S, register CompoundHeapNode *T)
#else
   (S, T) register CompoundHeapNode *S, *T;
#endif

{
   if (S == Nil)
      return T;
   else if (T == Nil)
      return S;
   else if (IsLessThan(S->Key, T->Key))
   {
      T->Up = S;
      T->Right = S->Left;
      if (T->Right != Nil)
         T->Right->Up = T;
      S->Left = T;
      return S;
   }
   else
   {
      S->Up = T;
      S->Right = T->Left;
      if (S->Right != Nil)
         S->Right->Up = S;
      T->Left = S;
      return T;
   }
}


/*
 * Delete -- Delete the root of a heap-ordered subtree
 */
static CompoundHeapNode *DeleteCompoundHeap

#ifdef Ansi
      (register CompoundHeapNode *N)
#else
   (N) register CompoundHeapNode *N;
#endif

{
   register CompoundHeapNode *Pairs, *Tree, *P, *Q, *R;

   /*
    * Delete the root of the subtree
    */
   R = N->Left;
   FreeCompoundHeapNode(N);
   if (R == Nil)
      return Nil;
   
   /*
    * Front-to-back, link the children of the root, in pairs
    */
   Pairs = Nil;
   while (R != Nil && R->Right != Nil)
   {
      P = R;
      Q = R->Right;
      R = Q->Right;
      P->Right = Q->Right = Nil;
      P = LinkCompoundHeap(P, Q);
      P->Right = Pairs;
      Pairs = P;
   }
   if (R != Nil)
   {
      R->Right = Pairs;
      Pairs = R;
   }

   /*
    * Back-to-front, link the paired trees, cumulatively
    */
   Tree = Pairs;
   R = Tree->Right;
   Tree->Right = Nil;
   while (R != Nil)
   {
      P = R;
      R = R->Right;
      P->Right = Nil;
      Tree = LinkCompoundHeap(Tree, P);
   }
   Tree->Up = Nil;
   
   /*
    * Return the resulting tree
    */
   return Tree;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Arbitrary access
 *
 */


/*
 * HeapNodeKey -- Key of a heap node
 *
 */
CompoundKey CompoundHeapNodeKey

#ifdef Ansi
   (register CompoundHeapNode *N)
#else
   (N) register CompoundHeapNode *N;
#endif

{
   return N->Key;
}


/*
 * HeapNodeValue -- Value of a heap node
 *
 */
Pointer CompoundHeapNodeValue

#ifdef Ansi
   (register CompoundHeapNode *N)
#else
   (N) register CompoundHeapNode *N;
#endif

{
   return N->Value;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Finding the minimum
 *
 */


/*
 * HeapMinKey -- Return the minimum key in a heap
 *
 */
CompoundKey CompoundHeapMinKey

#ifdef Ansi
   (register CompoundHeap *H)
#else
   (H) register CompoundHeap *H;
#endif

{
/*modified by August Woerner. was return Nil, but that doesn't work for compound values*/
    CompoundKey retValue = {-1, -1.0};
    return (H->Root != Nil) ? H->Root->Key : retValue;
}


/*
 * HeapMinValue -- Return the value associated with the minimum key in a heap
 *
 */
Pointer CompoundHeapMinValue

#ifdef Ansi
   (register CompoundHeap *H)
#else
   (H) register CompoundHeap *H;
#endif

{
   return (H->Root != Nil) ? H->Root->Value : Nil;
}


/*
 * HeapDeleteMin -- Delete the value with minimum key from a heap
 *
 */
Pointer CompoundHeapDeleteMin

#ifdef Ansi
   (register CompoundHeap *H)
#else
   (H) register CompoundHeap *H;
#endif

{
   Pointer V;
   
   if (H->Root != Nil)
   {
      V = H->Root->Value;
      H->Root = DeleteCompoundHeap(H->Root);
      return V;
   }
   else
      return Nil;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Miscellaneous
 *
 */


/*
 * HeapUnion -- Union two disjoint heaps destructively
 *
 */
CompoundHeap *CompoundHeapUnion

#ifdef Ansi
   (register CompoundHeap *A, register CompoundHeap *B)
#else
   (A, B) register CompoundHeap *A, *B;
#endif

{
   A->Root = LinkCompoundHeap(A->Root, B->Root);
   B->Root = Nil;
   DestroyCompoundHeap(B);
   
   return A;
}


/*
 * HeapDecreaseKey -- Decrease the key of a node in a heap
 *
 */
Void CompoundHeapDecreaseKey

#ifdef Ansi
   (register CompoundHeapNode *N, CompoundKey K, register CompoundHeap *H)
#else
   (N, K, H) register CompoundHeapNode *N; CompoundKey K; register CompoundHeap *H;
#endif

{
   N->Key = K;
   if (N != H->Root)
   {
      if (IsLeftChild(N))
         N->Up->Left = N->Right;
      else
         N->Up->Right = N->Right;
      if (N->Right != Nil)
         N->Right->Up = N->Up;
      N->Up = Nil;
      H->Root = LinkCompoundHeap(H->Root, N);
   }
}


/*
 * HeapIsEmpty -- Is the heap empty?
 *
 */
short CompoundHeapIsEmpty

#ifdef Ansi
   (register CompoundHeap *H)
#else
   (H) register CompoundHeap *H;
#endif

{
   return H->Root == Nil;
}
