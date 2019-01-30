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
#include "heap.h"


typedef HeapKey   Key;
typedef HeapValue Value;
typedef HeapNode  Node;


static Node *CreateNode
   Proto(( Key K, Value V ));

static Node *Link
   Proto(( Node *A, Node *B, int (*Compare) Proto((Key, Key)) ));
         
static Node *Delete
   Proto(( Node *N, int (*Compare) Proto((Key, Key)) ));


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Creation and destruction
 *
 */


static Heap *HeapPool = Nil; /* Pool of free heaps */
static Node *NodePool = Nil; /* Pool of free nodes */

/*
 * Heap pool maintenance
 */
#define FreeHeap(H) (((H)->Root = (Node *) HeapPool), HeapPool = (H))
#define NewHeap(H)  (((H) = HeapPool), HeapPool = (Heap *) HeapPool->Root)

/* 
 * Node pool maintenance
 */
#define FreeNode(N) (((N)->Up = NodePool), NodePool = (N))
#define NewNode(N)  (((N) = NodePool), NodePool = NodePool->Up)


/*
 * CreateHeap -- Create an empty heap
 *
 */
Heap *CreateHeap

#ifdef Ansi
   (int (*Compare) Proto((Key, Key)))
#else
   (Compare) int (*Compare) Proto((Key, Key));
#endif

{
   register Heap *H;
   
   if (HeapPool == Nil)
   {
      register Heap *Block, *P;
      
      /*
       * Allocate a block of heaps
       */
      Block = (Heap *) Allocate(sizeof(Heap) * HeapBlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(CreateHeap) Memory allocation failed.\n");
         Halt();
      }

      /*
       * Place the heaps in the block into the pool
       */
      for (P = Block; P - Block < HeapBlockSize; P++)
         FreeHeap(P);
   }
   
   NewHeap(H);
   H->Compare = Compare;
   H->Root = Nil;
   
   return H;
}


/*
 * DestroyHeap -- Destroy a heap
 *
 */
Void DestroyHeap

#ifdef Ansi
   (Heap *H)
#else
   (H) Heap *H;
#endif

{
   register Node *N, *P, *R, *L;

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
      
      FreeNode(N);
      N = P;
   }
   
   FreeHeap(H);
}


/*
 * CreateNode -- Create a new heap node
 *
 */
static Node *CreateNode

#ifdef Ansi
   (Key K, Value V)
#else
   (K, V) Key K; Value V;
#endif

{
   register Node *N;

   if (NodePool == Nil)
   {
      register Node *Block, *P;

      /*
       * Allocate a block of nodes
       */
      Block = (Node *) Allocate(sizeof(Node) * NodeBlockSize);
      if (Block == NULL)
      {
         fprintf(stderr, "(CreateHeapNode) Memory allocation failed.\n");
         Halt();
      }

      /*
       * Place the nodes in the block into the pool
       */
      for (P = Block; P - Block < NodeBlockSize; P++)
         FreeNode(P);
   }
   
   NewNode(N);
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
#define IsLessThan(A, B, F) ((*(F))(A, B) < 0)


/*
 * HeapInsert -- Insert a key-value pair into a heap
 *
 */
Node *HeapInsert

#ifdef Ansi
   (Key K, Value V, register Heap *H)
#else
   (K, V, H) Key K; Value V; register Heap *H;
#endif

{
   register Node *N;
   
   N = CreateNode(K, V);
   H->Root = Link(H->Root, N, H->Compare);
   return N;
}


/*
 * HeapDelete -- Delete a node from a heap
 */
Void HeapDelete

#ifdef Ansi
   (register Node *N, register Heap *H)
#else
   (N, H) register Node *N; register Heap *H;
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
      H->Root = Link(H->Root, Delete(N, H->Compare), H->Compare);
   }
   else
      H->Root = Delete(H->Root, H->Compare);
}


/*
 * Link -- Link together two heap-ordered trees
 *
 */
static Node *Link

#ifdef Ansi
   (register Node *S, register Node *T, int (*F) Proto((Key, Key)))
#else
   (S, T, F) register Node *S, *T; int (*F) Proto((Key, Key));
#endif

{
   if (S == Nil)
      return T;
   else if (T == Nil)
      return S;
   else if (IsLessThan(S->Key, T->Key, F))
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
static Node *Delete

#ifdef Ansi
   (register Node *N, int (*F) Proto((Key, Key)))
#else
   (N, F) register Node *N; int (*F) Proto((Key, Key));
#endif

{
   register Node *Pairs, *Tree, *P, *Q, *R;

   /*
    * Delete the root of the subtree
    */
   R = N->Left;
   FreeNode(N);
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
      P = Link(P, Q, F);
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
      Tree = Link(Tree, P, F);
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
Key HeapNodeKey

#ifdef Ansi
   (register Node *N)
#else
   (N) register Node *N;
#endif

{
   return N->Key;
}


/*
 * HeapNodeValue -- Value of a heap node
 *
 */
Value HeapNodeValue

#ifdef Ansi
   (register Node *N)
#else
   (N) register Node *N;
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
Key HeapMinKey

#ifdef Ansi
   (register Heap *H)
#else
   (H) register Heap *H;
#endif

{
   return (H->Root != Nil) ? H->Root->Key : Nil;
}


/*
 * HeapMinValue -- Return the value associated with the minimum key in a heap
 *
 */
Value HeapMinValue

#ifdef Ansi
   (register Heap *H)
#else
   (H) register Heap *H;
#endif

{
   return (H->Root != Nil) ? H->Root->Value : Nil;
}


/*
 * HeapDeleteMin -- Delete the value with minimum key from a heap
 *
 */
Value HeapDeleteMin

#ifdef Ansi
   (register Heap *H)
#else
   (H) register Heap *H;
#endif

{
   Value V;
   
   if (H->Root != Nil)
   {
      V = H->Root->Value;
      H->Root = Delete(H->Root, H->Compare);
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
Heap *HeapUnion

#ifdef Ansi
   (register Heap *A, register Heap *B)
#else
   (A, B) register Heap *A, *B;
#endif

{
   A->Root = Link(A->Root, B->Root, A->Compare);
   B->Root = Nil;
   DestroyHeap(B);
   
   return A;
}


/*
 * HeapDecreaseKey -- Decrease the key of a node in a heap
 *
 */
Void HeapDecreaseKey

#ifdef Ansi
   (register Node *N, Key K, register Heap *H)
#else
   (N, K, H) register Node *N; Key K; register Heap *H;
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
      H->Root = Link(H->Root, N, H->Compare);
   }
}


/*
 * HeapIsEmpty -- Is the heap empty?
 *
 */
short HeapIsEmpty

#ifdef Ansi
   (register Heap *H)
#else
   (H) register Heap *H;
#endif

{
   return H->Root == Nil;
}
