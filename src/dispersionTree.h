#ifndef DISPERSION_TREE_
#define DISPERSION_TREE_



#include <stdio.h>
#include <limits.h>

#include "portable.h"

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

// this is used in: writeDispersionTreeToHeaderFile routine
// specifically, I want to (losslessly) print out the radius of a subtree to file
// note that if the type of value used in the distance function changes (eg, to an int), CHANGE THIS
// if you want to be able to write out a header file...
#define RADIUS_PRINTF "%.9g"

// this keeps track of the number of times the distance function was invoked... (in a non-thread safe way!)
extern unsigned long numDistanceComputations;

// this is the prototype for what the callback for the distance function looks like.
typedef distancevalue (*Distance)(const Pointer, const Pointer);

// if the type unsigned is too small, you'll need to change all instances of the word unsigned to the 
// appropriate sized unsigned integer (think typedef), and NULL_POINT
#define NULL_POINT UINT_MAX
#define NULL_DISTANCE -1


// macros to get a node in the tree, and compute the distance 
#define GET_NODE(i, t) (t->root[i])
#define GET_NODE_REF(i, t) (&(t->root[i]))

// used SPECIFICALLY to compute distances between points used to construct the tree
#define DISTANCE(i,j,t) (t->dist(t->data[i], t->data[j]))

// used SPECIFICALLY to compute the distance between a query object q
// and a particular center (c) in the tree
#define QUERYDISTANCE(q,c,d,t) ( d(q, t->data[c]))

// gets the underlying point (used to compute distances)
// used to represent a center of a (sub)tree
// 
// takes in the tree (t) and a node in the tree (n)
#define GET_CENTER(n, t) ( t->data[n])

// node in the tree
typedef struct {
    distancevalue radius; // the max distance to any point below this point in the tree
    unsigned center; // the index in *data of the center of this tree
    unsigned left; // index in *root of the left child
    unsigned right; // ditto for the right child
} DispersionTreeNode;

// the tree itself
typedef struct {
    DispersionTreeNode *root; //array-based representation of the nodes in the tree. root[0] is the root node itself.
    unsigned treeSize; // number of points in the tree
    unsigned arraySize; // the capacity of the array used to represent the tree
    Pointer *data; // the original array of data on which distance computations are performed
    Distance dist; // and of course the distance function
} DispersionTree;


// this is the only non-thread safe routine.
// it merges two subtrees, creating the parent, and returns in the index in tree->root of said parent
// if we use threads to speed up construction, you'll need a mutex around this
// (in search I do increment a distance counter in a non-thread-safe way, but that doesn't influence correctness)
unsigned
mergeNodes(DispersionTree *tree, unsigned leftIndex, unsigned rightIndex, distancevalue radius, unsigned *nextFreeIndex);

DispersionTree*
allocDispersionTree(unsigned numPoints);

DispersionTree*
createDispersionTree(Pointer *data, Distance dist, unsigned numPoints, unsigned numLevels, unsigned flags);

void
destroyDispersionTree(DispersionTree *tree);


// binary IO
// this reads a dispersion tree from a file written by the write routine
DispersionTree*
readDispersionTreeFromFile(Pointer *data, Distance dist, unsigned numPoints, char *filename);

// this writes a dispersion tree to memory (trivial to do as it's implemented as an array of structs)
// note that the points are NOT written to file (as I have no idea what they are)
// thus to read the data structure back in, I need to know what the **data are, and they have to be the *same*
// as when the tree was created. there will be no testing of this presumption.
void
writeDispersionTreeToFile(char *filename, DispersionTree *tree);


// this writes a .h file with a dispersion tree in it
void
writeDispersionTreeToHeaderFile(char *headerFilename, char *dataStructureName, DispersionTree *tree);

void
initDispersionTreeFromHeaderFile(DispersionTree *tree, Distance dist, Pointer *data);


/*
  this does proximity search
   by which I mean either a nearest neighbor search or a range search.
   It is an implementation of the best-first lower bound search of Samet. (citeme)
*/
// Pointer is just a void*
unsigned
proximitySearch(DispersionTree *tree, Pointer query, unsigned numNeighbors, Pointer *nearestNeighbors, distancevalue max);


// wrappers to the proximitySearch function that allow you to do k-nearest neighbor searches
unsigned
knnSearch(DispersionTree *tree, Pointer query, unsigned numNeighbors, Pointer *nearestNeighbors);
// and range searches...
unsigned
rangeSearch(DispersionTree *tree, Pointer query, distancevalue range, Pointer *nearestNeighbors);


#endif // header guard

