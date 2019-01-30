#include "dispersionTree.h"
#include "partition.h"
#include "binaryHeapTieBreaker.h"
#include "compoundHeap.h"
#include "heap.h"
#include <math.h>
#include <ctype.h>
#include <string.h>


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

unsigned long numDistanceComputations = 0;

// allocates memory for a dispersion tree
DispersionTree*
allocDispersionTree(unsigned numPoints) {

// do a single alloc to get the tree and the nodes    
    DispersionTree *tree = (DispersionTree*) Allocate( 
	sizeof(DispersionTreeNode)*(numPoints + numPoints - 1) + sizeof(DispersionTree) );
    if (tree==NULL) {
	fprintf(stderr, "Failed to allocate memory for a dispersion tree\n");
	Halt();
    }
    char *mem = (char*)tree;
    tree->root = (DispersionTreeNode*) &(mem[ sizeof(DispersionTree)]);
    tree->treeSize = numPoints; // the number of distinct points (and leaves)
    tree->arraySize = numPoints + numPoints -1; // total number of nodes

    return tree;
}


/*
  This reads a dispersion tree written to file (in binary) with the corresponding write
  function. As dispersion trees are array-based, binary file IO is pretty trivial.
  Note that the write routine doesn't write the metric space, just the organization of
  the metric space. thus it still needs both the points (*data), and the distance function
  between points. if either one changes after the tree has been written, then the metric
  space will run fine, but give you the WRONG answer.
 */
DispersionTree*
readDispersionTreeFromFile(Pointer *data, Distance dist, unsigned numPoints, char *filename) {
    FILE *ptr;
    char *mem;

    ptr = fopen(filename,"rb");
    if (ptr==NULL) {
	fprintf(stderr, "Failed to open %s for reading\n", filename);
	Halt();
    }
    DispersionTree *tree = allocDispersionTree(numPoints);
    mem = (char*)tree;
    
    fread(mem,sizeof(DispersionTree) + sizeof(DispersionTreeNode)*(numPoints + numPoints - 1),1,ptr);    
    tree->root = (DispersionTreeNode*) &(mem[ sizeof(DispersionTree)]);

    tree->data = data;
    tree->dist = dist;

    fclose(ptr);
    return tree;
}

/*
  This writes a dispersion tree (in binary) to file *filename.
  Note that this doesn't write the metric space to file, just the organization of the points
  in the metric space (see readDispersionTreeFromFile)
 */

void
writeDispersionTreeToFile(char *filename, DispersionTree *tree) {
    FILE *ptr;
    ptr = fopen(filename,"wb");
    if (ptr == NULL) {
	fprintf(stderr, "Problem opening file: %s for writing\n", filename);
	Halt();
    }

    fwrite(tree,sizeof(DispersionTree) + sizeof(DispersionTreeNode)*(tree->arraySize),1,ptr);  // the tree itself is already laid out
    // contiguously in memory...
    fclose(ptr);
}

/*
  This is for use when the dispersion tree is written as a .h file!!
  when a dispersion tree is written to file
  I can't write the distance function or the data itself to disk
  this initializes them to the correct values...
*/
void
initDispersionTreeFromHeaderFile(DispersionTree *tree, Distance dist, Pointer *data) {
    tree->dist = dist;
    tree->data = data;
}

/*
  This writes a dispersion tree to a .h file. I don't use an extern, but instead just write the structure
  straight to the .h file. In general this isn't a great practice (namespace collisions being one problem),
  but for now that's how it's implemented.
  It writes a header file (.h file)
  in the header file is 2 data structures. A DispersionTree struct, and an array of DispersionTreeNode[] structs
  these have variable names dataStructureName and dataStructureNameRoot[].
  The contents of the .h file are taken from tree
 */
void
writeDispersionTreeToHeaderFile(char *headerFilename, char *dataStructureName, DispersionTree *tree) {
    FILE *ptr;
    ptr = fopen(headerFilename,"w");
    if (ptr == NULL) {
	fprintf(stderr, "Problem opening file: %s for writing\n", headerFilename);
	Halt();
    }
    unsigned i, strLength = strlen(dataStructureName);
    char *headerName = (char*)Allocate(strLength + 20);
    if (headerName==NULL) {
	fprintf(stderr, "malloc error in writeDispersionTreeToHeaderFile\n");
	Halt();
    }
    for (i=0; i < strLength; ++i) {
	headerName[i] = toupper( (int) dataStructureName[i]);
	if (isspace((int)headerName[i])) {
	    headerName[i] = '_';
	}
    }
    headerName[i] = '_';
    headerName[i+1] = '\0';

// write the header guard and minimalistic include statements
    fprintf(ptr, "#ifndef %s\n#define %s\n\n#include \"dispersionTree.h\"\n\n", headerName, headerName);

    fprintf(ptr, "DispersionTreeNode %sroot[] = {\n", dataStructureName);
    for (i=0; i < tree->arraySize; ++i) {
	fprintf(ptr, "{"); 
	fprintf(ptr, RADIUS_PRINTF, tree->root[i].radius);
	fprintf(ptr, ", %u, %u, %u}", tree->root[i].center, tree->root[i].left, tree->root[i].right);
	if (i < tree->arraySize-1) 
	    fprintf(ptr, ",\n");
	else
	    fprintf(ptr, "\n");
    }

    
    fprintf(ptr, "};\n\nDispersionTree %s = {%sroot, %u, %u, NULL, NULL};\n", 
	  dataStructureName, dataStructureName, tree->treeSize, tree->arraySize);

// header guard stop
    fprintf(ptr, "\n\n#endif\n");


    Free(headerName);
    fclose(ptr);
}



// This merges two subtrees (left and right). The left's center is reused, and the right subtree is subsumed
// this routine is NOT threadsafe. 
// it returns the index in tree->root[] of the center point created
unsigned
mergeNodes(DispersionTree *tree, unsigned leftIndex, unsigned rightIndex, distancevalue radius, unsigned *nextFreeIndex) {
    DispersionTreeNode *next = &(tree->root[*nextFreeIndex]);
    (*next).left = leftIndex;
    (*next).right = rightIndex;
    (*next).radius = radius;
    (*next).center = tree->root[leftIndex].center;

    --(*nextFreeIndex); 
    return *nextFreeIndex+1;
}


// this is a helper function (only used to debug / check accuracy of the radii of subtrees)
// it takes in a tree, an index in tree->root[]
// and a center, and it computes the max distance from that center to
// any point in that subtree (aka, the radius of the bounding ball centered around
// tree->root[index]

distancevalue
getRadius(DispersionTree *tree, unsigned index, DispersionTreeNode center) {
    distancevalue radius;
    if (index==NULL_POINT)
	return 0;

    ++numDistanceComputations;
    radius = DISTANCE(center.center, tree->root[index].center, tree);
    distancevalue d = getRadius(tree,  tree->root[index].left, center);
    if (d > radius)
	radius = d;

    d = getRadius(tree,  tree->root[index].right, center);
    if (d > radius)
	radius = d;

    return radius;
}


// another helper-- brute force verification that the radii are correct.
void
verifyRadii(DispersionTree *tree, DispersionTreeNode *node) {
    distancevalue radius = getRadius(tree, node - tree->root, *node);
    if (radius != node->radius) {
	fprintf(stderr, "problema! %f %f\n", radius, node->radius);
	Halt();
    }
    if (node->left != NULL_POINT) {
	verifyRadii(tree, GET_NODE_REF(node->left, tree));
	verifyRadii(tree, GET_NODE_REF(node->right, tree));
    }
}

// this is the recursive implementation of building a dispersion (sub)tree
// there are two cases-- one where we merge leaves of the tree (no partitioning required)
// and one where we merge internal nodes.
// 
// it returns the index of the subtree in the array of subtree nodes kept in *tree
unsigned
buildDispersionSubtree(DispersionTree *tree, DispersionTreeNode points[], unsigned pointCount, unsigned numLevels, unsigned flags, unsigned currLevel, unsigned *nextFreeIndex) {


    unsigned numPartitions = pow(tree->treeSize, 1.0/numLevels);
    unsigned partitionSize = pointCount / numPartitions;
    unsigned i, j, k, bestCenterIndex=0;
    distancevalue *inner, **bounds, maxD, d;
    distancevalue *sumDist2Center, bestCenterRadius=NULL_DISTANCE;

    Intpair pair;
    BinaryHeapTiebreaker *heap;
    unsigned *centerIndexes, numCenters=numPartitions;
    char *exclude;

    if (pointCount == 0) {
	fprintf(stderr, "Cannot happen! no points to partition?\n");
	Halt();
    } else if (pointCount == 1) {
// this is an unlikely case-- a partition of a single point...
// return the index of that point in tree->root[]
	return points - tree->root; 
    }

    ++currLevel;
    if (currLevel < numLevels && partitionSize >= 3) { // check for a constant number of points in each partition...
	unsigned *partitionSizes = (unsigned*) Allocate(sizeof(unsigned)*numPartitions);
	centerIndexes = (unsigned*) Allocate(sizeof(unsigned)*numPartitions);
	DispersionTreeNode *tmp;
	heap = initHeapTB( numPartitions * numPartitions - numPartitions, MINHEAP);

	if (partitionSizes==NULL || centerIndexes==NULL){
	    fprintf(stderr, "Failed to malloc in buildDispersionTree\n");
	    Halt();
	}

// creates numPartition partitions, each with a size at most 2*partitionSize
// the sizes of each of the subsets created is stored in partitionSizes
	partition(points, pointCount, partitionSize, numPartitions, partitionSizes, flags, tree);
	tmp = points;
// build the subtrees
	for(i=0; i < numPartitions; ++i) {
	    centerIndexes[i] = buildDispersionSubtree(tree, tmp, partitionSizes[i], numLevels, flags, currLevel, nextFreeIndex);
	    tmp += partitionSizes[i];
	}
// fill out the array bounds[][]
// bounds [i][j] is the radius of the subtree formed by merging the jth subtree into the ith, with the center
// being the ith's
				    
	bounds = (distancevalue**) Allocate(sizeof(distancevalue*)*numPartitions);
	inner = (distancevalue*) Allocate(sizeof(distancevalue)* (numPartitions*numPartitions + numPartitions));
	if (bounds==NULL || inner==NULL) {
	    fprintf(stderr, "Malloc error in partition\n");
	    Halt();
	}
	sumDist2Center = &(inner[ numPartitions*numPartitions]);

	for (i=0; i < numPartitions; ++i, inner += numPartitions) {
	    bounds[i]=inner;
	    sumDist2Center[i] = 0;
	    DispersionTreeNode subtreeRoot = GET_NODE(centerIndexes[i], tree); // the ith subtree
	    tmp = points;
	    // and iterate through the partitions, and compute all possible radii
	    for (k=0; k < numPartitions; ++k) {

		if (i==k ) { // this computation has already been done
		    bounds[i][i] = subtreeRoot.radius;
		    tmp += partitionSizes[k];
		} else {
		    maxD = subtreeRoot.radius; // the smallest the new subtree will be is the size of the current (ith) subtree

		    // compute the radius of the tree resulting from joining the ith and the jth subtrees
		    for (j=0; j < partitionSizes[k]; ++j) {
			++numDistanceComputations;
			d = DISTANCE( (*tmp).center, subtreeRoot.center, tree);
			if (d > maxD)
			    maxD = d;

			sumDist2Center[i] += d;
			++tmp;
		    }
		    bounds[i][k] = maxD;
		}
	    }
	    tmp = points;
	}

	distancevalue radius;
	for (i=0; i < numPartitions; ++i) {
	    radius = 0;

	    for (j=0; j < numPartitions; ++j) {
		if (i != j) {

		    pair.x = i;
		    pair.y = j;
                    heapInsertTB(heap, bounds[i][j],  sumDist2Center[i],  pair);
		    if (bounds[i][j] > radius)
			radius = bounds[i][j];
		}
	    }

	    if(bestCenterRadius == NULL_DISTANCE || radius < bestCenterRadius) {
		bestCenterIndex = i;
		bestCenterRadius = radius;
	    }
	}
	Free(partitionSizes);
    } else { // case of leaves (no partitions)
	bounds = (distancevalue**) Allocate(sizeof(distancevalue*)*pointCount);
	inner = (distancevalue*) Allocate(sizeof(distancevalue)*(pointCount*pointCount + pointCount) );
	centerIndexes = (unsigned*) Allocate(sizeof(unsigned)*pointCount);

	if (bounds==NULL || inner==NULL||centerIndexes==NULL) {
	    fprintf(stderr, "Malloc error in partition\n");
	    Halt();
	}

	numCenters = pointCount;
	sumDist2Center  = &(inner[ pointCount*pointCount]);
	heap = initHeapTB( pointCount * pointCount - pointCount, MINHEAP);

	for (i=0; i < pointCount; ++i) {
	    bounds[i] = inner;
	    sumDist2Center[i] = 0;
	    for (j=0; j < pointCount; ++j, ++inner) {
		if (i != j) {
		    ++numDistanceComputations;
		    *inner = DISTANCE( points[j].center, points[i].center, tree);
		    sumDist2Center[i] += *inner;
		} else
		    bounds[i][j] = 0;
	    }
	}
	distancevalue radius;

	for (i=0; i < pointCount; ++i) {
	    radius = 0;
	    centerIndexes[i] = points - tree->root + i;
	    for (j=0; j < pointCount; ++j, ++inner) {
		if (i != j) {

		    pair.x = i;
		    pair.y = j;
                    heapInsertTB(heap, bounds[i][j],  sumDist2Center[i],  pair);

// while we're here, compute the best radius for the whole subtree
		    if (bounds[i][j] > radius) 
			radius = bounds[i][j];

		}
	    }

	    if (bestCenterRadius == NULL_DISTANCE || bestCenterRadius > radius) {
		bestCenterIndex= i;
		bestCenterRadius = radius;
	    }
	}
    }


    exclude = (char*) Allocate(numCenters);
    if (exclude == NULL) {
	fprintf(stderr, "Malloc error in partition!\n");
	Halt();
    }

    for (i=0; i < numCenters; ++i)
	exclude[i]=0;

    unsigned bestj, bestk;
    for (i=2; i < numCenters; ++i) { // merge until there are two centers left.

// find the next pair of subtrees to merge
	while (1) {
            pair = getBestItemTB(heap);
            d = getBestDistanceTB(heap);
	    heapExtractTB(heap);
	    if (! exclude[ pair.x] && ! exclude[ pair.y]) { // neither partition has been deleted                                
                if (d < bounds[ pair.x][pair.y]) { // but the bounds are incorrect (stale in heap)                                      
		    heapInsertTB(heap, bounds[ pair.x][pair.y], sumDist2Center[pair.x],  pair);
		} else if (currLevel > 1 && pair.y == bestCenterIndex) { // do not let the best center get consumed (unless this is the top level of the tree... in which case it doesn't matter if we pass up the best center

		} else { // bounds are correct                                                                                   
		    break;
		}

	    }
	}
	bestj = pair.x;
	bestk = pair.y;

        exclude[bestk]=1; // this subtree is subsumed
// merge the nodes. the bestj one remains, but now its bounds may be wrong
        centerIndexes[bestj] = mergeNodes(tree, centerIndexes[bestj], centerIndexes[bestk], d, nextFreeIndex);

// now that we've merged the two subtrees, we recompute the radii of all remaining subtree radii that may be effected
        for (j=0; j < numCenters; ++j) {
            if (! exclude[j])
                bounds[j][bestj] = MAX( bounds[j][bestj], bounds[j][bestk]);
        }
// you do NOT need to do this for the rows involving bestj, eg bounds[bestj][*], b/c the we're picking the merge w/ minimal radius
// ie, the MAX operation would do nothing...
	
    }

// merge the last two balls
    bestj = bestk = NULL_POINT;
    for (i=0; i < numCenters; ++i) {
	if (! exclude[i]) {
	    if (bestj == NULL_POINT)
		bestj = i;
	    else if (bestk == NULL_POINT)
		bestk = i;
	    else {
		fprintf(stderr, "Cannot happen. error merging last two balls...\n");
		Halt();
	    }
	}
    }
    if (bestk == NULL_POINT) {
	fprintf(stderr, "Cannot happen. error merging last two balls (not enough free?)...\n");
	for (i=0; i < numCenters; ++i)
	    printf("%u\n", (unsigned) exclude[i]);

	Halt();
    }
// ensure that the best center isn't consumed (merged and center consumed)
    if ( bestk == bestCenterIndex) {
	bestk = bestj;
	bestj = bestCenterIndex;
    }

    j = mergeNodes(tree, centerIndexes[bestj], centerIndexes[bestk], bounds[bestj][bestk], nextFreeIndex);


    destroyHeapTB(heap);
    Free(*bounds);
    Free(bounds);
    Free(centerIndexes);
    Free(exclude);

    return j;
}



/*
  This creates a dispersion tree from scratch. It takes in a metric space (the points, data[], and a distance function between them)
  it also needs to know how many points are in the metric space, the number of levels used when the tree is constructed.
  I also have built in the ability to do a littler performance tuning with the flags argument (not implemented yet).
 */
DispersionTree*
createDispersionTree(Pointer *data, Distance dist, unsigned numPoints, unsigned numLevels, unsigned flags) {

    if (numLevels == 0) {
	fprintf(stderr, "To make a dispersion tree you need a positive number of levels. I recommend 2!\n");
	Halt();
    }

    DispersionTree *tree = allocDispersionTree(numPoints);
    unsigned i;
// put the leaves in the end of the array
    DispersionTreeNode *root = tree->root + numPoints - 1;
//    printf("%p %p\n", root, &(tree->root[ numPoints-1]));
//    Halt();
    for (i=0; i < numPoints; ++i) {
	(*root).center = i; // initializes all leaf nodes to an aribtrary point with no radius and no children...
	(*root).radius = 0.0;
	(*root).left = NULL_POINT;
	(*root).right = NULL_POINT;
	++root;
    }
    tree->data = data;
    tree->dist = dist;
// and we'll fill out the inner nodes in the tree going from left to right
// note this value is accessed in the mergeNodes routine, which is NOT thread safe
    unsigned nextFreeIndex = numPoints-2;
    buildDispersionSubtree(tree, tree->root + numPoints - 1, numPoints, numLevels, flags, 0, &nextFreeIndex);
    return tree;
}


// works as such b/c of how the allocDispersionTree was alloc'd
void
destroyDispersionTree(DispersionTree *tree) {
    Free(tree);
}

/*
################################ This section is for searching a dispersion tree! ################################ 
*/

/*
// the callback used in the minheap used for the best-first search
int
compareHeapkeysMinheap(distancevalue k1, distancevalue k2) {
    distancevalue diff = k1 - k2;
    if (diff < 0)
	return -1;
    else if (diff > 0)
	return 1;
    return 0;
}
*/

// the callback used to form the set of nearest neighbors (maxheap)

int
compareHeapkeysMaxheap(distancevalue k1, distancevalue k2) {
    distancevalue diff = k1 - k2;

    if (diff > 0)
	return -1;
    else if (diff < 0)
	return 1;
    return 0;
}


/* 
   This is a helper function for proximitySearch
   It's envoked whenever a new point is encountered in the proximitysearch function.
   *maxguess is used with a range search
   K is used with a KNN search
   (both can be used, too, wherein at most K neighbors are found w/ distance at most maxguess)
   The candidates are the current k-nearest neighbors, and their max distance to the query is *knnDistance
*/
void
considerPoint(Pointer candidate, distancevalue d, unsigned *numCandidates, Heap *candidates, unsigned K, distancevalue *knnDistance, distancevalue * maxguess) {

    if (*maxguess < 0 || d <= *maxguess) {

	if (*numCandidates < K) { // just take the first few entries until we get enough neighbors
	    ++(*numCandidates);
	    HeapInsert(d, candidate, candidates);
	    // for convenience I keep the distance to the kth nearest neighbor
	    if (*numCandidates == K && *knnDistance == NULL_DISTANCE) {
		*knnDistance = HeapMinKey(candidates); 
	    }
	    
	} else if (*knnDistance > d) { // the new candidate is closer, so lets add it
	    HeapDeleteMin(candidates);
	    HeapInsert(d, candidate, candidates);
	    *knnDistance = HeapMinKey(candidates);  
	    if (*maxguess >= 0) { //for use in searches that do both range and knn searching...
		*maxguess = *knnDistance;
	    }
	}
    }
}


/*
  This is the catch-all Knn/range search routine.
  It takes in a tree, a query point, a K (set to UINT_MAX incase of a range search)
  an array to write the nearest-neighbors too
  and a maxguess value for use in a range search (set to NULL_DISTANCE, or really any negative value in a KNN search)
  It is an implementation of the range-optimal search of:
  Hjaltason, Gísli R., and Hanan Samet. "Incremental similarity search in multimedia databases." (2000).

  Range-optimal means that when the KNN search is performed, no more distance computations are used then what
  would be required had a range search been performed using the true knn-distance as the range.
  This isn't quite optimal (in the case of ties, eg multiple neighbors at the K'th nearest neighbor distance),
  but it's really close.
 */
unsigned
proximitySearch(DispersionTree *tree, Pointer query, unsigned K, Pointer *nearestNeighbors, distancevalue maxguess) {
    unsigned i, numCandidates = 0;
    distancevalue lb, parentlb, d, knnDistance = NULL_DISTANCE;
    CompoundKey heapkey;
    Distance distanceFunction = tree->dist;
    DispersionTreeNode *parent, *left, *right;
// this will be the set of nearest neighbors (in a heap)
    Heap* candidates = CreateHeap( &compareHeapkeysMaxheap );

// and this is the heap used to search the tree
    CompoundHeap* priority = CreateCompoundHeap(); // the callback was removed (to make it run more quickly)

    numDistanceComputations=1; // count the number of times the distance function was invoked (for each query)
    d = QUERYDISTANCE(query, tree->root[0].center, distanceFunction, tree);

    considerPoint( GET_CENTER(tree->root[0].center, tree), d, &numCandidates, candidates, K, &knnDistance, &maxguess);

    heapkey.lowerboundDistance = 0; 
    heapkey.centerdistance = d;
    CompoundHeapInsert(heapkey, &(tree->root[0]), priority);

    while (! CompoundHeapIsEmpty(priority) ) {
	heapkey = CompoundHeapMinKey(priority);
	parent = (DispersionTreeNode*) CompoundHeapDeleteMin(priority);
	parentlb = heapkey.centerdistance - parent->radius; 

// knndistance may have shrunk since we put parent on the heap, 
// however if we can prune this node, we can prune ALL nodes on the heap

	if (knnDistance != NULL_DISTANCE && parentlb >= knnDistance) {// the case of a k-nearest neighbor search
	    break;
	} else if (maxguess >= 0 && parentlb > maxguess) {// the case of a range search (note > vs >=)
	    break;
	}


	if (parent->left != NULL_POINT) { // check to see if this is a leaf
// the case of the left subchild!
	    left = GET_NODE_REF(parent->left, tree);
	    lb = heapkey.centerdistance - left->radius; // centers of parent and left child are the same...
	    // and lb >= parentlb (it's the same center)

	    if ( (knnDistance == NULL_DISTANCE || lb <= knnDistance) 
		 || (maxguess >=  0 && lb <= maxguess  )  ){ // might this subtree yield a closer neighbor?
		heapkey.lowerboundDistance= lb;
		CompoundHeapInsert(heapkey, left, priority);	    
	    }
	}

	if (parent->right != NULL_POINT) { // check to see if this is a leaf (separate checks necessary if/when insert/deletes get implmeneted. otherwise left==NULL_POINT implies right==NULL_POINT
// the case of the right subtree
	    right = GET_NODE_REF(parent->right, tree);
	    heapkey.centerdistance = QUERYDISTANCE(query, right->center, distanceFunction, tree); // centers are different.
	    ++numDistanceComputations;
	    lb = heapkey.centerdistance - right->radius;
	    if (parentlb > lb) // both lower-bounds are correct-- pick the max to improve pruning
		lb = parentlb;

	    if ( (knnDistance == NULL_DISTANCE || lb <= knnDistance) 
		 || (maxguess >=  0 && lb <= maxguess  )  ){ // might this subtree yield a closer neighbor?
		heapkey.lowerboundDistance= lb;
// right subtrees always involve a new point
		considerPoint( GET_CENTER(right->center, tree), heapkey.centerdistance, &numCandidates, 
			       candidates, K, &knnDistance, &maxguess);
		CompoundHeapInsert(heapkey, right, priority);
	    }
	}
    } 

    for (i = numCandidates ; i > 0; --i) {
	nearestNeighbors[i-1] = HeapDeleteMin(candidates);
    }

    DestroyHeap(candidates);
    DestroyCompoundHeap(priority);

    return numCandidates; // the number of neighbors found. may be <= the requested (b/c of a the hybrid range search option)
}


/* 
   this is a simple wrapper to the proximity search function for use with KNN searching 
   nearestNeighbors[] is where the numNeighbors closest to query are in the tree
*/
unsigned
knnSearch(DispersionTree *tree, Pointer query, unsigned numNeighbors, Pointer *nearestNeighbors) {
    return proximitySearch(tree, query, numNeighbors, nearestNeighbors, NULL_DISTANCE);
}

// and another wrapper for doing a range search. 
// TODO: more gracefully handle how the nearest-neighbors are reported
unsigned
rangeSearch(DispersionTree *tree, Pointer query, distancevalue range, Pointer *nearestNeighbors) {
    return proximitySearch(tree, query, UINT_MAX, nearestNeighbors, range);
}
