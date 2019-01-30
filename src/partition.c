
#include "dispersionTree.h"
#include "partition.h"
#include <stdlib.h>
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

// computes the min-sum center (*medoid).
// ie, finds the point in sample[] that minimizes the sum (average) distance to points[]
distancevalue
getOneMedoidExactly(DispersionTreeNode sample[], unsigned sampleSize, DispersionTreeNode points[], unsigned pointCount,
                    DispersionTree *tree, DispersionTreeNode *medoid) {

    unsigned i, j, minIndex=0;
    distancevalue minD = NULL_DISTANCE, d;
    if (pointCount < 1 || sampleSize < 1) {
        fprintf(stderr, "Should never happen\n%s\n%d %d\n", __FUNCTION__, pointCount, sampleSize);
        Halt();
    }

    for (i=0; i < sampleSize; ++i) {
        d = 0.0;

        for (j=0; j < pointCount; ++j) {
	    ++numDistanceComputations;
            d += DISTANCE(sample[i].center, points[j].center, tree);
        }

        if (i==0 || d < minD) {
            minD = d;
            minIndex =i;
        }
    }

    *medoid = sample[minIndex];
    return minD;

}

// this is the center heuristic. it is an n log n heuristic to find a center in a metric space
// with n elements in it.
distancevalue
Center(DispersionTreeNode points[], unsigned pointCount, DispersionTree *tree, DispersionTreeNode *best) {
    
    unsigned pSize = pointCount / CENTER_HEURISTIC_NUM_PARTITIONS;
    unsigned j, x, mod = pointCount % CENTER_HEURISTIC_NUM_PARTITIONS;
    DispersionTreeNode candidates[ CENTER_HEURISTIC_NUM_PARTITIONS ];

// we're looking at a small number of points, so let's exhaustively search for the best center
    if (pointCount <= CENTER_HEURISTIC_STOP || pSize < 1) {
        if (pointCount > 0)
            return getOneMedoidExactly(points, pointCount, points, pointCount, tree, best);
        else
            return 0.0;

    } else {

        j = pSize;
        if (mod) {
            ++j;
            --mod;
        }

// form an artibrary partition, and find the best center within each partition
        DispersionTreeNode *head = points;
        for(x=0; x < CENTER_HEURISTIC_NUM_PARTITIONS; ++x) {
            Center(head, j, tree, &candidates[x]);

            head += j;
            j = pSize;
            if (mod) {
                ++j;
                --mod;
            }
        }
    }
// find the best center amongst the best centers found so far
    return getOneMedoidExactly(candidates, x, points, pointCount, tree, best);
}


unsigned
SparseSampleHelperSwap( unsigned sampleSize, DispersionTreeNode select_points[], DispersionTreeNode offer,
			distancevalue distance_sums[], DispersionTree *tree){
    unsigned ii, argmax;
    distancevalue max, EE;
    distancevalue distance_table[ sampleSize ];
    /* distance_table records the distance between offer and each                                                                         * select_point. */

    distancevalue dist_sum = 0;
    for( ii = 0; ii < sampleSize; ii++ ){
	++numDistanceComputations;
	distance_table[ ii ] = DISTANCE( offer.center, select_points[ ii ].center, tree);
	dist_sum += distance_table[ ii ];
    }

    argmax = 0;
    /* initialize the max value to the 0th entry */
    max = dist_sum - distance_table[ 0 ] - distance_sums[ 0 ];

    for( ii = 1; ii < sampleSize; ii++ ){
	EE = dist_sum - distance_table[ ii ] - distance_sums[ ii ];
	if( EE > max ){
	    argmax = ii;
	    max = EE;
	}
    }

    /* updated by August                                                                                                                
     I added the for loop to update the distance matrix w/ the new distance                                                           
     values that occur b/c of the addition of the new point (and removal of the old point)                                            
     to the sparse sample */

    if( max > 0 ){

	/* update the sum of distances for the new point */
	distance_sums[ argmax ] = dist_sum - distance_table[ argmax ];

	for (ii = 0; ii < sampleSize; ++ii) {
	    if (ii != argmax) { /*   new point's distance to i - old points distance to i */
		++numDistanceComputations;
		distance_sums[ ii ] += (distance_table[ii] - DISTANCE(select_points[ii].center, select_points[ argmax].center, tree));
	    }
	}

	return argmax;

    }
    return NULL_POINT;
}


/* this finds a sparse sample, and puts those points into the array sample[]                                                          
however, we assume that the sample[] array is populated w/ points already,                                                            
and the points are swapped out 
further, we assume that points[] and sample[] are nonoverlapping (but possibly adjacent) address spaces
*/
unsigned
findSparseSampleSwap(DispersionTreeNode points[], unsigned pointCount,
		     DispersionTreeNode sample[], unsigned sampleSize, DispersionTree *tree ){
    /* In order to speed up the calculation of our objective, we save                                                                   
       some of the sums and distances we calculate in two tables.*/
    distancevalue distance_sums[ sampleSize ];
    unsigned ii, jj;
    DispersionTreeNode tmp;

//    unsigned end = ( sampleSize < pointCount ) ? sampleSize : pointCount;
    unsigned end = sampleSize;


    /* Initialize the tables. */
    for( ii = 0; ii < end; ii++ ){
/*
        tmp = sample[ii];
        sample[ ii ] = points[ ii ];
        points[ ii ] = tmp;
*/
        distance_sums[ ii ] = 0.0;
    }


    /* If the number of input points is smaller than the sample size,                                                                          the code below isn't run. */

    for( ii = 0; ii < end; ii++ )
	for( jj = 0; jj < end; jj++ ) {
	    distance_sums[ ii ] += DISTANCE( sample[ ii ].center, sample[ jj ].center, tree );
	    ++numDistanceComputations;
	}

    for( ii = 0; ii < pointCount; ii++ ) {
	jj = SparseSampleHelperSwap( sampleSize, sample, points[ ii ], distance_sums, tree );
	if (jj != NULL_POINT) { // returns the index in the sample that should be swapped (or NULL if no swap)                                     
	    tmp = sample[jj];
	    sample[jj] = points[ii];
	    points[ii] = tmp;
	}
    }
    return end;
}


/* 
   This does a counting sort using the values in closestTo[], which is a parallel array
   to points[].
   This is used when you have points assigned to a partition number, where the values in
   closestTo[] is that assignment. this allows the partition to be done in place 
*/

void
inPlaceCountingSort(DispersionTreeNode points[], unsigned pointCount,
                    unsigned closestTo[],
                    unsigned partitionSize[], unsigned numPartitions) {

    unsigned cumulative[numPartitions], offsets[ numPartitions], i, j, k;
    DispersionTreeNode temp;

    cumulative[0] = partitionSize[0];
    offsets[0] = 0;
    for (i=1; i < numPartitions; ++i) {
        cumulative[i] = cumulative[i-1] + partitionSize[i];
        offsets[i] = offsets[i-1] + partitionSize[i-1];
    }

    j=0; /* this is the point in the sparse sample whose hyperplane we are creating */
    for (i=0; i < pointCount-1; ++i) {
/* is true when you finish partitioning 1 of the hyperplanes. skip all of the points which you know                                   
   have already been correctly partitioned */
        while (i >= cumulative[j]) {

            ++j;
            if (j == numPartitions)
                break;
            i = offsets[j];
        }

        if (j == numPartitions)
            break;

        k = closestTo[ i]; /* this point belongs to the kth hyperplane */

        /* got lucky, this point just happens to be correctly partitioned */
        if (k == j)
            continue;


        temp = points[i];


        while ( closestTo[ offsets[k]  ] == k) {
            ++offsets[k];
        }
/* and put kth point into that spot, swapping in the point that did not belong in this hyperplane                                     
   into the ith index (which may or may not be right) */
        closestTo[i] = closestTo[ offsets[k] ];
        closestTo[ offsets[k] ] = k; /* technically, we don't need this side of the swap. but it's handy for debugging */
        points[i] = points[ offsets[k] ];

        points[ offsets[k] ] = temp;

        ++offsets[k];
        --i;
    }
}


int 
compare(const void *a, const void *b) {
    Sorttrio *one = (Sorttrio*) a;
    Sorttrio *two = (Sorttrio*) b;

    if (one->radius < two->radius)
	return -1;
    else if (one->radius > two->radius)
	return 1;
    return 0;
}

/* 
   This does the Partition step in the paper
   It takes in a set of points[], the number of them, and 3 pieces of information on the partition.
   The size of the partition, the (minimum) partition size, and the modulo. 
 */
void
partition(DispersionTreeNode points[], unsigned pointCount, unsigned partitionSize, unsigned numPartitions, unsigned *partitionSizes, unsigned flags, DispersionTree *tree) {
    
    if (partitionSize >= pointCount) {
	fprintf(stderr, "Cannot happen...\n");
	Halt();
    }
    
    unsigned i, j, *closestTo, minJ;
    distancevalue d, minD;
    DispersionTreeNode *tmp;

// find a sparse sample of points. the initial values of the sparse sample are just an arbitrary (first numPartition) set of points
// in practice this works better than having the starting values be an approximate dispersed set
// we speculate that the points in the approximate dispersed set are mutually far apart, but make bad centers b/c of that (they're far from everyone)
    findSparseSampleSwap(&(points[numPartitions]), pointCount - numPartitions, points, numPartitions, tree); 
    
    // the first numPartition points are (initial) centers used to form the partition
    // the first partition step is a hyperplane partition (each point is assigned to its closest center)
    closestTo = (unsigned*) Allocate(sizeof(unsigned) * pointCount );
    if (closestTo==NULL) {
	fprintf(stderr, "Not enough memory in partition!\n");
	Halt();
    }


    for (i=0; i < numPartitions; ++i) {// centers are closest to themselves by definition
	closestTo[i] = i;
	partitionSizes[i] = 1;
//	fprintf(stderr, "%s\n", tree->data[points[i].center]);
    }

// assign the remaining points to whichever center is closest
    for ( ; i < pointCount; ++i) {
	++numDistanceComputations;
	minD = DISTANCE( points[i].center, points[0].center, tree);
	minJ = 0;
	for (j=1; j < numPartitions; ++j) {
	    ++numDistanceComputations;
	    d = DISTANCE(points[i].center, points[j].center, tree);
	    if (d < minD) {
		minD = d;
		minJ = j;
	    }
	}
	closestTo[i] = minJ;
	++partitionSizes[minJ];
    }

// forms an initial partition of the points
    inPlaceCountingSort(points, pointCount, closestTo, partitionSizes, numPartitions);

    tmp = points;
    for (i=0; i < numPartitions; ++i) {
	DispersionTreeNode center;
// within this initial partition, find a good min-sum center
	Center(tmp, partitionSizes[i], tree, &center);

// look for the center we just found
	for (j=0; j < partitionSizes[i]; ++j) {
	    if (tmp[j].center == center.center)
		break;
	}
	if (j == partitionSizes[i]) {
	    fprintf(stderr, "Cannot happen! Failed to find center with index %u\n", center.center);
	    Halt();
	}
// and put it along with the other centers at the beginning of the array!
	tmp[j] = points[i];
	points[i] = center;

	tmp += partitionSizes[i];
    }


    for (i=0; i < numPartitions; ++i) {// centers are closest to themselves by definition
	closestTo[i] = i;
	partitionSizes[i] = 1;
    }


    Sorttrio *closestCenterSort = (Sorttrio*) Allocate(sizeof(Sorttrio)*(pointCount - numPartitions));

    if (closestCenterSort == NULL) {
	fprintf(stderr, "Malloc error in partition!\n");
	Halt();
    }
    for (i=0; i < pointCount - numPartitions; ++i) {
	++numDistanceComputations;
	minD = DISTANCE( points[i+numPartitions].center, points[0].center, tree);
	minJ = 0;
	for (j=1; j < numPartitions; ++j) {
	    ++numDistanceComputations;
	    d = DISTANCE(points[i+numPartitions].center, points[j].center, tree);
	    if (d < minD) {
		minD = d;
		minJ = j;
	    }
	}
	closestCenterSort[i].pointIndex = i + numPartitions;
	closestCenterSort[i].centerIndex = minJ;
	closestCenterSort[i].radius = minD;
    }
// sort each point by its min distance to the set of centers
    qsort(closestCenterSort, i, sizeof(Sorttrio), compare);

    for (i=numPartitions; i < pointCount; ++i) {
	Sorttrio tri = closestCenterSort[i-numPartitions];
	minJ = tri.centerIndex; // only redo distance computations if the best center's
	// partition is already filled up
	if (partitionSizes[ minJ ] >= PARTITION_CONSTANT * partitionSize) {
	    minD = NULL_DISTANCE;
	    minJ = NULL_POINT;

	    for (j=0; j < numPartitions; ++j) {
		if (partitionSizes[j] < PARTITION_CONSTANT * partitionSize) {
		    ++numDistanceComputations;
		    d = DISTANCE(points[ tri.pointIndex ].center, points[j].center, tree);
		    if ( minD == NULL_DISTANCE || d < minD) {
			minD = d;
			minJ = j;
		    }
		}
	    }
	}
	// and assign points to the closest center's partition (that isn't already filled)
	++partitionSizes[minJ];
	closestTo[ tri.pointIndex ] = minJ;
    }


/*
    for ( ; i < pointCount; ++i) {
 
//   this forms a pseudo-hyperplane partition.
//   this greedily assigns each (non-center) point to whichever center it's closest to, given
//   that assignment doesn't cause the partition associated with that center to grow too large.
//   With normal hyperplane partitioning, the size of the subset created can be the same 
//   (ie, not decreasing by a constant fraction ) 
	//  size as the original superset. This avoids that problem.
	minD = NULL_DISTANCE;
	minJ = NULL_POINT;

	for (j=0; j < numPartitions; ++j) {
	    if (partitionSizes[j] < PARTITION_CONSTANT * partitionSize) {
		d = DISTANCE(points[i].center, points[j].center, tree);
		if ( minD == NULL_DISTANCE || d < minD) {
		    minD = d;
		    minJ = j;
		}
	    }
	}
	closestTo[i] = minJ;
	++partitionSizes[minJ];
    }
*/    

    Free(closestCenterSort);
    inPlaceCountingSort(points, pointCount, closestTo, partitionSizes, numPartitions);
    Free(closestTo);
}
