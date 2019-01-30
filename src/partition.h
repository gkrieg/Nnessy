
#ifndef PARTITION_H_
#define PARTITION_H_

#include "dispersionTree.h"

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
   Written by August Woerner
   This does everything needed for finding good centers
   in a metric space, and then finding a good partition
   with respect to these centers.
   note that with Dispersion Trees the
   centers are solely
   used to define what a "good" partition might look like.
   (usually by assigning points to their closest center)
*/

/* 
   for the center heuristic, at what point do we break down and do the quadratic # of comparisons
   this is the variable s in the paper
*/
#ifndef CENTER_HEURISTIC_STOP
#define CENTER_HEURISTIC_STOP 15
#endif

/* likewise, this is the number of partitions made in the Center heuristic (variable r in the paper) */
#ifndef CENTER_HEURISTIC_NUM_PARTITIONS
#define CENTER_HEURISTIC_NUM_PARTITIONS 3
#endif

/* 
   this is the constant c in the partition function.
   we form a "pseudo" hyperplane, by allowing the partition size to grow within a factor of c
   of the smallest value it can be (the value obtained if the points were exactly evenly distributed)
*/
#ifndef PARTITION_CONSTANT
#define PARTITION_CONSTANT 2
#endif


typedef struct {
    unsigned pointIndex;
    unsigned centerIndex;
    distancevalue radius;
} Sorttrio;


/* O(n*m) computation that finds the point in the sample[] 
   that minimizes the sum of distances to all other points[] */
distancevalue
getOneMedoidExactly(DispersionTreeNode sample[], unsigned sampleSize, DispersionTreeNode points[], unsigned pointCount, 
		    DispersionTree *tree, DispersionTreeNode *medoid);

/* O(n log n) heuristic to find a min sum (medoid) center. */
distancevalue
Center(DispersionTreeNode points[], unsigned pointCount, DispersionTree *tree, DispersionTreeNode *best);

void
partition(DispersionTreeNode points[], unsigned pointCount, unsigned partitionSize, unsigned numPartitions, unsigned *partitionSizes, unsigned flags, DispersionTree *tree);

/* 
   This is part of the Partition routine described in the paper. It takes in two disjoint sets.
   A set of (candidate) centers (sample), and a set of points (points[]), and it attempts
   to find a set of centers that are mutually far apart (by maximizing the average distance 
   between the points in sample[] */

unsigned
findSparseSampleSwap(DispersionTreeNode points[], unsigned pointCount,
		     DispersionTreeNode sample[], unsigned sampleSize, DispersionTree *tree );

void
inPlaceCountingSort(DispersionTreeNode points[], unsigned pointCount,
                    unsigned closestTo[],
                    unsigned partitionSize[], unsigned numPartitions);

#endif // header guard
