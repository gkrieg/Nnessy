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
 * metricspaces.h/c
 *
 * In the course of testing and profiling our nearest-neighbor data
 * structure, we found it necessary to implement a number of metric
 * spaces. Although one could define each metric space on its own,
 * there are some common (tedious) operations on metric spaces for
 * which it is advantageous to have a general solution. This includes
 * tasks like finding nearest neighbors by brute force, loading a set
 * of points from a file, etc.
 *
 * To that end, we introduced the metricspaces module to perform these
 * tasks. This module provides a collection of common metric spaces
 * (for testing purposes). To add a metric space to this collection,
 * the following are needed:
 *
 * 1) A procedure which evaluates the distance between two points.
 * 2) A procedure for computing a string representation of a point.
 * 3) A procedure for reading N points from a file stream.
 *
 * Given these three procedures, one can adjust the code below
 * (basically by adding case statements) to include the desired metric
 * space in the library.
 */

#ifndef METRICSPACES_H_
#define METRICSPACES_H_

#include <stdio.h>
#include "portable.h"
#include "dispersionTree.h"

typedef enum { R2_L1, R4_L1, R6_L1, 
	       R2_L2, R4_L2, R6_L2,
	       R8_L2, R16_L2,
	       R32_L2, R64_L2,
	       R2_LI, R4_LI, R6_LI,
	       PROTEIN_WORDS, HAMMING_WORDS, LEVENSHTEIN_WORDS, FAIL } SpaceType; 

SpaceType StringToSpaceType( char* str );

//typedef float (*Distance)(const Pointer, const Pointer);

/* A "ToString" procedure takes as input a point and a buffer and
   writes a  string representation of that point to the buffer. */
typedef void (*ToString)( Pointer data, char* buffer );

typedef Pointer (*Reader)( FILE* datastream, int numPts );

typedef struct {
  /* Fields for Users: */
  Distance dist;    /* Pointer to the distance function defined on the points. */ 
  ToString toString; /* Pointer to a string-defining function defined on points. */
  Pointer* points;  /* An array of pointers to the points in memory. */
  unsigned numpoints;    /* The number of points in this space. */
  SpaceType sp_type;     /* An enum describing the type of the space. */
  
  /* Fields for Internal Use: */
  size_t point_size;        /* The size in memory of a point. */ 
  Pointer point_mem_block;  /* Pointer to the block of memory holding the points. */
} MetricSpace;

/* From the given datastream (a text file), recover a set of points
 * and their associated metric space data. In all of our files, sets
 * of points are given in the following format:
 *
 * # Points: [ The number of points to follow, say N ] [ linebreak ]
 * [ A textual presentation of N points (strings, tuples, etc.) ]
 *
 */
MetricSpace* CreateMetricSpace( char* space_type, FILE* datastream );
MetricSpace** CreateMetricSpaces( char* space_type, FILE* datastream, int num_splits);

/* Free all memory associated with the input metric space. */
void DestroyMetricSpace( MetricSpace* toDestroy );

/* Occassionally, it will be useful to create a metric space
 * consisting of a single point. The procedure below accomplishes
 * this. It assumes the "# Points:" header does not appear, since it
 * is just loading a single point. */ 
Pointer LoadMetricSpacePoint( char* space_type, FILE* datastream );

/* Free all memory associated with the input point. */
void DestroyMetricSpacePoint( Pointer p );

/* Given a MetricSpace of points and a query, compute the N nearest
 * negihbors in the space to the query. This is accomplished with a
 * straightforward linear search of the points.*/
void NaiveNearestNeighborsSearch( MetricSpace* points, Pointer query, 
				  Pointer* neighbors, int N ); 

/* Given a MetricSpace of points and a query, compute an array of
 * distances which describes the distance from the query to each of
 * the points in the space. */
void ComputeAllDistances( MetricSpace* points, Pointer query, float* dist_array );

#endif /* METRICSPACES_H_ */
