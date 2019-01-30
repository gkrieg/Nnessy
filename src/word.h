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
/****************************************************************************
 * word.c/h
 *
 * These files define a number of procedures for working with elements
 * of the BLOSUM62-based metric space. In particular, it provides
 * procedures for reading/writing these elements in from a text/binary
 * file.
 ****************************************************************************/

#ifndef WORD_H_
#define WORD_H_

#include <stdio.h>

#include "portable.h"
/* #include "ptree.h" */

#define MAX_WORD_LENGTH 50

typedef char word[ MAX_WORD_LENGTH + 1];

void printepsilon();
void setEpsilonRho(float e, float aar, float ssr);
void setShift(float s);
void WordToString( const Pointer data, char* buffer );
float WordDist( const Pointer pp, const Pointer qq );
float WeightedWordDist( const Pointer pp, const Pointer qq );
float PrivateWordDist( const char* pp, const char* qq );
float CombinedWordDist(const Pointer p, const Pointer q);
float CombinedAlphaDistance(const Pointer p, const Pointer q);
float CombinedBetaDistance(const Pointer p, const Pointer q);
void readinweights(int lengthOfWord, const char* filename);

/* distance metric used to compute probabilities for protein secondary structure prediction */
float wordDist2Alpha(const Pointer p, const Pointer q);
float wordDist2Beta(const Pointer p, const Pointer q);
float wordDist2Loop(const Pointer p, const Pointer q);
float wordDist2Coil(const Pointer p, const Pointer q);
float WeightedAlphaDistance(const Pointer p, const Pointer q);
float WeightedBetaDistance(const Pointer p, const Pointer q);
float WeightedCoilDistance(const Pointer p, const Pointer q);
float PrivatewordDist2Alpha(const char* p, const char* q);
float PrivatewordDist2Beta(const char* p, const char* q);
float PrivateSecondaryStructureDistance(const char* w,const char* v,float ssoffset);
void setdfs(char*,char*);
int setcoildf(char*);

/* These procedures are used to read unit-test data from a
 * text-file. (That is, a SMALL number of words.) */
Pointer WordReader( FILE* datastream, int numPts );

void
adjustWeights(int wordLength, float dropoff);

/* These features are frozen while the PartitionTree module is under
 * development.
 *
 *  These procedures are used to serialize large amounts of data to a
 *  text file.

void WordWriteDataBinary( FILE* bitstream, word* words, int numPts );
word* WordReadDataBinary( FILE* bitstream, int* numwords );

PartitionTree* CreatePartitionTreeOfWords( char** points, int numPoints,
					   word** word_mem );

PartitionTree* LoadPartitionTreeOfWords( char* filename, word** word_mem, 
					 int* numWords );

void SavePartitionTreeOfWords( char* filename, word* word_mem, 
			       int num_words, PartitionTree* pt );

void DestroyWords( word* words );
*/

#endif /* WORD_H_ */
