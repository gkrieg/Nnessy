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

#include "metricspaces.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "word.h"
#include "r32_point.h"


static const char* SpaceTypeToString[ 20 ] = 
  { "R2_L1", "R4_L1", "R6_L1", 
    "R2_L2", "R4_L2", "R6_L2",
    "R8_L2", "R16_L2",
    "R32_L2", "R64_L2",
    "R2_LI", "R4_LI", "R6_LI",
    "PROTEIN_WORDS", "HAMMING_WORDS", "LEVENSHTEIN_WORDS" };


SpaceType StringToSpaceType( char* str ){
  int ii;
  for( ii = 0; ii < 16; ii++ ){
    if( strcmp( str, SpaceTypeToString[ ii ] ) == 0 ){
      return (SpaceType) ii;
    }
  }

  return FAIL;
}

MetricSpace* CreateMetricSpace( char* space_type, FILE* datastream ){
   
  MetricSpace* ms = malloc( sizeof( MetricSpace ) );
  if( ms == NULL ){
    fprintf( stderr, "Error: Unable to allocate memory for metric space. Aborting.\n" );
    Halt();
  }

  ms->sp_type = StringToSpaceType( space_type );

  Reader rr;
  switch( ms->sp_type ){
  case FAIL:
    fprintf( stderr, "Error: Invalid metric space type. Aborting.\n" );
    Halt();

  case R32_L2: 
    ms->dist = R32PointL2Dist;
    rr = R32PointReader;
    ms->point_size = sizeof( R32Point );
    ms->toString = R32PointToString;
    break;

  case PROTEIN_WORDS:

    rr = WordReader;
    ms->dist = WordDist;
    ms->point_size = sizeof( word );
    printf("word size is: %zu\n",ms->point_size);

    ms->toString = WordToString;
    break;

  default:
      fprintf(stderr, "Can't handle space %s\n", space_type);
      Halt();
  }   

  ms->numpoints = 0;
  fscanf( datastream, "# Points: %u\n", &( ms->numpoints ) );
  if( ms->numpoints == 0 ){
      fprintf( stderr, "Error: Illegal or absent number of points: %d.\n", ms->numpoints );
      Halt();
  }

  ms->point_mem_block = rr( datastream, ms->numpoints );


  ms->points = malloc( sizeof( Pointer ) * ms->numpoints );
  if( ms->points == NULL ){
    fprintf( stderr, "Error! Unable to allocate memory for metric space. Aborting.\n");
    Halt();
  }

  int ii;
  for( ii = 0; ii < ms->numpoints; ii++ ){
    ms->points[ ii ] = ms->point_mem_block + ii * ms->point_size;
    //printf("word after: %s\n",(char*) ms->points[ii]);
  }

  return ms;
}
MetricSpace** CreateMetricSpaces( char* space_type, FILE* datastream, int num_splits ){
  MetricSpace** spaces = malloc(num_splits * sizeof(MetricSpace*));
  int numpoints;
  fscanf( datastream, "# Points: %u\n", &( numpoints ) );
  int pointsperMS = numpoints / num_splits;
  int leftover = numpoints % num_splits;
  if( pointsperMS == 0 ){
      fprintf( stderr, "Error: Illegal or absent number of points: %d.\n", pointsperMS );
      Halt();
  }
   
  //MetricSpace* ms = malloc( sizeof( MetricSpace ) );
  if( spaces == NULL ){
    fprintf( stderr, "Error: Unable to allocate memory for metric spaces. Aborting.\n" );
    Halt();
  }
  int c = 0;
  for(c = 0;c < num_splits;c++){ 
      MetricSpace* ms = malloc(sizeof(MetricSpace));
      ms->sp_type = StringToSpaceType( space_type );

      Reader rr;
      switch( ms->sp_type ){
      case FAIL:
        fprintf( stderr, "Error: Invalid metric space type. Aborting.\n" );
        Halt();

      case R32_L2: 
        ms->dist = R32PointL2Dist;
        rr = R32PointReader;
        ms->point_size = sizeof( R32Point );
        ms->toString = R32PointToString;
        break;

      case PROTEIN_WORDS:

        rr = WordReader;
        ms->dist = WordDist;
        ms->point_size = sizeof( word );
        printf("word size is: %zu\n",ms->point_size);

        ms->toString = WordToString;
        break;

      default:
          fprintf(stderr, "Can't handle space %s\n", space_type);
          Halt();
      }   

      ms->numpoints = pointsperMS;
      if (c == 0)
          ms->numpoints += leftover;

      ms->point_mem_block = rr( datastream, ms->numpoints );


      ms->points = malloc( sizeof( Pointer ) * ms->numpoints );
      if( ms->points == NULL ){
        fprintf( stderr, "Error! Unable to allocate memory for metric space. Aborting.\n");
        Halt();
      }

      int ii;
      for( ii = 0; ii < ms->numpoints; ii++ ){
        ms->points[ ii ] = ms->point_mem_block + ii * ms->point_size;
        //printf("word after: %s\n",(char*) ms->points[ii]);
      }
      spaces[c] = ms;
  }

  return spaces;
}

void DestroyMetricSpace( MetricSpace* ms ){
  free( ms->points );
  free( ms->point_mem_block );
  free( ms );
}

Pointer LoadMetricSpacePoint( char* space_type, FILE* datastream ){
  SpaceType st = StringToSpaceType( space_type );

  Pointer retval;
  switch( st ){
  case FAIL:
    fprintf( stderr, "Error: Invalid metric space type. Aborting.\n" );
    Halt();

  case R32_L2:  
     retval = R32PointReader( datastream, 1 );
     break;

  case PROTEIN_WORDS: 
    retval = WordReader( datastream, 1 );
    break;

      default:
	  fprintf(stderr, "Can't load space...\n");
	  Halt();

  }

  return retval;
}

void DestroyMetricSpacePoint( Pointer p ){
  free( p );
}
