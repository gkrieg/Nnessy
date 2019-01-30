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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "r32_point.h"

float R32PointL1Dist( const Pointer pp, const Pointer qq ){
  float* p_ptr = *((R32Point*) pp);
  float* q_ptr = *((R32Point*) qq);

  float sum = 0;
  float diff;
  int ii;
  for( ii = 0; ii < 32; ii++ ){
    diff = p_ptr[ ii ] - q_ptr[ ii ];
    sum += fabs( diff );
  }

  return sum;
}

float R32PointL2Dist( const Pointer pp, const Pointer qq ){
  float* p_ptr = *((R32Point*) pp);
  float* q_ptr = *((R32Point*) qq);

  float sum = 0;
  float diff;
  int ii;
  for( ii = 0; ii < 32; ii++ ){
    diff = p_ptr[ ii ] - q_ptr[ ii ];
    sum += diff * diff;
  }

  return sqrt( sum );
}

float R32PointLIDist( const Pointer pp, const Pointer qq ){
  float* p_ptr = *((R32Point*) pp);
  float* q_ptr = *((R32Point*) qq);

  float max = -1;
  float diff;
  int ii;
  for( ii = 0; ii < 32; ii++ ){
    diff = fabs( p_ptr[ ii ] - q_ptr[ ii ] );
    if( max < diff ){
      max = diff;
    }
  }

  return max;
}

void R32PointToString( const Pointer data, char* buffer ){
  float* pdata = *( (R32Point*) data );
  sprintf( buffer, "(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f)",
	   pdata[0], pdata[1], pdata[2], pdata[3], pdata[4], pdata[5], pdata[6], pdata[7], pdata[8],
	   pdata[9], pdata[10], pdata[11], pdata[12], pdata[13], pdata[14], pdata[15],
	   pdata[16], pdata[17], pdata[18], pdata[19], 
	   pdata[20], pdata[21], pdata[22], pdata[23], pdata[24], pdata[25], pdata[26], pdata[27], pdata[28],
	   pdata[29], pdata[30], pdata[31]);
	   
}

Pointer R32PointReader( FILE* datastream, int numPts ){
  R32Point* data = malloc( sizeof(R32Point) * numPts );

  if( data == NULL ){
    fprintf( stderr, "Unable to allocate memory for data points. Aborting.\n");
    Halt();
  }

  int c, ii, j;
  for( ii = 0; ii < numPts; ii++ ){
    float* point = (float*)( data + ii );
    for (j=0; j < 32; ++j)
	point[j] = 0;
    do {
	c = getc(datastream);
    } while (c != '(' && c != EOF);
    if (c == EOF) {
	fprintf(stderr, "Premature ending!\n");
	Halt();
    }    
    /* Be cautious editing the line below, scanf parsing is tricky. */
    fscanf( datastream, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f) ", 
	    point + 0, point + 1, point + 2, point + 3, point + 4, point + 5,
	    point + 6, point + 7, point + 8, point + 9,
	    point + 10, point + 11, point + 12, point + 13, point + 14, point + 15,
	    point + 16, point + 17, point + 18, point + 19, 
	    point + 20, point + 21, point + 22, point + 23, point + 24, point + 25,
	    point + 26, point + 27, point + 28, point + 29,
	    point + 30, point + 31);

  }

  return data;
}
