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

#ifndef R32_POINT_H_
#define R32_POINT_H_

#include "portable.h"

#include <stdio.h>

typedef float R32Point[32];

float R32PointL1Dist( const Pointer pp, const Pointer qq );
float R32PointL2Dist( const Pointer pp, const Pointer qq );
float R32PointLIDist( const Pointer pp, const Pointer qq );

void R32PointToString( Pointer data, char* buffer );
Pointer R32PointReader( FILE* datastream, int numPts );

#endif /* R32_POINT_H_ */
