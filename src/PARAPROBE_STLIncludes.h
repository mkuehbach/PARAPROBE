/*
	Copyright Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Data structure, code design, parallel implementation:
	Markus K\"uhbach, 2017-2018

	Third-party contributions:
	Andrew Breen - sequential Matlab code snippets for reconstruction and EPOS
	Markus G\"otz et al. - HPDBScan
	Kartik Kukreja - path compressed union/find
	Lester Hedges - AABBTree

	PARAPROBE --- is an MPI/OpenMP/SIMD-parallelized tool for efficient scalable
	processing of Atom Probe Tomography data targeting back-end processing.
	
	This file is part of PARAPROBE.

	PARAPROBE is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

 	PARAPROBE is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with paraprobe.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __PARAPROBE_STLINCLUDES_H__
#define __PARAPROBE_STLINCLUDES_H__

//C++ STL
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <list>

#include <cassert>
#include <map>
#include <iterator>
#include <utility>


#include <random>
#include <set>

//#define NDEBUG
#include <assert.h>


#include <stdint.h>

//C header for querying file size
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//forward declaration for global scope
using namespace std;

//generic global functions to report state, warnings, and erros
void reporting( const int rank, const string what );
void reporting( const string what );
void complaining( const int rank, const string what );
void complaining( const string what );
void stopping( const int rank, const string what );
void stopping( const string what );


//##MK::add compile time check if target system is little endian
//this is probably not portable to Intel and should be improved to run watertight
//https://stackoverflow.com/questions/4239993/determining-endianness-at-compile-time
#if __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
	#error PARAPROBE requires LittleEndian Architecture
#endif



#endif
