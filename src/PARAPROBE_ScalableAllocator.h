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

#ifndef __PARAPROBE_SCALABLEALLOCATOR_H__
#define __PARAPROBE_SCALABLEALLOCATOR_H__

#include "PARAPROBE_NUMABinding.h"

//define if to use a scalable allocator
//MK::whats behind it? reduce fragmentation and improve memory placement,
//makes sense only if threads are pinned!

//#define USE_SCALABLE_ALLOCATOR

#ifdef USE_SCALABLE_ALLOCATOR
	//jemalloc
	#include <stdlib.h>
	#include <jemalloc/jemalloc.h>

	//overload new and delete operators
	void * operator new(size_t size);
	void operator delete(void * p);

#endif

#endif
