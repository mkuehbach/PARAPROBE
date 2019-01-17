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

#ifndef __PARAPROBE_PARALLELIZATION_H__
#define __PARAPROBE_PARALLELIZATION_H__

#include "PARAPROBE_Abinitio.h"

#include <mpi.h>
#include <omp.h>

#define MASTER									0
#define SINGLETHREADED							1
#define	SINGLEPROCESS							1

#define	MPI_COMM_WORLD_OMP_GET_NUM_THREADS		1 //12




//##MK::currently on the RWTH Aachen University cluster include HDF5 library
//#include "hdf5.h"

/*
//##MK::on the MAWS machine potential modification still necessary depending on where the HDF5 file is located user modification 
#include "/usr/local/hdf5/include/hdf5.h"
*/


//implicitly performance affecting choices

//file read ahead system related settings
#define SEQIO_READ_CACHE						((10)*(1024)*(1024)) //bytes
#define MPIIO_READ_CACHE						((10)*(1024)*(1024)) //bytes

#ifdef EMPLOY_SINGLEPRECISION
	#define SIMDREGISTER_WIDTH					(8) //elements assuming eight 32 bit floats to fit in 256bit wide SIMD register of contemporary processor
#else
	#define SIMDREGISTER_WIDTH					(4) //256bit can take four 64bit double at a time
#endif

#define MPIIO_EPOS_SIZE							((11)*(4)) //bytes per element
//MK::to enforce SIMD-friendly data chunking organize rawdata bucket size such that
//ideally pinned OpenMP threads can machine off buckets of integer multiple the simd width
//##MK>should be optimized for memory page size...


#define SPACEBUCKETING_BINWIDTH					(2.f) //nm

#endif
