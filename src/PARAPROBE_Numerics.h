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


#ifndef __PARAPROBE_NUMERICS_H__
#define __PARAPROBE_NUMERICS_H__

#include "PARAPROBE_Information.h"


//precision
//MK::we utiĺize f32 single precision by default
//for APT data the IVAS output is currently single precision making it the recommended and also memory leaner choice
//single vs double precision is usually faster due to indirect effect: less memory content has to be loaded thereby
//reducing memory latency costs better cache line utilization
#define EMPLOY_SINGLEPRECISION

#define EPSILON							(1.0e-6)
#define DOUBLE_EPSILON					(1.0e-12)
#define AABBINCLUSION_EPSILON			(1.0e-4)
typedef float apt_xyz;
typedef apt_xyz apt_real;
typedef size_t apt_int;
typedef double real_m33;


//type range
#define UCHARMX							(numeric_limits<unsigned char>::max())
#define UCHARMI							(numeric_limits<unsigned char>::lowest())
#define UINT64MX						(numeric_limits<size_t>::max())
#define UINT64MI						(numeric_limits<size_t>::lowest())
#define UINT32MX						(numeric_limits<unsigned int>::max())
#define UINT32MI						(numeric_limits<unsigned int>::lowest())
#define F32MX							(numeric_limits<apt_xyz>::max())
#define F32MI							(numeric_limits<apt_xyz>::lowest()) //MK::for floating point values numlimits::min is not ::lowest!
#define F64MX							(numeric_limits<double>::max())
#define F64MI							(numeric_limits<double>::lowest())
#define SIZETMX							(numeric_limits<size_t>::max())

#define RMIN							(F32MI)
#define RMAX							(F32MX)

//user-defined accuracy limits
#define MINIMUM_PDF_RESOLUTION			(0.1)		//in nanometer
#define TIPAABB_GUARDZONE				(0.1)		//in nanometer
#define MAXIMUM_NPC3D_BIN_RESOLUTION	(503)		//pixel


//handling of iontypes and flagging them for inclusion/exclusion in analyses
//##MK::unknown type needs to be zero-th type
#define UNKNOWNTYPE						0		    //must be as large to keep the largest type ##MK::(UINT32MX)
#define ION_IN_CLUSTER					1000000	    //MK::must be UINT32MX thereby we can flag ions for temporary excluding
#define ION_AT_CLUSTER					2000000		//MK::must be larger than ION_IN_CLUSTER and fit
//MK::marks are unsigned int, much larger than required to distinguish different ion types
//therefore use marks to flag them if they should be considered in the analysis or not

//MK::coordinate system is right-handed x,y,z
#define PARAPROBE_XAXIS					0
#define PARAPROBE_YAXIS					1
#define PARAPROBE_ZAXIS					2

//MK::Mersenne initialization
#define MT19937SEED						(-1)
#define MT19937WARMUP					(700000)


//##MK::Voronoi functionality
//if undefined we also output halo cells at the boundary
#define VALIDZONE_IONS_ONLY

#endif
