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

#ifndef __PARAPROBE_NUMABINDING_H__
#define __PARAPROBE_NUMABINDING_H__

#include "PARAPROBE_Parallelization.h"

#include <numa.h>


//define NUMA node type, these are architecture/CPU dependent
//use lstopo hwloc console command to identify target architecture layout

//some examples given here relevant for workstations at MPIE
//MK::this is tailoring the application to a specific workstation, i.e. in general not portable
struct NUMANodeMAWS30 {
	int num_cpus;
	int numa_cpus[10*2]; //a Xeon ten (hyper-threading pair) core CPU
};

struct NUMANodeMAWS15 {
    int num_cpus;
    int numa_cpus[18*2]; //a Xeon eightteen HT pair core CPU_SET
};

//##MK::add further NUMANode type structs specific to your machine
/*
struct NUMANode##### {
	int num_cpus;
	//either or
	int numa_cpus[####*2]; //if hyperthread core pair and Intel
	int numa_cpus[####*1]; //if no hyperthreading core pair and Intel
}
*/


//pick NUMANode type you want to use
typedef struct NUMANodeMAWS15 NUMANodeType;


unsigned int my_numa_bitmask_weight(const struct bitmask *mask);
/*{
	unsigned int weight = 0;
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			weight++;
		}
	}
	return weight;
}*/

#endif
