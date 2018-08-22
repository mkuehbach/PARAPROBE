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


#ifndef __PARAPROBE_MPIDATATYPES_H__
#define __PARAPROBE_MPIDATATYPES_H__

#include "PARAPROBE_Datatypes.h"


/*//MPI IO container for reading EPOS file

typedef struct
{
	float x;		//coordinates
	float y;
	float z;
	float m_q;		//mass-to-charge ratio
	float tof;		//time of flight
	float vdc;		//voltage DC
	float vpu;		//pulse voltage
	float detx;		//detector coordinates
	float dety;
	int nulls;		//nulls - is the time between each pulse
	int nat_pulse;	//###number of atoms?
} MPI_EPOS_EventIO;
*/

#endif