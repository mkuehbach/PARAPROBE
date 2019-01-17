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

#ifndef __PARAPROBE_H_SPACEBUCKETING_H__
#define __PARAPROBE_H_SPACEBUCKETING_H__

#include "PARAPROBE_AABBTree.h"

//implements a self-threadlocal-memory-allocating spatial index structure
//for amortized constant time queries of Moore neighborhood about spherical ball
//the container knows how the buckets have been mapped to threadlocal memory
//thus upon querying the struct sequentially by the callee using the range_rball function
//the callee will just read the data from shared memory and potentially has to query
//addresses in buckets which reference pieces of threadlocal memory
//if instead the callee is running threadparallel the individual threads will query
//the data structure and can make use of reading more likely local pieces of memory
//which will be faster at higher thread count on ccNUMA systems specifically on
//systems where multi-level hierarchies on top of the usual L1,L2,L3 cache hierarchy have
//were built

struct erase_log
{
	size_t ncleared;
	size_t nkept;
	erase_log() : ncleared(0L), nkept(0L) {}
	erase_log(const size_t _nc, const size_t _nk) :
		ncleared(_nc), nkept(_nk) {}
};


struct spacebucket
{
	sqb mdat;
	vector<vector<pos>*> buckets;

	spacebucket();
	~spacebucket();

	void initcontainer( const aabb3d & container );
	void freecontainer( void );
	void add_atom( const pos p );
	void range_rball_noclear_nosort( const pos p, apt_xyz r, vector<pos> & candidates );
	erase_log erase_rball( const p3d p, apt_xyz r );

	void write_occupancy_raw();
	size_t get_memory_consumption();
};

#endif
