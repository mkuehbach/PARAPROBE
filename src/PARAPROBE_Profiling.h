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


#ifndef __PARAPROBE_PROFILING_H__
#define __PARAPROBE_PROFILING_H__

#include "PARAPROBE_MPIDatatypes.h"


//program profiling should use double precision in general as
//MPI_Wtime() and omp_get_wtime() fires in double precision

//type of computational operations
#define APT_XX				0		//default, unspecified
#define APT_IO				1		//I/O
#define APT_RRR				2		//ranging
#define APT_REC				3		//reconstruction
#define APT_GEO				4		//computational geometry tip surface
#define APT_BVH				5		//spatial indexing of ion positions
#define APT_PPP				6		//descriptive spatial statistics
#define APT_CLU				7		//clustering
#define APT_UTL				8		//utility

class plog
{
public:
	plog() : dt(0.0), tstart(0.0), tend(0.0), what(""), typ(APT_XX) {}
	plog(const double _dt, const string s, const unsigned short t) :
		dt(_dt), tstart(0.0), tend(0.0), what(s), typ(t) {}
	plog(const double _ts, const double _te, const string s, const unsigned short t)
		: tstart(_ts), tend(_te), what(s), typ(t) {
		dt = _te - _ts;
	}
	~plog(){}

	double get_dt(){
		return dt;
	}
	double get_tstart(){
		return tstart;
	}
	double get_tend() {
		return tend;
	}
	string get_what() {
		return what;
	}
	unsigned short get_typ() {
		return typ;
	}

private:
	double dt;
	double tstart;
	double tend;
	string what;
	unsigned short typ;
};



class profiler
{
public:
	profiler() {};
	~profiler() {};

	void prof(const string whichenv, const unsigned short category, const double st, const double en);
	size_t get_nentries( void );
	void spit_profiling( const unsigned int simid, const int rank );

private:
	vector <plog> evn;
};


#endif
