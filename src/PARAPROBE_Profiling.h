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
#define APT_TES				8		//tessellation
#define APT_UTL				9		//utility

#define APT_IS_UNKNOWN		-1
#define APT_IS_SEQ			0
#define APT_IS_PAR			1		//information telling that parallelism is used

#define MEMORY_NOSNAPSHOT_TAKEN		-1


struct memsnapshot
{
	size_t virtualmem;
	size_t residentmem;
	memsnapshot() : virtualmem(MEMORY_NOSNAPSHOT_TAKEN),
			residentmem(MEMORY_NOSNAPSHOT_TAKEN) {}
	memsnapshot(const size_t _vm, const size_t _rm) :
		virtualmem(_vm), residentmem(_rm) {}
};


struct plog
{
	double dt;
	double tstart;
	double tend;
	size_t virtualmem;		//virtual memory consumption in bytes
	size_t residentmem;		//resident set size in bytes, i.e. number of pages process as in real memory times system specific page size
	string what;
	unsigned short typ;		//task identifier
	unsigned short pll;		//parallelism identifier
	unsigned int i;			//running number to identify the which-th snapshot
							//used because for easier utilizability of the result
							//we sort in ascending processing time thereby however having the mem data not as a time trajectory


	plog() : dt(0.0), tstart(0.0), tend(0.0), virtualmem(-1), residentmem(-1),
			what(""), typ(APT_XX), pll(APT_IS_SEQ), i(0) {}
	plog(const double _dt, const size_t _vm, const size_t _rm, const string _s,
			const unsigned short _t, const unsigned short _p, const unsigned int _i) :
				dt(_dt), tstart(0.0), tend(0.0), virtualmem(_vm), residentmem(_rm),
				what(_s), typ(_t), pll(_p), i(_i) {}
	plog(const double _ts, const double _te, const string _s,
			const unsigned short _t, const unsigned short _p, const unsigned int _i) :
				dt(_te - _ts), tstart(_ts), tend(_te), virtualmem(-1), residentmem(-1),
				what(_s), typ(_t), pll(_p), i(_i) {}	//version not tracking memory consumption
	plog(const double _ts, const double _te, const size_t _vm, const size_t _rm,
			const string _s, const unsigned short _t, const unsigned short _p,
			const unsigned int _i) :
				dt(_te - _ts), tstart(_ts), tend(_te), virtualmem(_vm), residentmem(_rm),
				what(_s), typ(_t), pll(_p), i(_i) {}	//version tracking memory consumption
};

class profiler
{
public:
	profiler() {};
	~profiler() {};

	//void prof(const string whichenv, const unsigned short category, const double st, const double en);
	void prof_elpsdtime_and_mem(const string whichenv, const unsigned short category,
			const unsigned short parallelism, memsnapshot const & mem, const double st, const double en);
	memsnapshot get_memoryconsumption( void );
	size_t get_nentries( void );
	void report_memory( pair<size_t,size_t> const & in );
	void spit_profiling( const unsigned int simid, const int rank );

private:
	vector <plog> evn;
};


#endif
