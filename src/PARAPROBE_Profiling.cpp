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


#include "PARAPROBE_Profiling.h"


void profiler::prof(const string whichenv, const unsigned short category,
		const double st, const double en)
{
	evn.push_back( plog(st, en, whichenv, category) );
}

size_t profiler::get_nentries( void ){
	return evn.size();
}

bool SortProfLogAscWallClock( plog & first, plog & second )
{
	bool comp = first.get_dt() < second.get_dt();
	return comp;
}

void profiler::spit_profiling( const unsigned int simid, const int rank )
{
	//##MK::further optimization aand convenience tasks: bundle all in one file, incr ID and so forth
	//##MK::suboptimal... one file per rank
	string fn = "PARAPROBE.SimID." + to_string(simid) + ".Rank." + to_string(rank) + ".MyProfiling.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "What;Category;WallClock;CumulatedWallClock;CDF;WallClockFraction\n";
		csvlog<< ";;s;s;1;1\n";
		csvlog << "What;Category;MPI_Wtime;CumulatedWallClock;CDF;WallClockFraction\n";

		//build map of categories
		map<unsigned int, string> categories;
		categories[APT_XX] = "APT_XX";
		categories[APT_IO] = "APT_IO";
		categories[APT_RRR] = "APT_RRR";
		categories[APT_REC] = "APT_REC";
		categories[APT_GEO] = "APT_GEO";
		categories[APT_BVH] = "APT_BVH";
		categories[APT_PPP] = "APT_PPP";
		categories[APT_CLU] = "APT_CLU";
		categories[APT_UTL] = "APT_UTL";

		//sort events increasing wallclock time
		sort( evn.begin(), evn.end(), SortProfLogAscWallClock);

		//compute total time
		double dt_total = 0.f;
		for(auto it = evn.begin(); it != evn.end(); ++it) { dt_total += it->get_dt(); }

		//report
		double dt_cumsum = 0.f;
		for (auto it = evn.begin(); it != evn.end(); ++it) {
			dt_cumsum += it->get_dt();
			auto cat = categories.find(it->get_typ());
			csvlog << it->get_what() << ";" << cat->second << ";" << it->get_dt();
			csvlog << ";" << dt_cumsum << ";" << (dt_cumsum / dt_total) << ";" << (it->get_dt() / dt_total) << endl;
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local processing files");
	}
}

//program profiling should use double precision in general as MPI_Wtime() and omp_get_wtime() fires in double precision
