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


/*void profiler::prof(const string whichenv, const unsigned short category,
		const double st, const double en)
{
	evn.push_back( plog(st, en, whichenv, category) );
}*/


void profiler::prof_elpsdtime_and_mem(const string whichenv,
		const unsigned short category, const unsigned short parallelism,
		memsnapshot const & mem, const double st, const double en )
{
	evn.push_back( plog(st, en, mem.virtualmem, mem.residentmem, whichenv, category, parallelism, evn.size()) );
}


memsnapshot profiler::get_memoryconsumption( void )
{
	//from http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
	//according to D. Morovoz/T. Peterka tess2 code https://github.com/diatomic/tess2/blob/90ca0536299cdc6bbd6b2d1547d939ffbea68539/examples/memory.h
	//'file' stat seems to give the most reliable results

	//ultimately this is file system interaction, so
	//MK::MUST NOT BE CALLED FROM WITHIN THREADED REGION
	memsnapshot out = memsnapshot();

	ifstream stat_stream("/proc/self/stat", ios_base::in);
	if ( stat_stream.good() == true ) {
		//dummies for leading entries in stat we don't care about
		string pid, comm, state, ppid, pgrp, session, tty_nr;
		string tpgid, flags, minflt, cminflt, majflt, cmajflt;
		string utime, stime, cutime, cstime, priority, nice;
		string O, itrealvalue, starttime;

		//the two fields we want
		unsigned long vsize;
		long rss;

		stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
					>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
					>> utime >> stime >> cutime >> cstime >> priority >> nice
					>> O >> itrealvalue >> starttime >> vsize >> rss;  //don't care about the rest

		stat_stream.close();

		long page_size_byte = sysconf(SYSTEMSPECIFIC_POSIX_PAGESIZE);
		out.virtualmem = static_cast<size_t>(vsize);
		out.residentmem = static_cast<size_t>(rss*page_size_byte);
	}
	return out;
}


size_t profiler::get_nentries( void ){
	return evn.size();
}


bool SortProfLogAscWallClock( plog & first, plog & second )
{
	return first.dt < second.dt;
}


void profiler::report_memory( pair<size_t,size_t> const & in )
{
	cout << "VM/RSS in Bytes\t\t" << in.first << "\t\t" << in.second << endl;

	cout << "VM/RSS in GB\t\t" << BYTE2GIGABYTE(static_cast<double>(in.first)) << "\t\t" << BYTE2GIGABYTE(static_cast<double>(in.second)) << endl;
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
		csvlog << "What;ID;Category;ParallelismInfo;ProcessVirtualMemory;ProcessResidentSetSize;WallClock;CumulatedWallClock;WallClockFraction\n";
		csvlog<< ";;;;B;B;s;s;1;1\n";
		csvlog << "What;ID;Category;ParallelismInfo;ProcessVirtualMemory;ProcessResidentSetSize;WallClock;CumulatedWallClock;WallClockFraction\n";

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
		categories[APT_TES] = "APT_TES";
		categories[APT_UTL] = "APT_UTL";
		map<unsigned short, string> parallelism;
		parallelism[APT_IS_PAR] = "PARALLEL";
		parallelism[APT_IS_SEQ] = "SEQUENTIAL";


		//sort events increasing wallclock time
		sort( evn.begin(), evn.end(), SortProfLogAscWallClock);

		//compute total time
		double dt_total = 0.f;
		for(auto it = evn.begin(); it != evn.end(); ++it) { dt_total += it->dt; }

		//report
		double dt_cumsum = 0.f;
		for (auto it = evn.begin(); it != evn.end(); ++it) {
			dt_cumsum += it->dt;
			auto cat = categories.find(it->typ);
			csvlog << it->what << ";" << it->i;
			if ( cat != categories.end() )
				csvlog << ";" << cat->second;
			else
				csvlog << ";" << "APT_XX";
			auto par = parallelism.find(it->pll);
			if ( par != parallelism.end() )
				csvlog << ";" << par->second;
			else
				csvlog << ";" << "APT_IS_UNKNOWN";

			if ( it->virtualmem != MEMORY_NOSNAPSHOT_TAKEN && it->residentmem != MEMORY_NOSNAPSHOT_TAKEN )
				csvlog << ";" << it->virtualmem << ";" << it->residentmem;
			else
				csvlog << ";" << "" << ";" << "";

			csvlog << ";" << it->dt << ";" << dt_cumsum << ";" << (it->dt / dt_total) << "\n";
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local processing files");
	}
}

//program profiling should use double precision in general as MPI_Wtime() and omp_get_wtime() fires in double precision
