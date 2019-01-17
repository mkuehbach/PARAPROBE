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

#ifndef __PARAPROBE_PERIODICTABLE_H__
#define __PARAPROBE_PERIODICTABLE_H__

//#include "PARAPROBE_Profiling.h"
#include "PARAPROBE_IntelMKL.h"


class TypeSpecDescrStat;

struct mqrange
{
	apt_real lo;
	apt_real hi;
	mqrange() :
		lo(0.f), hi(0.f) {}
	mqrange( const apt_real _lo, const apt_real _hi) : lo(_lo), hi(_hi) {}
	//~mqrange(){}

	inline apt_real width();
	inline bool inrange(const apt_real val);
};


inline bool SortMQRangeAscending( const mqrange &aa1, const mqrange &aa2)
{
	//promoting sorting in descending order
	return aa1.hi < aa2.lo;
}


struct nuclid
{
	unsigned int id; //##MK::so far a dummy only but can be extended by adding pieces of information by user
	nuclid() : id(0) {}
	nuclid(const unsigned int _id) : id(_id) {}
};


class PeriodicTable
{
	//MK::implements a smart periodic table with which to handle range files and ion type identification
	friend class TypeSpecDescrStat;
public:
	PeriodicTable();
	~PeriodicTable();

	void load_periodictable();
	unsigned int mq2ionname( const apt_real mq );
	unsigned int ion_name2type( const string name );
	string ion_type2name( const unsigned int type );
	unsigned int get_maxtypeid();

	void add_iontype_single_or_molecular( vector<string> const & in, mqrange const & ival );

	bool read_rangefile( string ascii_io_fn );
	bool read_rangefile2( string ascii_io_fn );

private:
	map<string, nuclid> NuclidTable;
	map<string, unsigned int> IonTypes; //two way referencing to ease ranging
	map<unsigned int, string> IonNames;
	unsigned int MaximumTypID;

	//vector<string> IonName;
	vector<vector<mqrange>> MQ2IonName;
	bool rangefile_loaded;
};


class TypeSpecDescrStat
{
	//class object to ease the processing of generic iontype-specific descriptive statistics characterization tasks
	//key idea: take for instance nearest neighbor distributions, APT folks are interested in studying not only how
	//an ion clusters is so and so far from another of any type but specific type
	//hence we should be able to process with the same function and code queries like
	//NN Al against all types
	//NN Al against only Sc
	//NN Al against only Sc and Zr
	//now lets abbreviate these queries via shorthand codes
	// Al,
	// Al,Sc
	// Al,Sc,Zr
	//##MK::add also negations as such Al,!Sc meaning Al against every type excluding Sc


public:
	TypeSpecDescrStat();
	~TypeSpecDescrStat();

	bool define_iontask3(const string command, PeriodicTable const & pse );
	bool define_iontask4(const string command, PeriodicTable const & pse );


	//pair<string, unsigned int> target;	 	 //##MK::deprecated, the typeID of the specimen for which do perform spatial statistics, i.e. Al e.g. 0
	map<string, unsigned int> trgcandidates; //the typeIDs of the central ions for which to perform spatial statistics
	map<string, unsigned int> envcandidates; //the corresponding typeIDs to test against in the environment of the central ions
};


void parse_tasks( const string command, vector<TypeSpecDescrStat> & these, PeriodicTable const & thispse ); //parses the command into individual iontype related tasks storing results as task class objects in these


#endif
