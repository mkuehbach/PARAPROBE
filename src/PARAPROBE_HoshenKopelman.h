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


#ifndef __PARAPROBE_HOSHENKOPELMAN_H__
#define __PARAPROBE_HOSHENKOPELMAN_H__

#include "PARAPROBE_KDTree.h"


//data structure for path/compressed union-find
class UF
{
	//##MK::potentially encapsulate
public:
	unsigned int* id;
	unsigned int* sz;

	// Create an empty union find data structure with N isolated sets.
	UF() : id(NULL), sz(NULL), stillgood(false) {}
	UF( const unsigned int N ) {
		stillgood = true;
		id = NULL;
		try { id = new unsigned int[N]; }
		catch (std::bad_alloc &hkexc) {
			stopping( "HoshenKopelman Union/Find Algorithm unable to hash container id" );
			stillgood = false;
		}
		sz = NULL;
		try { sz = new unsigned int[N]; }
		catch (std::bad_alloc &hkexc) {
			stopping( "HoshenKopelman Union/Find Algorithm unable to hash container sz" );
			stillgood = false;
		}
		for( unsigned int i = 0; i < N; ++i ) {
			id[i] = i;
			sz[i] = 1;
		}
	}
	~UF() {
		if ( id != NULL ) {	delete [] id; id = NULL; }
		if ( sz != NULL ) { delete [] sz; sz = NULL; }
		stillgood = false;
	}
	
	//Return the id of component corresponding to object p.
	unsigned int initialAssgn( const unsigned int i) {
		return i;
	}
	unsigned int find_root(unsigned int p)	{
		unsigned int root = p;
		while (root != id[root])
			root = id[root]; //traversal operation because for a root node i == id[i] holds
		while (p != root) { //path compression
			unsigned int newp = id[p];
			id[p] = root;
			p = newp;
		}
		return root;
	}

	//Replace sets containing x and y with their union, MK::was originally void
	unsigned int merge( unsigned int x, unsigned int y) {
		unsigned int i = find_root(x);
		unsigned int j = find_root(y);
		if ( i == j ) return i; //MK::return simply one...

		// make smaller root point to larger one
		if (sz[i] < sz[j])	{ //eq class j larger than eq class i
			id[i] = j; 
			sz[j] += sz[i];
			return j;
		}
		else { //eq class j smaller or equal to eq class i
			id[j] = i; 
			sz[i] += sz[j];
			return i;
		}
	}
	bool still_good( void ) {
		return stillgood;
	}

private:
	bool stillgood;
};


struct percAnalysis
{
	unsigned int nClusterTotal;
	unsigned int LargestClusterCnt;
	unsigned int LargestClusterID;
	percAnalysis() : nClusterTotal(0), LargestClusterCnt(0), LargestClusterID(0) {}
	percAnalysis(const unsigned int _nct, const unsigned int _lcnt, const unsigned int _lid) :
		nClusterTotal(_nct), LargestClusterCnt(_lcnt), LargestClusterID(_lid) {}
	//~percAnalysis(){}
};


struct loginfo_perc
{
	//MPI_Wtimers are double precision per se!
	double ProfInitializing;
	double ProfHoshenKopeling;
	double ProfCompactifying;
	double ProfCheckLabeling;
	double ProfCharacterizing;
	double ProfPercTotal;
	loginfo_perc() : ProfInitializing(0.0), ProfHoshenKopeling(0.0), ProfCompactifying(0.0), ProfCheckLabeling(0.0),
			ProfCharacterizing(0.0), ProfPercTotal(0.0) {}
	loginfo_perc(const double _init, const double _hk, const double _cmptfy, const double _chklbl, const double _chr, const double _tot) :
		ProfInitializing(_init), ProfHoshenKopeling(_hk), ProfCompactifying(_cmptfy), ProfCheckLabeling(_chklbl), 
			ProfCharacterizing(_chr), ProfPercTotal(_tot) {}
	//~loginfo_perc(){}
};


class percAnalyzer
{
public:
	bool initialize( const bool* inputdata, const unsigned int nxx, const unsigned int nyy, const unsigned int nzz );
	void hk3d_core_nonperiodic( const unsigned int x, const unsigned int y, const unsigned int z );
	bool hoshen_kopelman( void );
	bool compactify( void );
	bool checkLabeling( void );
	bool determine_clustersize_distr();
	bool rebinarize( const unsigned int target, const unsigned int nsz, bool* bitmap );

	percAnalyzer()
	{
		idd = NULL;
		finder = NULL;
		setGridsize( 0, 0, 0 );
	};

	~percAnalyzer()
	{
		delete [] idd; idd = NULL;
		delete finder; finder = NULL;
	};

	void setGridsize( const unsigned int nx, const unsigned int ny, const unsigned int nz )
	{
		NX = nx;
		NY = ny;
		NZ = nz;
		NXY = NX*NY;
		NXYZ = NXY*NZ;
	}

	inline unsigned int getNCluster( void )
	{
		return results.nClusterTotal;
	}

	inline unsigned int getClusterSize( void )
	{
		return results.LargestClusterCnt;
	}

private:
	unsigned int* idd;
	UF* finder;

	unsigned int NX;
	unsigned int NY;
	unsigned int NZ;
	unsigned int NXY;
	unsigned int NXYZ;

	struct percAnalysis results;
};


#endif
