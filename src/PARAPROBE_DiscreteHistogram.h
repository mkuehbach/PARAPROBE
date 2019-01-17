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


#ifndef __PARAPROBE_DISCRETE_HISTOGRAM_H__
#define __PARAPROBE_DISCRETE_HISTOGRAM_H__


#include "PARAPROBE_Histogram.h"


class discrhistogram
{
public :
	discrhistogram(){}
	~discrhistogram(){}

	//MK::implementation 1 : log(N) access time, but more memory efficient when bin IDs range is occupied sparse
	/*
	inline void add( const size_t b, const apt_real val )
	{
		auto it = cnts.find(b);
		if ( it != cnts.end() ) { //element exists
			it->second += val;
		}
		else { //element not yet exists, create and assign
			cnts[b] = val;
		}
	}

	inline void normalize()
	{
		apt_real sum = 0.f;
		for ( auto it = cnts.begin(); it != cnts.end(); ++it )
			sum += it->second;

		if ( sum >= 1.0 ) {
			for ( auto it = cnts.begin(); it != cnts.end(); ++it )
				it->second = it->second / sum;
		}
	}

	inline apt_real report( const size_t b )
	{
		auto it = cnts.find(b);
		if ( it != cnts.end() )
			it->second;
		else
			return 0.f;
	}

	map<size_t, apt_real> cnts;
	*/

	//MK::implementation 2 : O(1) constant access time but potentially memory inefficient if bins ID range are filled sparsely
	inline void add( const size_t b, const apt_real val )
	{
		//STL vector::resize(size_type n, value_type val --> if n is greater than the current container size,
		//the content is expanded by inserting at the end as many elements as needed to reach a size of n.
		//If val is specified, the new elements are initialized as copies of val, otherwise, they are value-initialized.
		if ( b < cnts.size() ) { //bin exists
			cnts.at(b) += val;
		}
		else { //bin with ID b does not yet exist, mind C++ style zero-indexing i.e. b=0 is first value of cnts, b=1 second and so forth
			try {
				cnts.resize( b + 1, 0.f ); //if for instance cnts.size() was 2 and b should be 4 we want to be able to access cnts[4] safely, i.e. cnts.size() needs to be 4+1
				cnts.at(b) += val;
			}
			catch (std::bad_alloc &discrhistcroak) {
				cout << "Allocation error during resizing discrete histogram" << endl;
			}
		}
	}

	inline void normalize()
	{
		apt_real sum = 0.f;
		for ( auto it = cnts.begin(); it != cnts.end(); ++it )
			sum += *it;

		if ( sum >= 1.f ) {
			for ( auto it = cnts.begin(); it != cnts.end(); ++it )
				*it = *it / sum;
		}
	}

	inline apt_real report( const size_t b )
	{
		if ( b < cnts.size() )
			return cnts.at(b);
		else
			return 0.f;
	}

	vector<apt_real> cnts;
};


#endif
