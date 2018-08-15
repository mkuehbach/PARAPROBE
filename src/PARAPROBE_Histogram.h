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


#ifndef __PARAPROBE_HISTOGRAM_H__
#define __PARAPROBE_HISTOGRAM_H__


#include "PARAPROBE_PeriodicTable.h"

class histogram
{
public :
	histogram()
	{
		cnts_lowest = 0.f;
		cnts_highest = 0.f;
		binstart = 0.f;
		binwidth = 0.f;
		binend = 0.f;
		tmp = 0.f;
		nbins = 0;
		valid = false;
	}
	~histogram() {}
	/*histogram( const apt_real _start, const apt_real _end, const unsigned int _nbins )
	{

		cnts_lowest = 0.f;
		cnts_highest = 0.f;
		if ( _start < _end && _nbins > 0) {
			binstart = _start;
			binwidth = 0.f;
			binend = _end;
			tmp = 0.f;
			nbins = _nbins;
			try {
				cnts.clear();
				cnts.reserve(nbins);
				for(unsigned int b = 0; b < nbins; ++b){
					cnts.push_back(0.f);
				}
			}
			catch (bad_alloc &croak) {
				stopping("Unable to allocate memory in histogram class object");
			}
			binwidth = (binend - binstart) / static_cast<apt_real>(nbins);
			tmp = 1.0 / (binend-binstart) * static_cast<apt_real>(nbins);
			valid = true; //knows limits a
		}
		else {
			binstart = 0.f;
			binwidth = 0.f;
			binend = 0.f;
			tmp = 0.f;
			nbins = 0;
			valid = false;
		}
	}*/

	histogram( const apt_real _start, const apt_real _width, const apt_real _end )
	{
		cnts_lowest = 0.f;
		cnts_highest = 0.f;
		if ( _width > EPSILON ) {
			size_t _nbins = (_end - _start) / _width;
//cout << "Histogram construction [" << _start << ":" << _width << ":" << _end << "] " << _nbins << " bins" << endl;
			if ( _nbins <= UINT32MX ) {
				binstart = _start;
				binwidth = _width;
				binend = _end;
				nbins = static_cast<unsigned int>(_nbins); //type cast from size_t to uint32 safe
				try {
					cnts.clear();
					cnts.reserve(nbins);
					for ( unsigned int b = 0; b < nbins; ++b )
						cnts.push_back(0.f);
				}
				catch (bad_alloc &croak) {
					stopping("Unable to allocate memory in histogram class object");
				}
				tmp = 1.0 / (binend-binstart) * static_cast<apt_real>(nbins);
				valid = true;
			}
		}
		else {
			binstart = 0.f;
			binwidth = 0.f;
			binend = 0.f;
			tmp = 0.f;
			nbins = 0;
			valid = false;
		}
	}
	//void construct( const apt_real _start, const apt_real _width, const apt_real _end );

	inline void add( const apt_real x )
	{
		//increment cnts in specific bin at correct bin
		if ( x >= binstart && x <= binend ) {
			apt_real rb = (x-binstart) * tmp;
			unsigned int b = static_cast<unsigned int>(floor(rb));
			if ( b < nbins )
				cnts.at(b) += 1.0;
			else
				cnts_highest += 1.0;
		}
		else {
			if ( x < binstart ) {
				cnts_lowest += 1.0;
			}
			if ( x > binend ) {
				cnts_highest += 1.0;
			}
		}
	}

	inline void add_nodump( const apt_real x )
	{
		//assumes increment cnts in specific bin at correct bin
		if ( x >= binstart && x <= binend ) {
			apt_real rb = (x-binstart) * tmp;
			unsigned int b = static_cast<unsigned int>(floor(rb));
			if ( b < nbins ) { //##MK::correct for RDF ?
				cnts.at(b) += 1.0;
//cout << b << ";" << setprecision(32) << x << ";" << rb << ";" << floor(rb) << endl;
			}
		}
	}

	inline void normalize()
	{
		apt_real tcnts = 0.0;
		for ( unsigned int b = 0; b < nbins; ++b )
			tcnts += cnts.at(b); //get total counts
		if ( tcnts >= 1.0 ) { //at least one count
			for ( unsigned int b = 0; b < nbins; ++b )
				cnts.at(b) /= tcnts;
		}
	}

	inline apt_real report( const unsigned int b )
	{
		if ( b < nbins )
			return cnts.at(b);
		else
			return 0.f;
	}
	inline apt_real start()
	{
		return binstart;
	}
	inline apt_real width()
	{
		return binwidth;
	}
	inline apt_real end()
	{
		return binend;
	}
	inline unsigned int bincount()
	{
		return nbins;
	}

	apt_real cnts_lowest;
	vector<apt_real> cnts;
	apt_real cnts_highest;

private:
	apt_real binstart;
	apt_real binwidth;
	apt_real binend;
	apt_real tmp;
	unsigned int nbins;
	bool valid;
};


#endif
