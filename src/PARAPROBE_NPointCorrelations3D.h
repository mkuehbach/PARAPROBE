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


#ifndef __PARAPROBE_NPOINTCORR_H__
#define __PARAPROBE_NPOINTCORR_H__


#include "PARAPROBE_DiscreteHistogram.h"


struct npnbor
{
	d3d diffvector;
	apt_xyz d;
	//unsigned int m;
	npnbor() : diffvector(0.f,0.f,0.f), d(RMAX) {} //, m(UNKNOWNTYPE) {}
	npnbor(const d3d _diff, const apt_xyz _d ) :
		diffvector(_diff), d(_d) {}
	//, const unsigned int _m) : diffvector(_diff), d(_d), m(_m) {}
};

ostream& operator<<(ostream& in, npnbor const & val);

inline bool SortNPNeighborsForAscDistance(const npnbor & a, const npnbor & b)
{
	return a.d < b.d;
}


class npc3d
{
public :
	npc3d();
	npc3d( const apt_real _vfeature_maxlen, const apt_real _width, const apt_real initval, const bool allocate );
	~npc3d();

	//void add( d3d where, const apt_real val );
	void add( npnbor const & where );
	//void normalize();
	void report();

	sqb get_support();
	p3d get_bincenter( const size_t ix, const size_t iy, const size_t iz );

	//vector<apt_real> cnts;

	//##MK::if PCF is sparse maybe map is a better option
	vector<unsigned int> cnts;

private:
	sqb support;
	apt_real fvec_maxlen;
	vxxl center_offset; //min edge of bin cnts[0] is at support.box.imi which is offset by (nx-1)/2
};

#endif
