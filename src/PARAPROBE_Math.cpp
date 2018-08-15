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

#include "PARAPROBE_Math.h"


/* 
	Distance Between Point and Triangle in 3D
	David Eberly
	Geometric Tools, LLC
	http://www.geometrictools.com/
	Copyright 
	c 1998-2016. All Rights Reserved.
	Created: September 28, 1999
	Last Modified: March 1, 2008
	vector3 closesPointOnTriangle( const vector3 *triangle, const vector3 &sourcePosition )
*/

apt_real mathHdl::lerp(const apt_real v0, const apt_real v1, const apt_real t)
{
	return (1.f - t)*v0 + t*v1;
}


vector<apt_real> mathHdl::quantiles_nosort( vector<apt_real> const & in, vector<apt_real> const & q )
{
	if ( in.size() > 1 ) { //most likely more than one or nothing
		//##MK::if desired use random sub-sampling here
		//vector<apt_real> data = in;

		//##MK::can be improved for instance if only quantiles << 1 are sought via using n_th element
		//sort(data.begin(), data.end());

		vector<apt_real> quants;
		for (size_t i = 0; i < q.size(); ++i)
		{
			apt_real poi = lerp(-0.5, static_cast<apt_real>(in.size()) - 0.5, q.at(i));

			size_t left = max( static_cast<int64_t>(floor(poi)), static_cast<int64_t>(0) );
			size_t right = min( static_cast<int64_t>(ceil(poi)), static_cast<int64_t>(in.size() - 1) );

			apt_real datLeft = in.at(left);
			apt_real datRight = in.at(right);

			apt_real quantile = lerp(datLeft, datRight, poi - static_cast<apt_real>(left) );

			quants.push_back(quantile);
		}

		return quants;
	}
	 else if ( in.size() == 1) {
		return vector<apt_real>(1, in[0]);
	 }
	 else {
		 return vector<apt_real>();
	 }
}


apt_real mathHdl::closestPointOnTriangle( const tri3d face, const p3d src )
{
	v3d edge0( face.x2-face.x1, face.y2-face.y1, face.z2-face.z1 );
	v3d edge1( face.x3-face.x1, face.y3-face.y1, face.z3-face.z1 );
	v3d v0( face.x1-src.x, face.y1-src.y, face.z1-src.z );

	apt_real a = SQR(edge0.u) + SQR(edge0.v) + SQR(edge0.w); //edge0.dot( edge0 );
	apt_real b = edge0.u*edge1.u + edge0.v*edge1.v + edge0.w*edge1.w; //edge0.dot( edge1 );
	apt_real c = SQR(edge1.u) + SQR(edge1.v) + SQR(edge1.w); //edge1.dot( edge1 );
	apt_real d = edge0.u*v0.u + edge0.v*v0.v + edge0.w*v0.w; //edge0.dot( v0 );
	apt_real e = edge1.u*v0.u + edge1.v*v0.v + edge1.w*v0.w; //edge1.dot( v0 );

	apt_real det = a*c - b*b;
	apt_real s = b*e - c*d;
	apt_real t = b*d - a*e;

	if ( s + t < det ) {
		if ( s < 0.f ) {
			if ( t < 0.f ) { //region 4
				if ( d < 0.f ) {
					s = CLAMP( -d/a, 0.f, 1.f );
					t = 0.f;
				}
				else {
					s = 0.f;
					t = CLAMP( -e/c, 0.f, 1.f );
				}
			}
			else {//region 3
				s = 0.f;
				t = CLAMP( -e/c, 0.f, 1.f );
			}
		}
		else if ( t < 0.f ) { //region 5
			s = CLAMP( -d/a, 0.f, 1.f );
			t = 0.f;
		}
		else { //region 0
			apt_real invDet = 1.f / det;
			s *= invDet;
			t *= invDet;
		}
	}
	else {
		if ( s < 0.f ) { //region 2
			apt_real tmp0 = b+d;
			apt_real tmp1 = c+e;
			if ( tmp1 > tmp0 ) {
				apt_real numer = tmp1 - tmp0;
				apt_real denom = a-2*b+c;
				s = CLAMP( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else {
				t = CLAMP( -e/c, 0.f, 1.f );
				s = 0.f;
			}
		}
		else if ( t < 0.f ) { //region 6
			if ( a+d > b+e ) {
				apt_real numer = c+e-b-d;
				apt_real denom = a-2*b+c;
				s = CLAMP( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else {
				s = CLAMP( -e/c, 0.f, 1.f );
				t = 0.f;
			}
		}
		else { //region 1
			apt_real numer = c+e-b-d;
			apt_real denom = a-2*b+c;
			s = CLAMP( numer/denom, 0.f, 1.f );
			t = 1.f - s;
		}
	}

	//closest point is cp
	p3d cp( face.x1 + s*edge0.u + t*edge1.u, face.y1 + s*edge0.v + t*edge1.v, face.z1 + s*edge0.w + t*edge1.w );

	//so return SQR of difference to avoid sqrt in the code
	return ( SQR(cp.x-src.x) + SQR(cp.y-src.y) + SQR(cp.z-src.z) );
}
