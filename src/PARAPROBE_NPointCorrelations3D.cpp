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


#include "PARAPROBE_NPointCorrelations3D.h"

//implements functionality to compute two-point correlation functions in spherical environment about ion location supporting
//by 3d binning with binwidth _width to identify directionality dependent features and perform methods as developed by Kalidindi et al.


ostream& operator<<(ostream& in, npnbor const & val)
{
	//in << "DiffVec/Distance/MarkValue = " << val.diffvector.u << ";" << val.diffvector.v << ";" << val.diffvector.w << "---" << val.d << "\t\t" << val.m << endl;
	in << "DiffVec/Distance = " << val.diffvector.u << ";" << val.diffvector.v << ";" << val.diffvector.w << "---" << val.d << endl;

	return in;
}



npc3d::npc3d()
{
	support = sqb(); //##MK::physically pointless values
	fvec_maxlen = 0.f;
	center_offset = vxxl();
}


npc3d::npc3d( const apt_real _vfeature_maxlen, const apt_real _width, const apt_real initval, const bool allocate )
{
	//probe first of all expected number of bins, bin aggregate is a cube with odd number of bins along edge, cubic bins of width _width
	//MK::uses right-handed coordinate system with origin centered at ((support.ni - 1) + 0.5)*support.width

	apt_real realhalfedgelen = ceil(_vfeature_maxlen / _width ) + 1.f; //MK: +1 for one layer guardzone
	size_t halfedgelen = realhalfedgelen;
	size_t edgelen = 1 + 2 * halfedgelen; //two halfedgelen plus center bin

cout << "Computing n-point correlations with halfedgelen/edgelen " << halfedgelen << "/" << edgelen << endl;

	support.nx = edgelen;
	support.ny = edgelen;
	support.nz = edgelen;
	support.nxy = support.nx * support.ny;
	support.nxyz = support.nx * support.ny * support.nz;
	support.width = _width;
	support.box.xmi = -1.f * (halfedgelen + 0.5) * support.width; //MK:: local continuum coordinates relative to center (0,0,0) about target/central ions
	support.box.xmx = +1.f * (halfedgelen + 0.5) * support.width;
	support.box.ymi = -1.f * (halfedgelen + 0.5) * support.width;
	support.box.ymx = +1.f * (halfedgelen + 0.5) * support.width;
	support.box.zmi = -1.f * (halfedgelen + 0.5) * support.width;
	support.box.zmx = +1.f * (halfedgelen + 0.5) * support.width; //cubic bin grid about central bin
	support.box.scale();

	fvec_maxlen = _vfeature_maxlen; //maximum length of feature vector

	center_offset = vxxl( 	static_cast<int>(halfedgelen),
							static_cast<int>(halfedgelen),
							static_cast<int>(halfedgelen)  );

	if ( allocate == true ) {
		try {
			cnts.reserve( support.nxyz );
			for( size_t i = 0; i < support.nxyz; ++i )
				cnts.push_back( 0 ); //initval );
		}
		catch (std::bad_alloc &npc3dcroak) {
			cout << "Allocation error in npc3d" << endl;
			return;
		}
	}
}


npc3d::~npc3d()
{
	//MK::vector clears after itself
}


/*
void npc3d::add( d3d where, const apt_real val )
{
	//MK::add only if within maximum length of feature vector
	//MK::where uvw coordinates can be negative translation into local bin coordinate system necessary
	apt_real tvec_len = where.diff_vector_len();

	int NX = support.nx;
	int NXY = support.nxy;
	apt_real _width = 1.f / support.width;

	if ( tvec_len <= fvec_maxlen ) {
		//MK::cast issues http://jkorpela.fi/round.html
		//MK::https://www.cs.cmu.edu/~rbd/papers/cmj-float-to-int.html
		//##MK::binning like so int discretelyhere = there, we get concentration at the 0 component values symmetrical to coordinates undesired

		int targetbin = 0;
		//apt_real there = (where.u < 0.f ) ? ceil(where.u / support.width) : floor(where.u / support.width); //floor(-2.3) is -3.0 but floor(2.3) is 2 so to make binning symmetric
		apt_real there = where.u * _width;
		//int discretelyhere = there;
		int discretelyhere = (there >= 0.f) ? (int)(there+0.5) : (int)(there-0.5);
		int xx = discretelyhere + center_offset.x;
		targetbin = targetbin + (discretelyhere + center_offset.x) * 1;
//cout << where.u << "\t\t" << there << "\t\t" << targetbin << endl;

		//there = (where.v < 0.f ) ? ceil(where.v / support.width) : floor(where.v / support.width);
		there = where.v * _width;
		//discretelyhere = there;
		discretelyhere = (there >= 0.f) ? (int)(there+0.5) : (int)(there-0.5);
		int yy = discretelyhere + center_offset.y;
		targetbin = targetbin + (discretelyhere + center_offset.y) * NX;
//cout << where.v << "\t\t" << there << "\t\t" << targetbin << endl;

		//there = (where.w < 0.f ) ? ceil(where.w / support.width) : floor(where.w / support.width);
		there = where.w * _width;
		//discretelyhere = there;
		discretelyhere = (there >= 0.f) ? (int)(there+0.5) : (int)(there-0.5);
		int zz = discretelyhere + center_offset.z;
		targetbin = targetbin + (discretelyhere + center_offset.z) * NXY;

		targetbin = xx + yy*NX + zz*NXY;
//cout << where.w << "\t\t" << there << "\t\t" << targetbin << endl;


//		if ( xx == (NX-1)/2 || yy == (NX-1)/2 || zz == (NX-1)/2 ) {
//	 	 	#pragma omp critical
//			{
//				cout << "uvw/xxyyzz/dhw/bin = " << where.u << ";" << where.v << ";" << where.w << "\t\t" << xx << ";" << yy << ";" << zz << "\t\t" << discretelyhere << "---->" << targetbin << endl;
//				//cout << cnts.at(targetbin) << endl;
//			}
//		}

		cnts[targetbin] = cnts[targetbin] + val;
	}
}
*/


void npc3d::add( npnbor const & where )
{
	//MK::add only if within maximum length of feature vector
	//MK::where uvw coordinates can be negative translation into local bin coordinate system necessary
	if ( where.d <= fvec_maxlen ) {

		int NX = support.nx;
		int NXY = support.nxy;
		apt_real _width = 1.f / support.width;

		//MK::cast issues http://jkorpela.fi/round.html
		//MK::https://www.cs.cmu.edu/~rbd/papers/cmj-float-to-int.html
		//##MK::binning like so int discretelyhere = there, we get concentration at the 0 component values symmetrical to coordinates undesired

		int targetbin = 0;
		//apt_real there = (where.u < 0.f ) ? ceil(where.u / support.width) : floor(where.u / support.width); //floor(-2.3) is -3.0 but floor(2.3) is 2 so to make binning symmetric
		apt_real there = where.diffvector.u * _width;
		//int discretelyhere = there;
		int discretelyhere = (there >= 0.f) ? (int)(there+0.5) : (int)(there-0.5);
		int xx = discretelyhere + center_offset.x;
		targetbin = targetbin + (discretelyhere + center_offset.x) * 1;
//cout << where.u << "\t\t" << there << "\t\t" << targetbin << endl;

		//there = (where.v < 0.f ) ? ceil(where.v / support.width) : floor(where.v / support.width);
		there = where.diffvector.v * _width;
		//discretelyhere = there;
		discretelyhere = (there >= 0.f) ? (int)(there+0.5) : (int)(there-0.5);
		int yy = discretelyhere + center_offset.y;
		targetbin = targetbin + (discretelyhere + center_offset.y) * NX;
//cout << where.v << "\t\t" << there << "\t\t" << targetbin << endl;

		//there = (where.w < 0.f ) ? ceil(where.w / support.width) : floor(where.w / support.width);
		there = where.diffvector.w * _width;
		//discretelyhere = there;
		discretelyhere = (there >= 0.f) ? (int)(there+0.5) : (int)(there-0.5);
		int zz = discretelyhere + center_offset.z;
		targetbin = targetbin + (discretelyhere + center_offset.z) * NXY;

		targetbin = xx + yy*NX + zz*NXY;
//cout << where.w << "\t\t" << there << "\t\t" << targetbin << endl;

		cnts[targetbin]++;
	}
}


/*
void npc3d::normalize()
{
	apt_real sum = 0.f;
	for( size_t b = 0; b < cnts.size(); ++b )
		sum += cnts[b];

	if ( sum >= 1.f ) { //at least one count
		for( size_t b = 0; b < cnts.size(); ++b )
			cnts[b] /= sum;
	}
}
*/


void npc3d::report()
{
	cout << "n-point correlations" << endl;
	cout << support << endl;
}


sqb npc3d::get_support()
{
	return support;
}


p3d npc3d::get_bincenter( const size_t ix, const size_t iy, const size_t iz )
{
	//the npc3d histogram bins belong to an 2n+1 aggregate where the center cell
	//(i.e. cell position n,n,n) contains the origin of a 3d Euclidean coordinate system
	//if ( ix < this->support.nx && iy < this->support.ny && iz < this->support.nz ) {
	return p3d(
			static_cast<apt_real>(ix)*this->support.width - static_cast<apt_real>((this->support.nx-1)/2),
			static_cast<apt_real>(iy)*this->support.width - static_cast<apt_real>((this->support.ny-1)/2),
			static_cast<apt_real>(iz)*this->support.width - static_cast<apt_real>((this->support.nz-1)/2)
			);
}

