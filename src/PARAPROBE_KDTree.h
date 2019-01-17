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

#ifndef __PARAPROBE_H_KDTREE_H__
#define __PARAPROBE_H_KDTREE_H__

#include "PARAPROBE_SpaceBucketing.h"


//##MK::implement KDTree

//using namespace std;
//MK::SINGLE_PRECISION for APT data suffices
//typedef float apt_xyz;


inline apt_xyz pi()
{
	return static_cast<apt_xyz>(PI);
}

inline apt_xyz mysqr(apt_xyz v)
{
	return v*v;
}


inline apt_xyz euclidean_sqrd(const p3dm1 & a, const p3dm1 & b) {

	return mysqr(a.x-b.x) + mysqr(a.y-b.y) + mysqr(a.z-b.z);
}



struct cuboid
{
	p3d min;
	p3d max;
	cuboid() : min(RMAX,RMAX,RMAX), max(RMIN,RMIN,RMIN) {}
	cuboid(const p3d & thiscorner, const p3d & thatcorner) :
		min(thiscorner), max(thatcorner) {}
	apt_xyz outside_proximity(const p3dm1 & p ) const
	{
		apt_xyz res = 0.0; //assuming we are inside
		if ( p.x < this->min.x )
			res += mysqr(this->min.x - p.x);
		else if (p.x > this->max.x )
			res += mysqr(p.x - this->max.x);

		if ( p.y < this->min.y )
			res += mysqr(this->min.y - p.y);
		else if (p.y > this->max.y )
			res += mysqr(p.y - this->max.y);

		if ( p.z < this->min.z )
			res += mysqr(this->min.z - p.z);
		else if (p.z > this->max.z )
			res += mysqr(p.z - this->max.z);
		return res;
	}
};

//MK::https://programmizm.sourceforge.io/blog/2011/a-practical-implementation-of-kd-trees
struct node
{
	size_t i0, i1, split; //, dimdebug;dimdebug(-1),
	apt_xyz splitpos;

	node() : i0(-1), i1(-1), split(-1), splitpos(numeric_limits<apt_xyz>::max()) {}

	friend bool is_leaf(const node & n) {
		return n.split == size_t(-1);
	}
	friend size_t npoints(const node & n) {
		return (n.i1 - n.i0) + 1;
	}
};


struct build_task
{
	size_t first, last, node, dim;
};


struct traverse_task
{
	cuboid bx;
	size_t node;
	size_t dim;
	apt_xyz d;
};


struct scan_task
{
	cuboid bx;
	size_t node;
	size_t dim;
};


struct kd_tree
{
	p3d min;
	p3d max;
	vector<node> nodes;

	kd_tree() : min(p3d()), max(p3d()) {}
	void build(const vector<p3d> & points, vector<size_t> & permutate );

	void pack_p3dm1_d( const vector<size_t> & permutate, const vector<p3dm1> & apt1, const vector<apt_xyz> & dist1, vector<p3dm1> & apt2, vector<apt_xyz> & dist2 );
	void pack_p3dm1( const vector<size_t> & permutate, const vector<p3d> & aptpos1, const vector<unsigned int> & aptlabel1, vector<p3dm1> & apt2 );
	void pack_p3dm1_dist( vector<size_t> const & permutate, vector<p3dm1> const & in1, vector<apt_xyz> const & in2, vector<p3dm1> & out1, vector<apt_xyz> & out2 );

	p3dm1 nearest_external(const p3dm1 target, const vector<p3dm1> & sortedpoints, apt_xyz epsball ) const;
	p3dm1 nearest(const size_t idx, const vector<p3dm1> & sortedpoints, apt_xyz epsball ) const;

	void range_rball_noclear_nosort_external(const p3dm1 target, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<nbor> & result );
	void range_rball_noclear_nosort(const size_t idx, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<nbor> & result );
	void range_rball_noclear_nosort_external_p3dm1(const p3dm1 target, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<p3dm1> & result );
	void range_rball_noclear_nosort_p3dm1(const size_t idx, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<p3dm1> & result );
	void range_rball_noclear_nosort_p3d( const p3dm1 target, vector<p3dm1> const & sortedpoints, const apt_xyz radius_sqrd, vector<p3dm1> & result );


	void range_rball_noclear_nosort_indices(const size_t idx, const vector<p3dm1> & sortedpoints, const apt_xyz radius_sqrd, vector<size_t> & result );

	inline p3d get_min();

	inline p3d get_max();

	bool verify(const vector<p3dm1> & sortedpoints);

	void get_allboundingboxes( vector<scan_task> & out );

	void display_nodes();

	size_t get_treememory_consumption();
};


inline void expand(p3d & min, p3d & max, const p3dm1 & p)
{
	if (p.x < min.x)	min.x = p.x;
	if (p.x > max.x) 	max.x = p.x;
	if (p.y	< min.y) 	min.y = p.y;
	if (p.y > max.y) 	max.y = p.y;
	if (p.z	< min.z) 	min.z = p.z;
	if (p.z > max.z) 	max.z = p.z;
}

inline void expand(p3d & min, p3d & max, const p3d & p)
{
	if (p.x < min.x)	min.x = p.x;
	if (p.x > max.x) 	max.x = p.x;
	if (p.y	< min.y) 	min.y = p.y;
	if (p.y > max.y) 	max.y = p.y;
	if (p.z	< min.z) 	min.z = p.z;
	if (p.z > max.z) 	max.z = p.z;
}


struct lower_x
{
	const vector<p3d> & points;

	lower_x(const vector<p3d> & p) : points(p) {}

	bool operator()(size_t i1, size_t i2) const
	{
		return points[i1].x < points[i2].x;
	}
};

struct lower_y
{
	const vector<p3d> & points;

	lower_y(const vector<p3d> & p) : points(p) {}

	bool operator()(size_t i1, size_t i2) const
	{
		return points[i1].y < points[i2].y;
	}
};

struct lower_z
{
	const vector<p3d> & points;

	lower_z(const vector<p3d> & p) : points(p) {}

	bool operator()(size_t i1, size_t i2) const
	{
		return points[i1].z < points[i2].z;
	}
};


/* HOW TO USE
int main()
{
//g++ -Wall kdtree01.cpp -o kdtree
	cout << "Hello, World!" << setprecision(18) << scientific << endl;

//generate random test data on unit cube
	mt19937 rd(12345);
	uniform_real_distribution<apt_xyz> unifrnd(0,1);

	size_t N = NN;
	//MK::separate range labels and ion positions to improve cache utilization on aptpos1 during KDTree build
	vector<p3d> aptpos1;
	vector<unsigned int> aptlabel1;
	vector<size_t> permutations;
	vector<p3dm1> apt2;
	try {
		aptpos1.reserve(N);
		aptlabel1.reserve(N);
		permutations.reserve(N);

		apt2.reserve(N);
	}
	catch (bad_alloc &croak) {
		return 0;
	}

	for (size_t i = 0; i < N; ++i) {
		aptpos1.push_back( p3d(unifrnd(rd), unifrnd(rd), unifrnd(rd) ) );
		aptlabel1.push_back( 1 );
//cout << aptpos1.back().x << ";" << aptpos1.back().y << ";" << aptpos1.back().z << "\t\t" << aptlabel1.back() << endl;
	}

	cout << "Fake APT done " << aptpos1.size() << endl;

	kd_tree mytree( aptpos1, permutations );

	mytree.pack_p3dm1( permutations, aptpos1, aptlabel1, apt2 );

	//permutation array no longer necessary
	//##MK::permutations.swap( vector<size_t>() );

	if ( mytree.verify( apt2 ) == true )
		cout << "KDTree is valid!" << endl;
	else
		cout << "KDTree has overlapping indices!" << endl;


cout << "Utilizing tree for nearest neighbor" << endl;
	for(size_t ii = 0; ii < N; ++ii) {
		p3dm1 candidate = mytree.nearest( ii, apt2, 1.0e-4 );
		p3dm1 thisone = apt2[ii];
cout << ii << "\t\t" << sqrt(euclidean_sqrd(thisone, candidate)) << endl;

		vector<nbor> these;
		apt_xyz RadiusSQR = sqr(0.05);

		mytree.range_rball_noclear_nosort(ii, apt2, RadiusSQR, these );

		sort( these.begin(), these.end(), nbsort );
cout << ii << "\t\t" << these.size() << "\t\t" << these.at(0).d << "\t\t" << these.back().d << endl;

		//for(size_t j = 0; j < these.size(); ++j) {
		//	cout << these[j].d << " ";
		//} cout << endl;

//
//		//##MK::DEBUG BRUTE FORCE VALIDATION
//		//brute force
//		apt_xyz bruteforce = numeric_limits<apt_xyz>::max();
//
//		for(size_t j = 0; j < ii; ++j) {
//			apt_xyz current = euclidean_sqrd(thisone, aptpos2[j]);
//			if ( current > bruteforce )
//				continue;
//			else
//				bruteforce = current;
//		}
//		for(size_t j = ii + 1; j < N; ++j) {
//			apt_xyz current = euclidean_sqrd(thisone, aptpos2[j]);
//			if ( current > bruteforce )
//				continue;
//			else
//				bruteforce = current;
//		}
//		//if ( (euclidean_sqrd(thisone, candidate) - bruteforce) > 1.0e-12 ) {
//cout << ii << "\t\t" << euclidean_sqrd(thisone, candidate) << "\t\t" << bruteforce << endl;
//		//}
//		//##MK::DEBUG BRUTE FORCE VALIDATION
//
	}


	//mytree.display_nodes();

//###MK::open problems::finalize nearest
//###MK::

	//##https://programmizm.sourceforge.io/blog/2011/nearest-neighbor-search-using-kd-trees

	return 0;
}
*/

#endif
