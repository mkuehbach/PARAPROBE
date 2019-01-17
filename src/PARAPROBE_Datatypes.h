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


#ifndef __PARAPROBE_DATATYPES_H__
#define __PARAPROBE_DATATYPES_H__

//#include "PARAPROBE_BSIMD.h"
#include "PARAPROBE_EPOSEndianness.h"


struct epos_event
{
	//##MK::float is 32bit by definition in C/C++
	float x;		//coordinates
	float y;
	float z;
	float m_q;		//mass-to-charge ratio
	float tof;		//time of flight
	float vdc;		//voltage DC
	float vpu;		//voltage pulse
	float detx;		//detector coordinates
	float dety;
	int nulls;		//nulls - is the time between each pulse
	int nat_pulse;	//dead count related to pulse

	epos_event() :
		x(0.f),
		y(0.f),
		z(0.f),
		m_q(0.f),
		tof(0.f),
		vdc(0.f),
		vpu(0.f),
		detx(0.f),
		dety(0.f),
		nulls(0),
		nat_pulse(0) {}
	epos_event( const float _x, const float _y, const float _z, const float _mq,
		const float _tf, const float _vdc, const float _vpu, const float _dx, const float _dy,
		const int _nu, const int _napu) : x(_x), y(_y), z(_z), m_q(_mq), tof(_tf), vdc(_vdc), vpu(_vpu), 
		detx(_dx), dety(_dy), nulls(_nu), nat_pulse(_napu) {}
};


struct runparm
{
	apt_real eta;
	apt_real kf;
	apt_real icf;
	int targetrank;
	int targetthread;
	unsigned int jobid;

	runparm() :
		eta(0.f), kf(0.f), icf(0.f), targetrank(MASTER), targetthread(MASTER), jobid(0) {}
	runparm( const apt_real _eta, const apt_real _kf, const apt_real _icf,
			const int _tr, const int _tt, const unsigned int _jid) :
			eta(_eta), kf(_kf), icf(_icf),
			targetrank(_tr), targetthread(_tt), jobid(_jid) {}
};



struct p2d
{
	apt_xyz x;
	apt_xyz y;

	p2d() : x(0.f), y(0.f) {}
	p2d(const apt_xyz _x, const apt_xyz _y) : x(_x), y(_y) {}
};

ostream& operator<<(ostream& in, p2d const & val);



struct p3d
{
	apt_xyz x;
	apt_xyz y;
	apt_xyz z;

	p3d() : x(0.f), y(0.f), z(0.f) {}
	p3d(const apt_xyz _x, const apt_xyz _y, const apt_xyz _z) :
		x(_x), y(_y), z(_z) {}
};

ostream& operator<<(ostream& in, p3d const & val);


struct p3d64
{
	double x;
	double y;
	double z;
	p3d64() : x(0.f), y(0.f), z(0.f) {}
	p3d64(const double _x, const double _y, const double _z) :
		x(_x), y(_y), z(_z) {}
};


struct p6d64
{
	/*double xmi;
	double xmx;
	double ymi;
	double ymx;
	double zmi;
	double zmx;*/
	apt_real zmi;
	apt_real zmx;
	p6d64() : zmi(F32MX), zmx(F32MI) {}
	p6d64(const apt_real _zmi, const apt_real _zmx) :
		zmi(_zmi), zmx(_zmx) {}
};


struct p3i
{
	int x;
	int y;
	int z;

	p3i() : x(0), y(0), z(0) {}
	p3i(const int _x, const int _y, const int _z) :
		x(_x), y(_y), z(_z) {}
};

ostream& operator<<(ostream& in, p3i const & val);


struct p6i
{
	int xmi;
	int xmx;
	int ymi;
	int ymx;
	int zmi;
	int zmx;
	p6i() : xmi(0), xmx(0), ymi(0), ymx(0), zmi(0), zmx(0) {}
	p6i( const int _xmi, const int _xmx, const int _ymi, const int _ymx, const int _zmi,
		const int _zmx ) : xmi(_xmi), xmx(_xmx), ymi(_ymi), ymx(_ymx), zmi(_zmi), zmx(_zmx) {}
};

ostream& operator<<(ostream& in, p6i const & val);


struct p3ui64
{
	size_t nx;
	size_t ny;
	size_t nz;
	p3ui64() : nx(0), ny(0), nz(0) {}
	p3ui64( const size_t _nx, const size_t _ny, const size_t _nz ) :
		nx(_nx), ny(_ny), nz(_nz) {}
};


struct vxxl
{
	int x;
	int y;
	int z;
	vxxl() : x(0), y(0), z(0) {}
	vxxl(const int _x, const int _y, const int _z) :
		x(_x), y(_y), z(_z) {}
};

ostream& operator<<(ostream& in, vxxl const & val);

struct p1dm1
{
	apt_xyz pos;
	unsigned int m;

	p1dm1() : pos(0.f), m(UNKNOWNTYPE) {}
	p1dm1(const apt_xyz _pos, const unsigned int _m ) :
		pos(_pos), m(_m) {}
};


ostream& operator<<(ostream& in, p1dm1 const & val);


struct p2dm1
{
	apt_xyz x;
	apt_xyz y;
	size_t m;

	p2dm1() : x(0.f), y(0.f), m(-1) {}
	p2dm1(const apt_xyz _x, const apt_xyz _y, const size_t _m ) :
		x(_x), y(_y), m(_m) {}
};

ostream& operator<<(ostream& in, p2dm1 const & val);


struct p3dm1
{
	apt_xyz x;
	apt_xyz y;
	apt_xyz z;
	unsigned int m;

	p3dm1() : x(0.f), y(0.f), z(0.f), m(UNKNOWNTYPE) {}
	p3dm1(const apt_xyz _x, const apt_xyz _y, const apt_xyz _z, const unsigned int _m ) :
		x(_x), y(_y), z(_z), m(_m) {}
};


ostream& operator<<(ostream& in, p3dm1 const & val);


struct p3dm3
{
	//3*8+4+2+2, 32B per ion, if apt_real double
	apt_real x;
	apt_real y;
	apt_real z;
	unsigned int id;
	unsigned short sgn;
	unsigned short iontype;
	p3dm3() : x(0.f), y(0.f), z(0.f), id(UINT32MX), sgn(0), iontype(UNKNOWNTYPE) {}
	p3dm3( const apt_real _x, const apt_real _y, const apt_real _z,
			const unsigned int _id, const unsigned short _sgn, const unsigned short _it ) :
			x(_x), y(_y), z(_z), id(_id), sgn(_sgn), iontype(_it) {}
};


ostream& operator<<(ostream& in, p3dm3 const & val);



struct d3d
{
	apt_real u;
	apt_real v;
	apt_real w;
	d3d() : u(0.f), v(0.f), w(0.f) {}
	d3d( const apt_real _u, const apt_real _v, const apt_real _w ) :
		u(_u), v(_v), w(_w) {}

	apt_real diff_vector_len();
};

ostream& operator<<(ostream& in, d3d const & val);

struct v3d
{
	apt_real u;
	apt_real v;
	apt_real w;
	apt_real SQR_len;
	v3d() : u(0.f), v(0.f), w(0.f), SQR_len(0.f) {}
	v3d( const apt_real _u, const apt_real _v, const apt_real _w ) :
		u(_u), v(_v), w(_w), SQR_len( SQR(_u)+SQR(_v)+SQR(_w) ) {}

	inline apt_real len() const;
};

ostream& operator<<(ostream& in, v3d const & val);


inline bool SortSQRLenAscending( const v3d &aa1, const v3d &aa2)
{
	return aa1.SQR_len < aa2.SQR_len;
}




struct aabb2d
{
	apt_real xmi;
	apt_real xmx;
	apt_real ymi;
	apt_real ymx;
	apt_real xsz;
	apt_real ysz;
	aabb2d() :
		xmi(F32MX), xmx(F32MI), ymi(F32MX), ymx(F32MI), xsz(0.f), ysz(0.f)  {}
	aabb2d(const apt_real _xmi, const apt_real _xmx, const apt_real _ymi, const apt_real _ymx) :
		xmi(_xmi), xmx(_xmx), ymi(_ymi), ymx(_ymx), xsz(_xmx-_xmi), ysz(_ymx-_ymi) {}

	void scale();
	void blowup( const apt_real f );
	apt_real diag();
};

ostream& operator<<(ostream& in, aabb2d const & val);


struct aabb3d
{
	apt_real xmi;
	apt_real xmx;
	apt_real ymi;
	apt_real ymx;
	apt_real zmi;
	apt_real zmx;
	apt_real xsz;
	apt_real ysz;
	apt_real zsz;
	aabb3d() :
		xmi(F32MX), xmx(F32MI), ymi(F32MX), ymx(F32MI), zmi(F32MX), zmx(F32MI), xsz(0.f), ysz(0.f), zsz(0.f)  {}
	aabb3d(const apt_real _xmi, const apt_real _xmx, const apt_real _ymi, const apt_real _ymx, const apt_real _zmi, const apt_real _zmx) :
		xmi(_xmi), xmx(_xmx), ymi(_ymi), ymx(_ymx), zmi(_zmi), zmx(_zmx), xsz(_xmx-_xmi), ysz(_ymx-_ymi), zsz(_zmx-_zmi) {}

	void add_epsilon_guard();
	void scale();
	void blowup( const apt_real f );
	apt_real diag();
	bool inside( const p3dm1 test );
	bool is_inside_box_xy( aabb3d const & reference, apt_real guard );
	p3d center();
};

ostream& operator<<(ostream& in, aabb3d const & val);


struct cuboidgrid3d
{
	//implements metadata to a 3d aggregate of cuboids in the positive octant of Euclidean space
	p6i gridcells;
	p3d binwidths;
	aabb3d ve;				//the cuboid bounding the aggregate

	cuboidgrid3d() : gridcells(p6i()), binwidths(p3d()), ve(aabb3d()) {}

	//void init( const p3i extend, const p3d dims, aabb3d const & roi );
	p3d where( const int ix, const int iy, const int iz );
	//size_t whichcell( const p3d p );
};

ostream& operator<<(ostream& in, cuboidgrid3d const & val);


struct cylinder
{
	aabb3d aabb;
	apt_real H;
	apt_real R;
	cylinder() :
		aabb(aabb3d()), H(0.f), R(0.f) {}
	cylinder(const aabb3d & _box, const apt_real _H, const apt_real _R ) :
		aabb(_box), H(_H), R(_R) {}

	bool inside( const p3dm1 test );
};

ostream& operator<<(ostream& in, cylinder const & val);


struct jobreceipt
{
	double wallclock;
	int rank;
	int thread;
	unsigned int jobid;

	jobreceipt() : wallclock(0.f), rank(MASTER), thread(MASTER), jobid(0) {}
	jobreceipt( const double _wcl, const int _r, const int _t, const unsigned int _jid) :
		wallclock(_wcl),rank(_r), thread(_t), jobid(_jid) {}
};

ostream& operator<<(ostream& in, jobreceipt const & val);



struct nbor
{
	apt_xyz d;
	unsigned int m;
	nbor() : d(RMAX), m(UNKNOWNTYPE) {}
	nbor(const apt_xyz _d, const unsigned int _m) : d(_d), m(_m) {}
};

ostream& operator<<(ostream& in, nbor const & val);

inline bool SortNeighborsForAscDistance(const nbor & a, const nbor & b)
{
	return a.d < b.d;
}



struct tri3d
{
	apt_real x1;
	apt_real y1;
	apt_real z1;

	apt_real x2;
	apt_real y2;
	apt_real z2;

	apt_real x3;
	apt_real y3;
	apt_real z3;
	tri3d() : x1(0.f), y1(0.f), z1(0.f), x2(0.f), y2(0.f), z2(0.f), x3(0.f), y3(0.f), z3(0.f) {}
	tri3d( const apt_real _x1, const apt_real _y1, const apt_real _z1,
		const apt_real _x2, const apt_real _y2, const apt_real _z2,
		const apt_real _x3, const apt_real _y3, const apt_real _z3 ) :
		x1(_x1), y1(_y1), z1(_z1),
		x2(_x2), y2(_y2), z2(_z2),
		x3(_x3), y3(_y3), z3(_z3) {}
	p3d barycenter();
};

ostream& operator<<(ostream& in, tri3d const & val);


struct triref3d
{
	int v1; //##MK::consider changing to size or long int...
	int v2;
	int v3;
	triref3d() : v1(0), v2(0), v3(0) {}
	triref3d( const int _v1, const int _v2, const int _v3 ) :
		v1(_v1), v2(_v2), v3(_v3) {}
};

ostream& operator<<(ostream& in, triref3d const & val);



struct sqb
{
	//MK::add reject if binning is too fine thereby exceeding UINT32
	size_t nx;
	size_t ny;
	size_t nz;

	size_t nxy;
	size_t nxyz;
	
	apt_xyz width;
	aabb3d box;

	sqb() : nx(1), ny(1), nz(1), nxy(1), nxyz(1), width(F32MX), box(aabb3d()) {}
	sqb(const size_t _nx, const size_t _ny, const size_t _nz, const apt_xyz _w, const aabb3d _bx) :
		nx(_nx), ny(_ny), nz(_nz), nxy(_nx*_ny), nxyz(_nx*_ny*_nz), width(_w), box(_bx) {}

	size_t where( const p3dm1 p );
};

ostream& operator<<(ostream& in, sqb const & val);


struct pdist
{
	size_t pid;		//ion point ID
	apt_xyz dist;	//distance value
	pdist() : pid(-1), dist(0.f) {}
	pdist( const size_t _p, const apt_xyz _d) : pid(_p), dist(_d) {}
};


struct t3x1
{
	real_m33 a11;			//a column vector
	real_m33 a21;
	real_m33 a31;
	t3x1() : 	a11(static_cast<real_m33>(0.0)),
				a21(static_cast<real_m33>(0.0)),
				a31(static_cast<real_m33>(0.0)) {} //initialize to neutral column vector
	t3x1(const real_m33 _a11, const real_m33 _a21, const real_m33 _a31) :
				a11(_a11),
				a21(_a21),
				a31(_a31) {}
};

ostream& operator<<(ostream& in, t3x1 const & val);


struct t3x3
{
	real_m33 a11;				//a second order rank tensor with row-column indexing
	real_m33 a12;
	real_m33 a13;
	real_m33 a21;
	real_m33 a22;
	real_m33 a23;
	real_m33 a31;
	real_m33 a32;
	real_m33 a33;
	t3x3() :	a11(1.0), a12(0.0), a13(0.0),
				a21(0.0), a22(1.0), a23(0.0),
				a31(0.0), a32(0.0), a33(1.0) {}	//initialize to identity tensor
	t3x3( const real_m33* matrix3x3 ) :
				a11(matrix3x3[0]), a12(matrix3x3[1]), a13(matrix3x3[2]),
				a21(matrix3x3[3]), a22(matrix3x3[4]), a23(matrix3x3[5]),
				a31(matrix3x3[6]), a32(matrix3x3[7]), a33(matrix3x3[8]) {}
	t3x3(	const real_m33 _a11, const real_m33 _a12, const real_m33 _a13,
			const real_m33 _a21, const real_m33 _a22, const real_m33 _a23,
			const real_m33 _a31, const real_m33 _a32, const real_m33 _a33 ) :
				a11(_a11), a12(_a12), a13(_a13),
				a21(_a21), a22(_a22), a23(_a23),
				a31(_a31), a32(_a32), a33(_a33) {}
	void add( const t3x3 & increase, const real_m33 weight );
	void div( const real_m33 divisor );
};

std::ostream& operator << (std::ostream& in, t3x3 const & val);


struct plane3d
{
	apt_real n0;	//plane unit normal vector ##MK::so far do not assume consistent sign!
	apt_real n1;
	apt_real n2;
	apt_real d;		//signed distance from world origin
	plane3d() : n0(0.f), n1(0.f), n2(0.f), d(0.f) {} //init with invalid plane
	plane3d( const apt_real _n0, const apt_real _n1, const apt_real _n2, const apt_real _d ) :
		n0(_n0), n1(_n1), n2(_n2), d(_d) {}
	void normaldistance_parameterization( t3x1 const & normal, p3d const & ponplane );
};

ostream& operator<<(ostream& in, plane3d const & val);



//crystallite specific data types
//get planned tip geometry
struct geomodel
{
	apt_real crB;		//relative scaling factors
	apt_real crT;
	apt_real chcapB;
	apt_real chcapT;
	apt_real a;
	size_t N;

	apt_real H;			//nanometer real world dimensions
	apt_real rB;
	apt_real rT;
	apt_real hcapB;
	apt_real hcapT;

	aabb3d mybox;

	geomodel() : crB(0.f), crT(0.f), chcapB(0.f), chcapT(0.f), a(0.f), N(0),
			H(0.f), rB(0.f), rT(0.f), hcapB(0.f), hcapT(0.f), mybox(aabb3d()) {}
	geomodel(const apt_real _crb, const apt_real _crt, const apt_real _chcb, const apt_real _chct,
			const apt_real _a, const size_t _N );
	bool is_inside( const p3d p );
};

struct speci
{
	apt_real c;			//target global composition at %
	apt_real m2q;		//dummy mass to charge, ##MK::i know can be multiple because of isotopes here not considered yet...
	unsigned int typid;

	speci() : c(0.f), m2q(0.f), typid(UNKNOWNTYPE) {}
	speci(const apt_real _c, const apt_real _mc, const unsigned int _t) :
		c(_c), m2q(_mc), typid(_t) {}
};

class solutemodel
{
public:
	solutemodel();
	~solutemodel();

	unsigned int get_random_speci_typid();
	apt_real get_random_speci_mq();

	vector<speci> composition;
	mt19937 urng;			//a random number generator to sample from distribution curve types
};


class unitcellaggr
{
public:
	unitcellaggr();
	unitcellaggr(const apt_real _a, aabb3d unitbox, const unsigned int model );
	~unitcellaggr();

	p3d get_atom(const size_t b, const int u, const int v, const int w);

	apt_real a;
	int umin;
	int umax;
	int vmin;
	int vmax;
	int wmin;
	int wmax;

	//orthogonal base vectors
	v3d a1;
	v3d a2;
	v3d a3;

	//base atoms
	vector<p3d> base;
};


struct msphere
{
	p3d center;
	apt_real radius;
	msphere() : center(p3d()), radius(0.0) {}
	msphere( const p3d _c, const apt_real _r) : center(_c), radius(_r) {}

	void get_atoms( vector<pos> & out );
};


class secondphasemodel
{
public:
	secondphasemodel();
	secondphasemodel(const geomodel & geom,
			const size_t N, const apt_real rm, const apt_real rvar );

	void reportParticles( const unsigned int simid, const int rank );
	void reportParticlesVTK( const unsigned int simid, const int rank );


	vector<msphere> particles;
	mt19937 urng;			//a random number generator to sample from distribution curve types

};


struct occupancy
{
	//characterize how a vxlgrid is occupied by certain types
	size_t ntotal;
	size_t nvacuum;
	size_t nsurface;
	size_t ninside;

	size_t nions_surface;
	size_t nions_inside;

	apt_real volume_inside;
	occupancy() : ntotal(0), nvacuum(0), nsurface(0), ninside(0),
			nions_surface(0), nions_inside(0), volume_inside(0.f) {}
	occupancy( const size_t _ntot, const size_t _nvc, const size_t _nsrf, const size_t _nin,
			const size_t _nisrf, const size_t _niin, const apt_real _vol) :
				ntotal(_ntot), nvacuum(_nvc), nsurface(_nsrf), ninside(_nin),
				nions_surface(_nisrf), nions_inside(_niin), volume_inside(_vol) {}
};

ostream& operator<<(ostream& in, occupancy const & val);


struct zlim
{
	apt_real zmin;
	apt_real zmax;
	zlim() : zmin(F32MX), zmax(F32MI) {}
	zlim(const apt_real _zmi, const apt_real _zmx ) : zmin(_zmi), zmax(_zmx) {}
};


struct localmaximum
{
	float position;
	float strength;
	localmaximum() : position(0.f), strength(0.f) {}
	localmaximum(const float _p, const float _s) : position(_p), strength(_s) {}
};

ostream& operator<<(ostream& in, localmaximum const & val);


inline bool SortLocalMaximaForStrength(const localmaximum & a, const localmaximum & b)
{
	return a.strength < b.strength;
}


#endif
