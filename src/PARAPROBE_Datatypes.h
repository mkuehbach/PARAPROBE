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

	void scale();
	void blowup( const apt_real f );
	apt_real diag();
};

ostream& operator<<(ostream& in, aabb3d const & val);


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


#endif
