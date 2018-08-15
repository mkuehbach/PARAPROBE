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

#ifndef __PARAPROBE_TAPSIMHDL_H__
#define __PARAPROBE_TAPSIMHDL_H__

//#include "PARAPROBE_VTKIO.h"
#include "PARAPROBE_XDMF.h"

//##MK::float required as TAPSim v1.0b outputs float
typedef float tapsim_real;

#define TAPSIM_EPSILON		(1.0e-9)

struct tapsim_mvec3d
{
	float x;
	float y;
	float z;
	tapsim_mvec3d() : x(0.f), y(0.f), z(0.f) {}
	tapsim_mvec3d(const float _x, const float _y, const float _z) :
		x(_x), y(_y), z(_z) {}
};

ostream& operator<<(ostream& in, tapsim_mvec3d const & val);

struct tapsim_node
{
	float x;
	float y;
	float z;
	short id;
	short pad;
	tapsim_node() : x(0.f), y(0.f), z(0.f), id(numeric_limits<short>::max()), pad(0) {}
	tapsim_node(const float _x, const float _y, const float _z, const short _id) :
		x(_x), y(_y), z(_z), id(_id), pad(0) {}
};

ostream& operator<<(ostream& in, tapsim_node const & val);

struct tapsim_phaseVector
{
	tapsim_real t; 		//[time] = s
	tapsim_real px;		//[position] = m
	tapsim_real py;
	tapsim_real pz;
	tapsim_real vx;		//[velocity] = m/s
	tapsim_real vy;
	tapsim_real vz;		
	int tetIndex;		//[tetrahedron index] = 1
	//##MK::in case of float there will be no padding

	tapsim_phaseVector() :
		t(0.f), px(0.f), py(0.f), pz(0.f), vx(0.f), vy(0.f), vz(0.f), tetIndex(-1) {}
	tapsim_phaseVector(
		const tapsim_real _t, 
		const tapsim_real _px,  const tapsim_real _py,  const tapsim_real _pz,
		const tapsim_real _vx,  const tapsim_real _vy,  const tapsim_real _vz, const int _tidx) :
			t(_t), px(_px), py(_py), pz(_pz), vx(_vx), vy(_vy), vz(_vz), tetIndex(_tidx) {}
};

ostream& operator<<(ostream& in, tapsim_phaseVector const & val);


class tapsimHdl
{
	//top-level construct implementing the worker instance at process rank within the MPI process level parallelism
	//which can post-process TAPSim simulation output

public:
	tapsimHdl();
	~tapsimHdl();

	string get_filename( const string prefix, const unsigned int id );
	bool read_binary_log( const string fn );
	bool read_binary_pathdata( const string fprefix, const unsigned int logID_incr, const unsigned int logID_e );
	bool read_node_geometry( const string fn );
	void write_vtk_pathdata( const unsigned int ionID_s, const unsigned int ionID_e );
	void write_vtk_nodes();
	void characterize_pathdensity();

	//void set_mpidatatypes( void );
	//inline int get_rank( void ) { return myRank; }
	//inline int get_nranks( void ) { return nRanks; }
	///void set_rank( const int rr ) { myRank = rr; }
	//void set_nranks( const int nnrr ) { nRanks = nnrr; }
	
	vector<tapsim_node> nodes;
	vector<vector<tapsim_phaseVector>*> trajectories;

private:
	//MPI related
	//int myRank;											//my MPI ID in the MPI_COMM_WORLD
	//int nRanks;
};

#endif
