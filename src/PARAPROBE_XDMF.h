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

#ifndef __PARAPROBE_XDMF_H__
#define __PARAPROBE_XDMF_H__

#include "PARAPROBE_HDF5.h"

/*
void reconstruction_xdmf( vector<vector<p3d>*> const & ppp,
		vector<vector<unsigned int>*> const & lll, const string xdmf_io_fn );
*/

#define XDMF_HEADER_LINE1				"<?xml version=\"1.0\" ?>"
#define XDMF_HEADER_LINE2				"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>"
#define XDMF_HEADER_LINE3				"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">"

#define WRAPPED_XDMF_SUCCESS				+1
#define WRAPPED_XDMF_IOFAILED				-1


class xdmfHdl
{
//coordinating instance handling all (sequential) writing of XDMF metafiles detailing HDF5 additional metadata
//for visualization for Paraview or VisIt
public:
	xdmfHdl();
	~xdmfHdl();

	//file generation and closing
	int create_volrecon_file( const string xmlfn, const size_t nions, const string h5ref );
	int create_iondistance_file( const string xmlfn, const size_t nions, const string h5ref );
	int create_tipsurface_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const string h5ref );
	int create_voronoicell_vis_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref );
	int create_voronoicell_vol_file( const string xmlfn, const size_t ncells, const string h5ref );

	int create_voronoicell_debug_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref );

private:
	ofstream xdmfout;
};

#endif
