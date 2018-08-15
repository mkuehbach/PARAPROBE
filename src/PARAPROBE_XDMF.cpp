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

#include "PARAPROBE_XDMF.h"

void reconstruction_xdmf( vector<vector<p3d>*> const & ppp,
		vector<vector<unsigned int>*> const & lll, const string xdmf_io_fn )
{
	//writes an XDMF meta data file with which to access heavy data in HDF5 file
	double tic, toc;
	tic = MPI_Wtime();

	if ( ppp.size() != lll.size() ) {
		reporting("Point cloud data and label arrays have dissimilar size!");
		return;
	}

	string mess = "XDMFIO writing meta data file for ion location in reconstruction space to " + xdmf_io_fn;
	reporting( mess );

	toc = MPI_Wtime();
	cout << "XDMF Reconstruction meta data file written in " << (toc - tic) << " seconds" << endl;
}
