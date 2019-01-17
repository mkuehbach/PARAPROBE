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

/*
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
*/


xdmfHdl::xdmfHdl()
{
}


xdmfHdl::~xdmfHdl()
{
}


int xdmfHdl::create_volrecon_file( const string xmlfn, const size_t nions, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << nions << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << 3*nions << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_XYZ << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Iontype\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 1\" DataType=\"UInt\" Precision=\"1\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLRECON_IONTYPE_IDS << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"SqrdDist\">" << "\n";
		xdmfout << "	    <DataItem Dimensions=\"" << nions << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_SURFDISTSQR << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
 		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_iondistance_file( const string xmlfn, const size_t nions, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"volrecon\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"PolyVertex\" Dimensions=\"1 1 1\"/>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << nions << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_XYZ << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Distance\">" << "\n";
		xdmfout << "	    <DataItem Dimensions=\"" << nions << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLRECON_SURFDISTSQR << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
 		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_tipsurface_file( const string xmlfn, const size_t topo_nelements,
		const size_t topo_dims, const size_t geom_dims, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"alphashape\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << topo_nelements << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << topo_dims << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_SURFRECON_ASHAPE_HULL_TOPO << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << geom_dims << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" <<  PARAPROBE_SURFRECON_ASHAPE_HULL_GEOM << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_voronoicell_vis_file( const string xmlfn, const size_t topo_nelements,
		const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"voronoi\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << topo_nelements << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << topo_dims << "\" NumberType\"UInt\" Precision=\"8\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_TOPOLOGY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << geom_dims << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_GEOMETRY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}



int xdmfHdl::create_voronoicell_vol_file( const string xmlfn, const size_t ncells, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"voronoicells\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Polyvertex\" NodesPerElement=\"" << ncells << "\">" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ncells << " 3\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_DESCRSTATS_CELLPOS << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
		xdmfout << "      <Attribute Name=\"Volume\" AttributeType=\"Scalar\" Center=\"Node\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ncells << " 1\" DataType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_DESCRSTATS_VOL << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "      <Attribute Name=\"ThreadID\" AttributeType=\"Scalar\" Center=\"Node\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << ncells << "\" NumberType=\"UChar\" Precision=\"1\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_DESCRSTATS_THREADID << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}


int xdmfHdl::create_voronoicell_debug_file( const string xmlfn, const size_t topo_nelements,
			const size_t topo_dims, const size_t geom_dims, const size_t attr_dims, const string h5ref )
{
	xdmfout.open( xmlfn.c_str() );
	if ( xdmfout.is_open() == true ) {
		xdmfout << XDMF_HEADER_LINE1 << "\n";
		xdmfout << XDMF_HEADER_LINE2 << "\n";
		xdmfout << XDMF_HEADER_LINE3 << "\n";

		xdmfout << "  <Domain>" << "\n";
		xdmfout << "    <Grid Name=\"voronoi\" GridType=\"Uniform\">" << "\n";
		xdmfout << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << topo_nelements << "\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << topo_dims << "\" NumberType=\"UInt\" Precision=\"8\" Format=\"HDF\">" << "\n";
		xdmfout << "          " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_TOPOLOGY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Topology>" << "\n";
		xdmfout << "      <Geometry GeometryType=\"XYZ\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << geom_dims << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_GEOMETRY << "\n";
		xdmfout << "        </DataItem>" << "\n";
		xdmfout << "      </Geometry>" << "\n";
#ifndef VALIDZONE_IONS_ONLY
		xdmfout << "      <Attribute Name=\"ThreadID\" AttributeType=\"Scalar\" Center=\"Cell\">" << "\n";
		xdmfout << "        <DataItem Dimensions=\"" << attr_dims << " 1\" DataType=\"Int\" Precision=\"2\" Format=\"HDF\">" << "\n";
		xdmfout << "           " << h5ref << ":" << PARAPROBE_VOLTESS_CELLS_THREADIDATTR << "\n";
		xdmfout << "       </DataItem>" << "\n";
		xdmfout << "      </Attribute>" << "\n";
#endif
		xdmfout << "    </Grid>" << "\n";
		xdmfout << "  </Domain>" << "\n";
		xdmfout << "</Xdmf>" << "\n";

		xdmfout.flush();
		xdmfout.close();

		return WRAPPED_XDMF_SUCCESS;
	}
	else {
		return WRAPPED_XDMF_IOFAILED;
	}
}
