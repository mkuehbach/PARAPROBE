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

#include "PARAPROBE_VTKIO.h"

void reconstruction_vtk( vector<vector<p3d>*> const & ppp,
		vector<vector<unsigned int>*> const & lll, runparm const & parms, const string vtk_io_fn )
{
	//MK::write VTK file showing the positions of all ions in reconstructed space
	double tic, toc;
	tic = MPI_Wtime();

	if ( ppp.size() != lll.size() ) {
		reporting("Point cloud data and label arrays have dissimilar size!");
		return;
	}

	string mess = "VTKIO writing ion location in reconstruction space to " + vtk_io_fn;
	reporting( mess );

	//##MK::add VTK write flags

	ofstream vtk;
	vtk.open( vtk_io_fn.c_str() );
	if ( vtk.is_open() == true ) {
		vtk << "# vtk DataFile Version 2.0\n"; //header
		vtk << "PARAPROBE Reconstruction_Eta_ " << parms.eta << "_KF_ " << parms.kf << "_ICF_" << parms.icf << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		//count how many ions in ppp
		size_t pcontrol = 0;
		size_t lcontrol = 0;

		for(size_t b = 0; b < ppp.size(); ++b) {
			if ( ppp.at(b) != NULL )
				pcontrol += ppp.at(b)->size();
		}
		for(size_t b = 0; b < lll.size(); ++b) {
			if ( lll.at(b) != NULL )
				lcontrol += lll.at(b)->size();
		}
		if ( pcontrol != lcontrol ) {
			reporting("Internal content of point cloud data and/or label arrays dissimilar size!");
			vtk.flush();
			vtk.close();
			return;
		}

		size_t nevents = pcontrol;
		vtk << "POINTS " << nevents << " double\n";
		for(size_t b = 0; b < ppp.size(); ++b) {
			if ( ppp.at(b) != NULL ) {
				size_t n = ppp.at(b)->size();
				for(size_t i = 0; i < n; ++i)
					vtk << ppp.at(b)->at(i).x << " " << ppp.at(b)->at(i).y << " " << ppp.at(b)->at(i).z << "\n";
			}
		}
		vtk << "\n";
		vtk << "VERTICES " << nevents << " " << 2*nevents << "\n";
		for ( size_t e = 0; e < nevents; ++e ) {
			vtk << 1 << " " << e << "\n";
		}
		//MK::ranged ion type as field data for coloring in Paraview
		vtk << "POINT_DATA " << nevents << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "IonType 1 " << nevents << " float\n"; //do not make int because numeric_limits<unsigned int>::max() is not readable then
		for(size_t b = 0; b < lll.size(); ++b) {
			if ( lll.at(b) != NULL ) {
				size_t n = lll.at(b)->size();
				for(size_t i = 0; i < n; ++i)
					vtk << lll.at(b)->at(i) << "\n";
			}
		}
		vtk << endl;
		vtk.flush();
		vtk.close();
	}

	toc = MPI_Wtime();
	cout << "VTK Reconstruction ion locations written in " << (toc - tic) << " seconds" << endl;
}



void positions_vtk( vector<p3d> const & ppp, const string vtk_io_fn )
{
	//MK::write VTK file showing positions in reconstructed space
	double tic, toc;
	tic = MPI_Wtime();

	string mess = "VTKIO generic writing ion location in reconstruction space to " + vtk_io_fn;
	reporting( mess );

	//##MK::add VTK write flags

	ofstream vtk;
	vtk.open( vtk_io_fn.c_str() );
	if ( vtk.is_open() == true ) {
		vtk << "# vtk DataFile Version 2.0\n"; //header
		vtk << "PARAPROBE Reconstruction Ion Locations" << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		//count how many ions in ppp
		size_t nevents = ppp.size();
		vtk << "POINTS " << nevents << " double\n";
		for(size_t i = 0; i < ppp.size(); ++i) {
			vtk << ppp.at(i).x << " " << ppp.at(i).y << " " << ppp.at(i).z << "\n";
		}
		vtk << "\n";
		vtk << "VERTICES " << nevents << " " << 2*nevents << "\n";
		for ( size_t e = 0; e < nevents; ++e ) {
			vtk << 1 << " " << e << "\n";
		}
		/*
		//MK::add scalar value per point for coloring in Paraview
		vtk << "POINT_DATA " << nevents << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "Dist2ConvexHull 1 " << nevents << " double\n";
		for ( unsigned int e = 0; e < nevents; ++e ) {
			//vtk << dist->at(e) << " ";
		}
		*/
		vtk << endl;
		vtk.flush();
		vtk.close();
	}

	toc = MPI_Wtime();
	cout << "VTK Reconstruction generic ion locations written in " << (toc - tic) << " seconds" << endl;
}


void triangulation_vtk_naive( vector<tri3d> const & hull, const string vtk_io_fn )
{
	//MK::write VTK triangle hull
	double tic, toc;
	tic = MPI_Wtime();

	string mess = "Tip triangle hull written to VTK file " + vtk_io_fn;
	reporting( mess );

	ofstream vtk;
	vtk.open( vtk_io_fn.c_str() );
	if ( vtk.is_open() == true ) {
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "PARAPROBE TipSurface" << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		size_t nverts = hull.size() * 3; //##MK::naively storing duplicates
		vtk << "POINTS " << nverts << " double\n"; //order preserving
		for ( size_t i = 0; i < hull.size(); ++i ) {
			vtk << hull.at(i).x1 << " " << hull.at(i).y1 << " " << hull.at(i).z1 << "\n";
			vtk << hull.at(i).x2 << " " << hull.at(i).y2 << " " << hull.at(i).z2 << "\n";
			vtk << hull.at(i).x3 << " " << hull.at(i).y3 << " " << hull.at(i).z3 << "\n";
		}
		vtk << "\n";
		//count how many triangles in total
		size_t ntris = hull.size();
		vtk << "POLYGONS " << ntris << " " << 4*ntris << "\n";
		for ( size_t i = 0; i < ntris; ++i ) {
			size_t triid = i*3;
			vtk << 3 << " " << triid+0 << " " << triid+1 << " " << triid+2 << "\n";
		}
		vtk << "\n";
	}
	vtk.flush();
	vtk.close();

	//map and vector temporaries vunique, vuni and tri cleared by default destructor

	toc = MPI_Wtime();
	cout << "Triangle hull written to VTK file " << vtk_io_fn << " in " << (toc-tic) << " seconds" << endl;
}
