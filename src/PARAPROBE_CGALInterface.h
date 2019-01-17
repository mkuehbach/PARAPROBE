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

#ifndef __PARAPROBE_H_CGALINTERFACE_H__
#define __PARAPROBE_H_CGALINTERFACE_H__

#include "PARAPROBE_Settings.h"

//comment this line out to remove all compilation dependencies on and capability of utilizing the CGAL library
#define UTILIZE_CGAL


#ifdef UTILIZE_CGAL

	//required for Convex_hull_3 functionality
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Polyhedron_items_with_id_3.h>
	#include <CGAL/Surface_mesh.h>
	#include <CGAL/convex_hull_3.h>

	typedef CGAL::Exact_predicates_inexact_constructions_kernel			K;
	typedef CGAL::Polyhedron_3<K,CGAL::Polyhedron_items_with_id_3>		Polyhedron_3;
	typedef K::Point_3													Point_3;
	typedef CGAL::Surface_mesh<Point_3>									Surface_mesh;


	//required for Alpha_shape_3 functionality
	//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Delaunay_triangulation_3.h>
	#include <CGAL/Alpha_shape_3.h>

	typedef CGAL::Alpha_shape_vertex_base_3<K>							Vb;
	typedef CGAL::Alpha_shape_cell_base_3<K>							Fb;
	typedef CGAL::Triangulation_data_structure_3<Vb,Fb>					Tds;
	typedef CGAL::Delaunay_triangulation_3<K,Tds,CGAL::Fast_location>	Delaunay;
	typedef CGAL::Alpha_shape_3<Delaunay>								Alpha_shape_3;

	typedef K::Point_3													Point;
	typedef Alpha_shape_3::Alpha_iterator								Alpha_iterator;
	typedef Alpha_shape_3::NT											NT;

#endif


#endif
