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


#ifndef __PARAPROBE_ABINITIO_H__
#define __PARAPROBE_ABINITIO_H__

#include "PARAPROBE_Numerics.h"

//natural constants
#define PI								(3.141592653589793238462643383279502884197169399375105820974)
#define kboltzmann						(1.3806488e-23)		//J/Kmol
#define echarge							(1.602176565e-19)	//Coulomb
#define Navogadro						(6.022140857e+23)	//1/mol
#define RGAS							(8.31446154) 		//(Navogadro)*(kboltzmann)) //J/K/mol

//natural beauty
#define SQR(a)							((a)*(a))
#define CUBE(a)							((a)*(a)*(a))
#define MIN(X,Y)						(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)						(((X) < (Y)) ? (Y) : (X))
#define CLAMP(x,lo,hi)					(MIN(MAX((x), (lo)), (hi)))


//unit conversions
#define CELSIUS2KELVIN(T)				((273.15) + (T))
#define KELVIN2CELSIUS(T)				((T) - (273.15))

#define DEGREE2RADIANT(theta)			((theta)/(180.0)*(PI))
#define RADIANT2DEGREE(rad)				((rad)/(PI)*(180.0))

//scaling conversions
#define NANOMETER2ANGSTROEM(nm)			((10.0)*(nm))
#define ANGSTROEM2NANOMETER(ang)		((ang)/(10.0))
#define METER2NANOMETER(m)				((1.0e9)*(m))
#define NANOMETER2METER(nm)				((1.0e-9)*(nm))


//crystal lattice
#define ATOMS_PER_UNITCELL_FCC			(4)

#define FCC								(0)
#define AL3SC							(1)


#endif
