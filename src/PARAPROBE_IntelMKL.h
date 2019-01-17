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

#ifndef __PARAPROBE_INTELMKL_H__
#define __PARAPROBE_INTELMKL_H__

//#include "PARAPROBE_Parallelization.h"
#include "PARAPROBE_Profiling.h"



//utilize Intel MKL library
#include "mkl.h"
#include "mkl_dfti.h"

//eigen decomposition
#define MYIMKL_EIG_SUCCESS				1
#define MYIMKL_EIG_INCONSISTENT			-1
#define MYIMKL_EIG_NOTCONVERGED			-2


int eig(t3x3 const &in, t3x1 &e_real, t3x1 &e_img, t3x3 &evr_real, t3x3 &evr_img);


struct vicsfftstatus
{
	bool fft_allc_success;
	bool fft_plac_success;
	bool fft_init_success;
	bool fft_comp_success;
	bool fft_free_success;
	bool pad1;
	bool pad2;
	bool pad3;
	vicsfftstatus() :
		fft_allc_success(false),
		fft_plac_success(false),
		fft_init_success(false),
		fft_comp_success(false),
		fft_free_success(false),
		pad1(false), pad2(false), pad3(false) {}
};


class rfftn
{
public:
	rfftn();
	~rfftn();

	void init( const size_t NN, vector<float> const & in );
	void fill( vector<float> const & in );
	void forwardFFT();
	void getMagnitude( const size_t NNHalf, vector<float> & out );

	vector<char> fftTemp;		//forward transformation results carrier for real and imaginary parts
	vector<float> m_input; 		//input signal, i.e. tensor component

	DFTI_CONFIG_VALUE precision;
	DFTI_DESCRIPTOR_HANDLE m_f_handle;

	MKL_LONG NX;
	MKL_LONG dimensions;
	MKL_Complex8* fftTempP;
};


#endif
