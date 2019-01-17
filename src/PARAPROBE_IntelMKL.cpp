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


#include "PARAPROBE_IntelMKL.h"


int eig(t3x3 const & in, t3x1 & e_real, t3x1 & e_img, t3x3 & evr_real, t3x3 & evr_img)
{
	//eigenvalues and right handed matrix of eigenvectors
	e_real = t3x1();
	e_img = t3x1();
	evr_real = t3x3();
	evr_img = t3x3();

	MKL_INT n = 3, lda = 3, ldvl = 3, ldvr = 3, info;
	/* Local arrays */
	double wr[3], wi[3], vl[3*3], vr[3*3];
	double a[3*3];
	a[0]= static_cast<double>(in.a11);	a[1]= static_cast<double>(in.a12); 	a[2] = static_cast<double>(in.a13);
	a[3]= static_cast<double>(in.a21);	a[4]= static_cast<double>(in.a22); 	a[5] = static_cast<double>(in.a23);
	a[6]= static_cast<double>(in.a31);	a[7]= static_cast<double>(in.a32); 	a[8] = static_cast<double>(in.a33);

	/*
	//##MK::DEBUG
	cout << "in original" << endl;
	cout << a[0] << "," << a[1] << "," << a[2] << endl;
	cout << a[3] << "," << a[4] << "," << a[5] << endl;
	cout << a[6] << "," << a[7] << "," << a[8] << endl;
	cout << "LAPACKE_dgeev (row-major, high-level) Example Program Results" << endl;
	*/

	//solve eigen problem
	info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'N', 'V', n, a, lda, wr, wi, vl, ldvl, vr, ldvr ); //do not compute left eigenvectors

	if( info == 0 ) { //converged
		e_real.a11 = static_cast<real_m33>(wr[0]);	e_img.a11 = static_cast<real_m33>(wi[0]);//column vector
		e_real.a21 = static_cast<real_m33>(wr[1]);	e_img.a21 = static_cast<real_m33>(wi[1]);
		e_real.a31 = static_cast<real_m33>(wr[2]);	e_img.a31 = static_cast<real_m33>(wi[2]);

		/*
		//##MK::DEBUG
		cout << "three complex eigenvalues" << endl;
		cout << wr[0] << " / " << wi[0] << endl;
		cout << wr[1] << " / " << wi[1] << endl;
		cout << wr[2] << " / " << wi[2] << endl;
		*/

		//unroll for loop to get eigenvalues for (MKL_INT i = 0; i < n; ++i) { //only possible either all eigenvector real then three reals or one real and one complex conjugate pair
		MKL_INT	i = 0;
		MKL_INT	j = 0;
		if ( wi[j] == static_cast<double>(0.0) ) { //##MK::numerically not stable, first eigenvector real?
			evr_real.a11 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a11 = static_cast<real_m33>(0.0);		j++;
			if ( wi[j] == static_cast<double>(0.0) ) { //probe if next eigenvalue also real... , if not it is the complex conjugate pair
				evr_real.a12 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a12 = static_cast<real_m33>(0.0);		j++;
				if ( wi[j] == static_cast<double>(0.0) ) {
					evr_real.a13 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a13 = static_cast<real_m33>(0.0);
				}
				else {
					return MYIMKL_EIG_INCONSISTENT; 	//cout << "After two real, eigenvalue inconsistent" << endl;
				}
			}
			else { //...if not it is the complex conjugate pair
				evr_real.a12 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a12 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
				evr_real.a13 = static_cast<real_m33>(vr[i*ldvr+j]); 	evr_img.a13 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); 	j+=2;
			}
		}
		else { //first eigenvalue is complex therefore first two eigenvector are a conjugate pair, followed by a real valued
			evr_real.a11 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a11 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
			evr_real.a12 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a12 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); 	j+=2;
			if ( wi[j] == static_cast<double>(0.0) ) {
				evr_real.a13 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a13 = static_cast<real_m33>(0.0);			j++;
			}
			else {
				return MYIMKL_EIG_INCONSISTENT; //cout << "After conjugate pair, third eigenvalue inconsistent" << endl;
			}
		}
		i = 1;
		j = 0;
		if ( wi[j] == static_cast<double>(0.0) ) { //first eigenvalue real?
			evr_real.a21 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a21 = static_cast<real_m33>(0.0);		j++;
			if ( wi[j] == static_cast<double>(0.0) ) { //probe if next eigenvalue also real... , if not it is the complex conjugate pair
				evr_real.a22 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a22 = static_cast<real_m33>(0.0);		j++;
				if ( wi[j] == static_cast<double>(0.0) ) {
					evr_real.a23 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a23 = static_cast<real_m33>(0.0);	//j++;
				}
				else {
					return MYIMKL_EIG_INCONSISTENT; //cout << "After two real, eigenvalue inconsistent" << endl;
				}
			}
			else { //...if not it is the complex conjugate pair
				evr_real.a22 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a22 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
				evr_real.a23 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a23 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
			}
		}
		else { //first eigenvalue is complex therefore first two eigenvector are a conjugate pair, followed by a real valued
			evr_real.a21 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a21 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
			evr_real.a22 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a22 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
			if ( wi[j] == static_cast<double>(0.0) ) {
				evr_real.a23 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a23 = static_cast<real_m33>(0.0);			j++;
			}
			else {
				return MYIMKL_EIG_INCONSISTENT; //cout << "After conjugate pair, third eigenvalue inconsistent" << endl;
			}
		}
		i = 2;
		j = 0;
		if ( wi[j] == static_cast<double>(0.0) ) { //first eigenvalue real?
			evr_real.a31 = static_cast<real_m33>(vr[i*ldvr+j]);		evr_img.a31 = static_cast<real_m33>(0.0);		j++;
			if ( wi[j] == static_cast<double>(0.0) ) { //probe if next eigenvalue also real... , if not it is the complex conjugate pair
				evr_real.a32 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a32 = static_cast<real_m33>(0.0);		j++;
				if ( wi[j] == static_cast<double>(0.0) ) {
					evr_real.a33 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a33 = static_cast<real_m33>(0.0);	//j++;
				}
				else {
					return MYIMKL_EIG_INCONSISTENT; //cout << "After two real, eigenvalue inconsistent" << endl;
				}
			}
			else { //...if not it is the complex conjugate pair
				evr_real.a32 = static_cast<real_m33>(vr[i*ldvr+j]);	evr_img.a32 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
				evr_real.a33 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a33 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
			}
		}
		else { //first eigenvalue is complex therefore first two eigenvector are a conjugate pair, followed by a real valued
			evr_real.a31 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a31 = static_cast<real_m33>(+vr[i*ldvr+(j+1)]);
			evr_real.a32 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a32 = static_cast<real_m33>(-vr[i*ldvr+(j+1)]); j+=2;
			if ( wi[j] == static_cast<double>(0.0) ) {
				evr_real.a33 = static_cast<real_m33>(vr[i*ldvr+j]); evr_img.a33 = static_cast<real_m33>(0.0);			j++;
			}
			else {
				return MYIMKL_EIG_INCONSISTENT; //cout << "After conjugate pair, third eigenvalue inconsistent" << endl;
			}
		}

		//##MK::DEBUG
		/*
		cout << "three right-handed eigenvector" << endl;
		//##order of these?
		cout << "first complex column vector" << endl;
		cout << evr_real.a11 << " / " << evr_img.a11 << endl;
		cout << evr_real.a21 << " / " << evr_img.a21 << endl;
		cout << evr_real.a31 << " / " << evr_img.a31 << endl;
		cout << "second complex column vector" << endl;
		cout << evr_real.a12 << " / " << evr_img.a12 << endl;
		cout << evr_real.a22 << " / " << evr_img.a22 << endl;
		cout << evr_real.a32 << " / " << evr_img.a32 << endl;
		cout << "third complex column vector" << endl;
		cout << evr_real.a13 << " / " << evr_img.a13 << endl;
		cout << evr_real.a23 << " / " << evr_img.a23 << endl;
		cout << evr_real.a33 << " / " << evr_img.a33 << endl;
		*/
		return MYIMKL_EIG_SUCCESS;

	}
	else { //not converged
		cout << "MKL EIG dgeev did not converge" << endl;
		return MYIMKL_EIG_NOTCONVERGED;
	}
}




rfftn::rfftn()
{
	NX = 0;
	dimensions = 0;
	fftTempP = NULL;
}


rfftn::~rfftn()
{
	//MK::do not delete fftTempP, only reinterpretation pointer vehicle to the fftTemp vector!
	DftiFreeDescriptor(&m_f_handle);
}


void rfftn::init( const size_t NN, vector<float> const & in )
{
	//single precision real to complex domain fft
	precision = DFTI_SINGLE; //DFTI_DOUBLE;

	//prepare forward fft input buffer with real-valued deformation gradient data
	NX = static_cast<MKL_LONG>(NN);
	dimensions = NX;

	//pull over data to perform an FFT on
	m_input = vector<float>( NX, 0.f); //is equivalent to resize(NX) on empty vector resize sets new elements to default init val
	for( size_t i = 0; i < NN; ++i ) {
		m_input[i] = in[i];
	}

	//define size of OUTPLACE OUTPUT BUFFER reinterpreting a char array
	size_t requiredSize = static_cast<size_t>(NX) * sizeof(MKL_Complex8);
	fftTemp.resize(requiredSize);

	//redefine how to interpret this carrier buffer
	fftTempP = (MKL_Complex8*) &fftTemp[0];
}


void rfftn::fill( vector<float> const & in )
{
	for(size_t i = 0; i < in.size(); ++i) {
		m_input[i] = in[i];
	}
}


void rfftn::forwardFFT()
{
	DftiCreateDescriptor(&m_f_handle, precision, DFTI_REAL, 1, dimensions);
	DftiSetValue(m_f_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	//DftiSetValue(m_f_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
	//DftiSetValue(m_f_handle, DFTI_INPUT_STRIDES, input_strides);
	//DftiSetValue(m_f_handle, DFTI_OUTPUT_STRIDES, output_strides);
	DftiCommitDescriptor(m_f_handle);
	DftiComputeForward(m_f_handle, &m_input[0] , fftTempP);
}


void rfftn::getMagnitude( const size_t NNHalf, vector<float> & out )
{
	for( size_t i = 0; i < NNHalf; ++i) { //##bSIMD
		out[i] = sqrt( SQR(fftTempP[i].real) + SQR(fftTempP[i].imag) );
	}
}
