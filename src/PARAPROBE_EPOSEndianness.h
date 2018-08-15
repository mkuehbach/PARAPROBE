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

#ifndef __PARAPROBE_EPOSENDIANNESS_H__
#define __PARAPROBE_EPOSENDIANNESS_H__

//#include "PARAPROBE_HoshenKopelman.h"
#include "PARAPROBE_BSIMD.h"

//MK::smart struct encapsulates detector events and converts from bigendian-encoded IVAS EPOS binary files into C-compliant little endian
//MK::DO NOT CHANGE THE FLOATING POINT VARIABLE TYPE FROM FLOAT TO DOUBLE IN AN ATTEMPT TO ENFORCE CONSISTENCE WITH REAL!
//THEN THE ENDIANNESS TRANSFORMATION WILL VERY LIKELY NOT WORK! instead we accept that the APT delivers us single precision floats which we bring
//into little endian format and ones being stored in a struct epos we can do with them what we want, for instance
//keeping as floats when EMPLOY_SINGLEPRECISION or promote them implicitly to double when utilizing DOUBLEPRECISION


//native reading of pos files
struct pos_be_2_stl_le {
	//the four BigEndian floats of an IVAS POS file
	unsigned int f1; //x
	unsigned int f2; //y
	unsigned int f3; //z
	unsigned int f4; //mq
};


struct pos {
	float x;
	float y;
	float z;
	float mq;

	inline float my_bswap_uint_2_float( unsigned int ui )
	{
		float retVal = 0.f; //keep it deterministic
		char *pVal = (char*) &ui;
		char* pRetVal = (char*) &retVal;
		for (int i=0; i<4; ++i) {
			pRetVal[4-1-i] = pVal[i];
		}
		return retVal;
	}

	pos() : x(0.f),y(0.f),z(0.f), mq(0.f) {}

	pos(const struct pos_be_2_stl_le in)
	{ //upon storing the detector event within the PARAPROBE data structure we pass the bigendian values and change endianness on the fly
		x = my_bswap_uint_2_float(in.f1);
		y = my_bswap_uint_2_float(in.f2);
		z = my_bswap_uint_2_float(in.f3);
		mq = my_bswap_uint_2_float(in.f4);
	}

	pos(const float _x, const float _y, const float _z, const float _mq) :
		x(_x), y(_y), z(_z), mq(_mq) {}
};


//native reading of epos files
struct epos_be_2_stl_le {
	//the nine IEEE754 BigEndian floats of an IVAS EPOS file
	unsigned int f1; //x;
	unsigned int f2; //y;
	unsigned int f3; //z;

	unsigned int f4; //mq;
	unsigned int f5; //tof;
	unsigned int f6; //vdc;

	unsigned int f7; //vpu;
	unsigned int f8; //detx;
	unsigned int f9; //dety;

	//the two preceeding 32bit BigEndian ints
	int i1;
	int i2;
};


struct epos {
	float x;
	float y;
	float z;

	float mq;
	float tof;
	float vdc;

	float vpu;
	float detx;
	float dety;

	int i1;
	int i2;

	inline float my_bswap_uint_2_float( unsigned int ui ) //;
	//inline float epos::my_bswap_uint_2_float( unsigned int ui )
	{
		float retVal = 0.f; //keep it deterministic
		char *pVal = (char*) &ui;
		char* pRetVal = (char*) &retVal;
		for (int i=0; i<4; ++i) {
			pRetVal[4-1-i] = pVal[i];
		}
		return retVal;
	}

	inline int32_t my_bswap_int32( int32_t val )
	{
		val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF );
		return (val << 16) | ((val >> 16) & 0xFFFF);
	}

	epos() : x(0.f),y(0.f),z(0.f),
				mq(0.f), tof(0.f), vdc(0.f),
					vpu(0.f), detx(0.f), dety(0.f),
						i1(0), i2(0) {}

	epos(const struct epos_be_2_stl_le in)
	{ //upon storing the detector event within the PARAPROBE data structure we pass the bigendian values and change endianness
		x = my_bswap_uint_2_float(in.f1);
		y = my_bswap_uint_2_float(in.f2);
		z = my_bswap_uint_2_float(in.f3);
		mq = my_bswap_uint_2_float(in.f4);

		tof = my_bswap_uint_2_float(in.f5);
		vdc = my_bswap_uint_2_float(in.f6);
		vpu = my_bswap_uint_2_float(in.f7);
		detx = my_bswap_uint_2_float(in.f8);
		dety = my_bswap_uint_2_float(in.f9);

		i1 = my_bswap_int32(in.i1);
		i2 = my_bswap_int32(in.i2);
	}
};

ostream& operator<<(ostream& in, epos const & val);

ostream& operator<<(ostream& in, pos const & val);

#endif
