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



#ifndef __PARAPROBE_MPIIOEPOS_H__
#define __PARAPROBE_MPIIOEPOS_H__

/*
#define MPIIO_OFFSET_INCR_FLOAT32		4
#define MPIIO_OFFSET_INCR_INT32			4
*/

/*
#define MPI_READ_BLOCKLENGTH		((1)*(1024)*(1024))
bool solverHdl::load_eposfile( void ) 
{
	//check size of Settings::EPOSFilename via operating system call
	double rtimer1 = MPI_Wtime();

	double filesize[1] = {0.0};
	unsigned long nevents = 0;
	struct stat buf;
	string fn = Settings::EPOSFilename;
	if ( stat( fn.c_str() , &buf ) != -1 ) {
		filesize[0] = buf.st_size;
		nevents = (filesize[0] / (double) sizeof(MPI_EPOS_EventIO)); //filesize is in Byte and of character integer, thus
		if ( nevents * sizeof(MPI_EPOS_EventIO) != ((unsigned long) filesize[0]) ) {
			cout << "ERROR::EPOS File size is not an integer multiple of the elementary datatype MPI_EPOS_EventIO!" << endl;
			return false;
		}
	}
	else { 
		cout << "ERROR::EPOS file " << fn << " not readable" << endl; return false;
	}

	//file ais an as expected an integer multiple of the elementary datatype

	//so read blockwise, via MPI_COMM_SELF, native, byte order
	unsigned long ElementsTotal = nevents;
	unsigned long TotalBlocksToRead = 1; //( filesize[0] / MPI_READ_BLOCKLENGTH ) + 1;
	unsigned long ElementsPerBlock = MPI_READ_BLOCKLENGTH / sizeof(MPI_EPOS_EventIO); 
	if ( ElementsPerBlock >= std::numeric_limits<unsigned int>::max() - 1) {
		cout << "ERROR::The number of elements to read per block is too much, reduce MPI_READ_BLOCKLENGTH!" << endl; return false;
	}
	unsigned long ElementsRead = 0;
	unsigned long ElementsRemaining = 0;
	unsigned int ElementsNow = 0; //must be unsinged int

	MPI_File ioReadFileHdl;
	MPI_Status ioReadFileStatus;
	MPI_File_open(MPI_COMM_SELF, fn.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &ioReadFileHdl);
	MPI_File_set_view(ioReadFileHdl, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

	unsigned char* tc = new unsigned char[4];
	MPI_File_read( ioReadFileHdl, tc, 4, MPI_CHAR, &ioReadFileStatus );
	unsigned char tmpc = tc[0];
	tc[0] = tc[3];
	tc[3] = tmpc;
	tmpc = tc[1];
	tc[1] = tc[2];
	tc[2] = tmpc;
	float be2le = *tc;
	cout << setprecision(32) << be2le << endl;


	MPI_File_set_view(ioReadFileHdl, 0, MPI_BYTE, MPI_BYTE, "external32", MPI_INFO_NULL);
	//MPI_File_set_view(ioReadFileHdl, 0, MPI_BYTE, MPI_BYTE, "external32", MPI_INFO_NULL);
	MPI_File_seek(ioReadFileHdl, 0, MPI_SEEK_SET);

	float* testbuf = NULL;
	testbuf = new float[9];
	MPI_File_read( ioReadFileHdl, testbuf, 9, MPI_FLOAT, &ioReadFileStatus);
	for ( unsigned int i = 0; i < 9; i++ ) {
		floatSwap( &testbuf[i] );
		cout << "\t\t" << setprecision(32) << testbuf[i] << endl;
	}
	int* intbuf = NULL;
	intbuf = new int[2];
	
	MPI_File_read( ioReadFileHdl, intbuf, 2, MPI_INT, &ioReadFileStatus);
	cout << "\t\t" << setprecision(32) << intbuf[0] << endl;
	cout << "\t\t" << setprecision(32) << intbuf[1] << endl;

	delete testbuf;
	delete intbuf;
	
	for ( unsigned long b = 0; b < TotalBlocksToRead; b++ ) {
		ElementsNow = ElementsPerBlock;
		ElementsRemaining = ElementsTotal - ElementsRead;
		if ( ElementsRemaining < ((unsigned long) ElementsNow) ) { //partition remaining part
			if ( ElementsRemaining < ((unsigned long) std::numeric_limits<unsigned int>::max()) ) {
				ElementsNow = ElementsRemaining;
			}
			else {
				cout << "ERROR::Remaining number of elements exceeds capability of unsigned int!" << endl; 
				MPI_File_close(&ioReadFileHdl); return false;
			}
		}

		if ( ElementsNow > 0 ) {
			MPI_EPOS_EventIO* rbuf = NULL;
			try {
				rbuf = new MPI_EPOS_EventIO[ElementsNow];
			}
			catch (std::bad_alloc &exc) { 
				cout << "Allocation error in read_mpiio_gnu_binary_2d!" << endl; MPI_File_close(&ioReadFileHdl); return false;
			}

			MPI_File_read( ioReadFileHdl, rbuf, ElementsNow, MPI_EPOS_EventIO_Type, &ioReadFileStatus);

			//explicit datatransfer
			for ( unsigned int e = 0; e < ElementsNow; ++e ) {
				//rawdata.push_back( epos_event(rbuf[e].x, rbuf[e].y, rbuf[e].z, rbuf[e].m_q, rbuf[e].tof, rbuf[e].vdc, 
					rbuf[e].vpu, rbuf[e].detx, rbuf[e].dety, rbuf[e].nulls, rbuf[e].nat_pulse) );
cout << setprecision(32) << rbuf[e].x << ";" << rbuf[e].y << ";" << rbuf[e].z << ";" << rbuf[e].m_q << ";" << rbuf[e].tof << ";" << rbuf[e].vdc << ";" << rbuf[e].vpu << ";" << rbuf[e].detx << ";" << rbuf[e].dety << ";" << rbuf[e].nulls << ";" << rbuf[e].nat_pulse << endl;
			}

			ElementsRead = ElementsRead + ElementsNow;

//cout << "\t\tBlockID/elementsRead/elementsNow/elementsTotal--time = " << b << "/" << elementsRead << "/" << elementsNow << "/" << elementsTotal << "\t\t\t" << (MPI_Wtime() - rtimer) << "\t\tseconds" << endl;
			delete [] rbuf; rbuf = NULL;
		}
	} //data read

	MPI_File_close(&ioReadFileHdl);

	double rtimer2 = MPI_Wtime();
	cout << "...Worker " << this->get_rank() << " read " << fn.c_str() << " of size " << filesize[0] << " Bytes with " << nevents << " successfully in " << setprecision(6) << (rtimer2 - rtimer1) << " seconds" << endl;
	prof.logev( "ReadingEPOS", (rtimer2 - rtimer1) );

	return true;
}*/

#endif
