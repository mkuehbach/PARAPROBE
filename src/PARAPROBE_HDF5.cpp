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

#include "PARAPROBE_HDF5.h"

#ifdef UTILIZE_HDF5

void debug_hdf5( void )
{
#ifndef UTILIZE_HDF5
	cout << "PARAPROBE was not compiled with HDF5 support" << endl;
	return;
#else
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* Copyright by The HDF Group.                                               *
	* Copyright by the Board of Trustees of the University of Illinois.         *
	* All rights reserved.                                                      *
	*                                                                           *
	* This file is part of HDF5.  The full HDF5 copyright notice, including     *
	* terms governing use, modification, and redistribution, is contained in    *
	* the COPYING file, which can be found at the root of the source code       *
	* distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
	* If you do not have access to either file, you may request a copy from     *
	* help@hdfgroup.org.                                                        *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	/*
	*  This example illustrates how to create a dataset that is a 4 x 6
	*  array.  It is used in the HDF5 Tutorial.
	*/

	//#include "hdf5.h"
	#define FILE "dset.h5"
	//int main() {
	hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
	hsize_t     dims[2];
	herr_t      status;

	/* Create a new file using default properties. */
	file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Create the data space for the dataset. */
	dims[0] = 4;
	dims[1] = 6;
	dataspace_id = H5Screate_simple(2, dims, NULL);

	/* Create the dataset. */
	dataset_id = H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id,
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);

	/* Terminate access to the data space. */
	status = H5Sclose(dataspace_id);

	/* Close the file. */
	status = H5Fclose(file_id);
#endif
}


void write_pos_hdf5( vector<vector<pos>*> & in, const string h5_io_fn )
{
#ifndef UTILIZE_HDF5
	cout << "PARAPROBE was not compiled with HDF5 support" << endl;
	return;
#else
	//##MK::add group names
	//MK::write HDF5 file showing the positions and typ of all ions in reconstructed space
	double tic, toc;
	tic = MPI_Wtime();

	//MK::why HDF5 and chunking? i) enable compression, ii) avoid allocating monolithic block of additional buffer for ppp and lll iii) I/O efficiency
	//take for instance MATLAB loading a file takes intermediately twice the memory, the same with Paraview
	//sure this work nicely for files in the order of 1GB or so but when facing 30GB and more, read will not work because 60GB++ may exceed main memory
	//which is completely unnecessary
	//instead read file successively in packages of several Megabytes

	//##MK::HDF5 chunking requires dataset size to be integer multiple of chunk size therefore round to nearest integer multiple

	size_t Atoms = 0;
	for(size_t b = 0; b< in.size(); ++b) {
		if ( in.at(b) != NULL )
			Atoms += in.at(b)->size();
	}
	if ( Atoms >= UINT32MX ) {
		stopping("Attempting to export more than 4.2 billion atoms is not implemented!"); return;
	}

	size_t BytesPerChunk = 10*1024*1024; //##MK::must be integer multiple of 4 !
	size_t AtomsPerChunk = BytesPerChunk / (4 * 4); //10MB chunk size, four float32 sized 4B each
	size_t NChunks = static_cast<size_t>(ceil(static_cast<double>(Atoms)/static_cast<double>(AtomsPerChunk)));

	cout << "HDF5 Atoms/BytesPerChunk/AtomsPerChunk/NChunks\t\t" << Atoms << ";" << BytesPerChunk << ";" << AtomsPerChunk << ";" << NChunks << endl;

	//allocate a write buffer on stack memory
	float* wbuf = NULL;
	try {
		wbuf = new float[4*AtomsPerChunk];          //##MK::() initializes values to zero
		for(size_t i = 0; i < 4*AtomsPerChunk; ++i) { wbuf[i] = F32MI; }
	}
	catch (bad_alloc &h5croak) {
		stopping("Allocation error for write buffer in write HDF5");
		return;
	}

	//write total number of Atoms in the HDF5 dataset given that the last chunk might contain additional dummy data
	unsigned int wdata[1] = { static_cast<unsigned int>(Atoms) };
	//initialize the HDF5 file, does file with same name exist already overwrite, if not create, in every case open it
	herr_t status;
	hid_t fileid = H5Fcreate( h5_io_fn.c_str() , H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	hsize_t sdims[1] = {1};
	int rank = 1;
	hid_t dspcid = H5Screate_simple(rank, sdims, NULL);
	hid_t dsetid = H5Dcreate2(fileid, "/NumberOfIons", H5T_STD_U32LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dsetid, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata );
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);

	cout << "HDF5 wrote NumberOfIons to file " << wdata[0] << endl;

	//##MK::add versioning within HDF5 file via group ID

	//start processing chunks into buffer successively and write to HDF5 file
	hsize_t dims[2] = { AtomsPerChunk, 4};
	hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
	hsize_t offset[2] = {0, 0};
	hsize_t dimsnew[2] = {0, 4};

	dspcid = H5Screate_simple(2, dims, maxdims);
	hid_t cparms = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk( cparms, 2, dims);
	float fillval = 0.f;
	status = H5Pset_fill_value(cparms, H5T_IEEE_F32LE, &fillval );
	dsetid = H5Dcreate2(fileid, "/XYZMQ", H5T_IEEE_F32LE, dspcid, H5P_DEFAULT, cparms, H5P_DEFAULT);

	hid_t fspcid;

	//prepare HDF5 file for chunked data set
	size_t ChunkCurrPos = 0;
	size_t ChunkMaxiPos = AtomsPerChunk;
	for(size_t b = 0; b < in.size(); ++b) { //successive walk from buckets and collect result until chunk buffer is full then write, refresh buffer and proceed
		if( in.at(b) != NULL ) {
			vector<pos>* these = in.at(b);
			size_t BucketCurrPos = 0;
			size_t BucketMaxiPos = these->size(); //one last last element

			//cout << b << "\t\t\t[ " << BucketCurrPos << " ; " << BucketMaxiPos << " ) " << endl;
			while ( BucketCurrPos < BucketMaxiPos ) { //bucket can be much larger than chunk cache
				size_t ChunkRemaining = ChunkMaxiPos - ChunkCurrPos;
				size_t BucketRemaining = BucketMaxiPos - BucketCurrPos;
				size_t BucketStopPos = (BucketRemaining <= ChunkRemaining) ? BucketMaxiPos : (BucketCurrPos + ChunkRemaining);

				for(     ; BucketCurrPos < BucketStopPos; BucketCurrPos++ ) { //write content into cache
					wbuf[4*ChunkCurrPos+0] = these->at(BucketCurrPos).x;
					wbuf[4*ChunkCurrPos+1] = these->at(BucketCurrPos).y;
					wbuf[4*ChunkCurrPos+2] = these->at(BucketCurrPos).z;
					wbuf[4*ChunkCurrPos+3] = these->at(BucketCurrPos).mq;
					ChunkCurrPos++;
				}

				if ( ChunkCurrPos < ChunkMaxiPos ) { //chunk buffer not full yet
					//cout << "--->fill--->BucketCurrPos/MaxiPos/ChunkCurrPos/MaxiPos\t\t" << BucketCurrPos << "\t\t" << BucketMaxiPos << "\t\t" << ChunkCurrPos << "\t\t" << ChunkMaxiPos << endl;
					continue;
				}
				else { //chunk buffer full, empty write content to HDF5, ##MK::catch H5 errors
					//cout << "--->dump--->BucketCurrPos/MaxiPos/ChunkCurrPos/MaxiPos\t\t" << BucketCurrPos << "\t\t" << BucketMaxiPos << "\t\t" << ChunkCurrPos << "\t\t" << ChunkMaxiPos << endl;
					dimsnew[0] = dimsnew[0] + AtomsPerChunk;
					status = H5Dset_extent(dsetid, dimsnew);
					fspcid = H5Dget_space(dsetid);
					status = H5Sselect_hyperslab(fspcid, H5S_SELECT_SET, offset, NULL, dims, NULL);
					status = H5Dwrite(dsetid, H5T_IEEE_F32LE, dspcid, fspcid, H5P_DEFAULT, wbuf);
					offset[0] = offset[0] + AtomsPerChunk; //MK::second dimensions remains as is!

					cout << "Chunk intermediate " << ChunkCurrPos << " written to HDF5 file" << endl;

					//reset allocated chunk buffer
					for( size_t i = 0; i < 4*AtomsPerChunk; ++i ) { wbuf[i] = F32MI; }
					ChunkCurrPos = 0;

					//cout << "--->dump--->dimsnew[0]\t\t" << dimsnew[0] << endl;
				}

				//##MK::BucketCurrPos increases so somewhen we will get out of the while loop

			} //process all element from the bucket
		} //...if it contains elements
	} //do so for all buckets

	//cout << "Buckets done now ChunkCurrPos\t\t" << ChunkCurrPos << endl;

	//dont forget to write the last chunk as it might not have been filled completely even when all buffer were processed
	if ( ChunkCurrPos > 0 ) {
		dimsnew[0] = dimsnew[0] + AtomsPerChunk;
		status = H5Dset_extent(dsetid, dimsnew);
		fspcid = H5Dget_space(dsetid);
		status = H5Sselect_hyperslab(fspcid, H5S_SELECT_SET, offset, NULL, dims, NULL);
		status = H5Dwrite(dsetid, H5T_IEEE_F32LE, dspcid, fspcid, H5P_DEFAULT, wbuf);
		offset[0] = offset[0] + AtomsPerChunk;
	}

	//done writing HDF5 file, release resources
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Sclose(fspcid);
	status = H5Pclose(cparms);
	status = H5Fclose(fileid);

	cout << "Chunk last " << ChunkCurrPos << " written to HDF5 file" << endl;

	//deallocate chunk buffer in any caseas no longer needed
	delete [] wbuf; wbuf = NULL;

	string mess = "HDF5IO writing ion location in reconstruction space to " + h5_io_fn;
	reporting( mess );

	toc = MPI_Wtime();
	cout << "HDF5 Reconstruction ion locations written in " << (toc - tic) << " seconds" << endl;
#endif
}

#endif


/*
if ( ppp.size() != lll.size() ) {
	reporting("Point cloud data and label arrays have dissimilar size!");
	return;
}

//get total number of elements to write from buffers ppp and lll
size_t nn = 0;
for( size_t b = 0; b < ppp.size(); ++b ) {
	if ( ppp.at(b) != NULL && lll.at(b) != NULL ) { //point coordinate and label data exist
		if ( ppp.at(b)->size() == lll.at(b)->size() ) //dataset is consistent
			nn += ppp.at(b)->size();
		else {
			reporting("Point cloud data and label arrays have dissimilar size!"); return;
		}
	}
	else {
		reporting("Point cloud data and label arrays have dissimilar size!"); return;
	}
}
*/
