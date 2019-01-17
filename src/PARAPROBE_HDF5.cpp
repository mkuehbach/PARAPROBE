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

//#ifdef UTILIZE_HDF5


ostream& operator<<(ostream& in, h5iometa const & val)
{
	in << val.dsetnm << "__" << val.nr << "__" << val.nc << endl;
	return in;
}

ostream& operator<<(ostream& in, h5offsets const & val)
{
	in << "nr0/nr1__" << val.nr0 << "__" << val.nr1 << "__nc0/nc1__" << val.nc0 << "__" << val.nc1 << "__nrmx/ncrmx__" << val.nrmax << "__" << val.ncmax << endl;
	return in;
}



bool h5offsets::is_within_bounds( const size_t _nr0, const size_t _nr1, const size_t _nc0, const size_t _nc1 )
{
	if ( _nr0 >= this->nr0 && _nr1 <= (this->nr0 + this->nr1) &&
			_nc0 >= this->nc0 && _nc1 <= (this->nc0 + this->nc1) )
		return true;
	else
		return false;
}


/*
size_t h5offsets::interval_length( const size_t _n0, const size_t _n1 )
{
	if ( this->is_within_bounds( _n0, _n1) == true )
		return _n1-_n0;
	else
		return -1;
}
*/


void debug_hdf5( void )
{
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
}


/*
void write_pos_hdf5( vector<vector<pos>*> & in, const string h5_io_fn )
{
//#ifndef UTILIZE_HDF5
	cout << "PARAPROBE was not compiled with HDF5 support" << endl;
	return;
//#else
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
//#endif
}
*/


h5Hdl::h5Hdl()
{
	status = 0;
	fileid = 0;
	groupid = 0;
	mspcid = 0;
	dsetid = 0;
	dspcid = 0;
	cparms = 0;
	fspcid = 0;

	dims[0] = 0; 		dims[1] = 0;
	maxdims[0] = 0; 	maxdims[1] = 0;
	offset[0] = 0; 		offset[1] = 0;
	dimsnew[0] = 0; 	dimsnew[1] = 0;

	nrows = 0;
	ncols = 0;
	BytesPerChunk = 0;
	RowsPerChunk = 0;
	ColsPerChunk = 0;
	RowsColsPerChunk = 0;
	NChunks = 0;
	TotalWriteHere = 0;
	TotalWritePortion = 0;
	TotalWriteAllUsed = 0;
	TotalWriteAllChunk = 0;

	h5resultsfn = "";

	u32le_buf = NULL;
	f64le_buf = NULL;
}


h5Hdl::~h5Hdl()
{
	delete [] u32le_buf; u32le_buf = NULL;
	delete [] f64le_buf; f64le_buf = NULL;
}


void h5Hdl::reinitialize()
{
	status = 0;
	fileid = 0;
	groupid = 0;
	dsetid = 0;
	dspcid = 0;
	cparms = 0;
	fspcid = 0;

	dims[0] = 0; 		dims[1] = 0;
	maxdims[0] = 0; 	maxdims[1] = 0;
	offset[0] = 0; 		offset[1] = 0;
	dimsnew[0] = 0; 	dimsnew[1] = 0;

	nrows = 0;
	ncols = 0;
	BytesPerChunk = 0;
	RowsPerChunk = 0;
	ColsPerChunk = 0;
	RowsColsPerChunk = 0;
	NChunks = 0;
	TotalWriteHere = 0;
	TotalWritePortion = 0;
	TotalWriteAllUsed = 0;
	TotalWriteAllChunk = 0;

	h5resultsfn = "";

	delete [] u32le_buf;
	delete [] f64le_buf;
}


int h5Hdl::create_file( const string h5fn )
{
	h5resultsfn = h5fn;
	fileid = H5Fcreate( h5resultsfn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_paraprobe_results_file( const string h5fn )
{
	//generate PARAPROBE reporting file for general results referring to the reconstruction and surface
	//H5F_ACC_TRUNC truncate content of all files in current working directory with same name as h5fn
	h5resultsfn = h5fn;
	fileid = H5Fcreate( h5resultsfn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	//Ranging is always done, all types and ranges

	//VolumeRecon,xyz,mq
	if ( Settings::IOReconstruction == true || Settings::IOIonTipSurfDists == true ) {
		groupid = H5Gcreate2(fileid, PARAPROBE_VOLRECON, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
	}

	//SurfaceRecon, alphaShape, distances
	if ( Settings::SurfaceMeshingAlgo != E_NOSURFACE ) {
		groupid = H5Gcreate2(fileid, PARAPROBE_SURFRECON, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_SURFRECON_ASHAPE, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_SURFRECON_ASHAPE_HULL, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
	}

	//descriptive statistics
	if ( Settings::SpatialDistributionTask != E_NOSPATSTAT ) {
		groupid = H5Gcreate2(fileid, PARAPROBE_DESCRSTATS, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);

		if ( Settings::SpatStatDoNPCorr == true ) {
			groupid = H5Gcreate2(fileid, PARAPROBE_DESCRSTATS_NCORR, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Gclose(groupid);
		}
	}

	//close file
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_paraprobe_clust_file( const string h5fn )
{
	//generate PARAPROBE reporting file for clustering results
	//H5F_ACC_TRUNC truncate content of all files in current working directory with same name as h5fn
	h5resultsfn = h5fn;
	fileid = H5Fcreate( h5resultsfn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	//clustering
	groupid = H5Gcreate2(fileid, PARAPROBE_CLUST, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(groupid);

	if ( Settings::ClusteringTask == E_MAXSEPARATION ) {
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP_SZOUT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP_SZINN, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP_XYZOUT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP_XYZINN, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP_SZALL_CDF, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
		groupid = H5Gcreate2(fileid, PARAPROBE_CLUST_MAXSEP_SZINN_CDF, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Gclose(groupid);
	}

	//close file
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_paraprobe_cryst_file( const string h5fn )
{
	//generate PARAPROBE reporting file for crystallography results
	h5resultsfn = h5fn;
	fileid = H5Fcreate( h5resultsfn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	//Crystallography, sampling point cloud, elev/azim pair cloud, three strongest signals
	//cloud of points in corresponding spaces, not necessarily regular nd grids!
	groupid = H5Gcreate2(fileid, PARAPROBE_CRYSTALLO, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(groupid);
	groupid = H5Gcreate2(fileid, PARAPROBE_CRYSTALLO_THREESTRONGEST, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(groupid);

	//close file
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


/*
int h5Hdl::create_paraprobe_voronoi_file( const string h5fn )
{
	if ( Settings::VolumeTessellation != E_NOTESS ) {
		//generate PARAPROBE default results reporting file
		//H5F_ACC_TRUNC truncate content of all files in current working directory with same name as h5fn
		h5resultsfn = h5fn;
		fileid = H5Fcreate( h5resultsfn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

		//close file
		status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	return WRAPPED_HDF5_SUCCESS;
}
*/


int h5Hdl::create_group( const string h5fn, const string grpnm )
{
	//##MK::check if file exists
	fileid = H5Fopen(h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	groupid = H5Gcreate2(fileid, grpnm.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(groupid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_contiguous_matrix_u8le( h5iometa const & h5info )
{
	fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hsize_t dims[2] = { h5info.nr, h5info.nc };
	hsize_t maxdims[2] = { h5info.nr, h5info.nc };
	dspcid = H5Screate_simple( 2, dims, maxdims );
	dsetid = H5Dcreate2( fileid, h5info.dsetnm.c_str(), H5T_STD_U8LE, dspcid,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_contiguous_matrix_u32le( h5iometa const & h5info )
{
	fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hsize_t dims[2] = { h5info.nr, h5info.nc };
	hsize_t maxdims[2] = { h5info.nr, h5info.nc };
	dspcid = H5Screate_simple( 2, dims, maxdims );
	dsetid = H5Dcreate2( fileid, h5info.dsetnm.c_str(), H5T_STD_U32LE, dspcid,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_contiguous_matrix_u64le( h5iometa const & h5info )
{
	fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hsize_t dims[2] = { h5info.nr, h5info.nc };
	hsize_t maxdims[2] = { h5info.nr, h5info.nc };
	dspcid = H5Screate_simple( 2, dims, maxdims );
	dsetid = H5Dcreate2( fileid, h5info.dsetnm.c_str(), H5T_STD_U64LE, dspcid,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_contiguous_matrix_i16le( h5iometa const & h5info )
{
	fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hsize_t dims[2] = { h5info.nr, h5info.nc };
	hsize_t maxdims[2] = { h5info.nr, h5info.nc };
	dspcid = H5Screate_simple( 2, dims, maxdims );
	dsetid = H5Dcreate2( fileid, h5info.dsetnm.c_str(), H5T_STD_I16LE, dspcid,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_contiguous_matrix_f32le( h5iometa const & h5info )
{
	fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hsize_t dims[2] = { h5info.nr, h5info.nc };
	hsize_t maxdims[2] = { h5info.nr, h5info.nc };
	dspcid = H5Screate_simple( 2, dims, maxdims );
	dsetid = H5Dcreate2( fileid, h5info.dsetnm.c_str(), H5T_IEEE_F32LE, dspcid,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::create_contiguous_matrix_f64le( h5iometa const & h5info )
{
	fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hsize_t dims[2] = { h5info.nr, h5info.nc };
	hsize_t maxdims[2] = { h5info.nr, h5info.nc };
	dspcid = H5Screate_simple( 2, dims, maxdims );
	dsetid = H5Dcreate2( fileid, h5info.dsetnm.c_str(), H5T_IEEE_F64LE, dspcid,
							H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}



int h5Hdl::create_scalar_u64le( const string h5fn, const string grpnm, const size_t val )
{
	fileid = H5Fopen(h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	int rank = 1;
	hsize_t sdims[1] = { 1};
	dspcid = H5Screate_simple( rank, sdims, NULL);
	dsetid = H5Dcreate2(fileid, grpnm.c_str(), H5T_STD_U64LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	size_t wdata[1] = { val };
	status = H5Dwrite(dsetid, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &wdata );
	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Fclose(fileid);
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::init_chunked_matrix_u32le( const string h5fn, const string dsetnm, const size_t nr, const size_t nc )
{
	//customizes h5Hdl properties for writing a matrix of U32LE nr x nc elements to an existent H5 file
	//we used a chunked data layout in the HDF5 file to allow for inplace compression
	//we caveat is that for arbitrary data the total number of elements nr*nc may not be integer multiple
	//of chunk size, hence we have to fill with buffer
	//##MK::in the future use new partial chunking feature of HDF5 1.10.2 to avoid section of trailing dummies on last chunk
	//we can still though use current concept for visualizing for instance Voronoi cells using XDMF
	//for this we need at least two implicit arrays, one encoding the explicit cell topology
	//another the cell vertex coordinates, however if nr*nc is not necessarily integer multiple of chunk size
	//this is ugly but works as within the XDMF file we can just tell the Dimension of the dataset
	//here we can use the trick that the dimension of the implicit 1d array can be defined shorter
	//than the actual size of the chunked data set

	//we use a several MB write buffer fill before writing the first chunk with dummy values
	//this is useful then fusing successively threadlocal data because threads know their local
	//write entry and exit positions on the global array thereby allowing to handle the complexity
	//of having threadlocal datasets potentially crossing chunk boundaries or generating partially filled chunks
	//which the next thread continues filling and writes to file
	nrows = nr; //real dataset size in number of elements on 1d implicit buffer of unsigned int is nr*nc
	ncols = nc;
	BytesPerChunk = static_cast<size_t>(1*1024*1024); //MK::use default chunk size 1MB, ##MK::performance study
	//##MK::value check needs to be integer multiple of sizeof(unsigned int)*ncols !
	RowsPerChunk = BytesPerChunk / (sizeof(unsigned int) * ncols);
	ColsPerChunk = ncols;
	RowsColsPerChunk = RowsPerChunk*ColsPerChunk;
	NChunks = static_cast<size_t>( ceil(static_cast<double>(nrows*ncols) / static_cast<double>(RowsColsPerChunk)) );

cout << "HDF5 RealDataRows/Cols_BytesPerChunk_Rows/Cols/Rows*ColsPerChunk_NChunksTotal\t\t" << nrows << ";" << ncols << ";" << BytesPerChunk << ";" << RowsPerChunk << ";" << ColsPerChunk << ";" << RowsColsPerChunk << ";" << NChunks << endl;

	//allocate a read/write buffer to pass data to the H5 library call functions
	if ( u32le_buf != NULL ) {
		delete [] u32le_buf;
		u32le_buf = NULL;
	}
	else {
		try {
			u32le_buf = new unsigned int[RowsPerChunk*ColsPerChunk];
		}
		catch (bad_alloc &h5croak) {
			stopping("Allocation error for write buffer in init_chunked_matrix_u32le HDF5");
			return WRAPPED_HDF5_ALLOCERR;
		}
	}

	//open an existent H5 file with read/write access
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	//initialize in this file a chunked data set in an existent group

	//start processing chunks into buffer successively and write to HDF5 file
	dims[0] = RowsPerChunk; 		dims[1] = ColsPerChunk;
	maxdims[0] = H5S_UNLIMITED;		maxdims[1] = H5S_UNLIMITED;
	offset[0] = 0;					offset[1] = 0;
	dimsnew[0] = 0;					dimsnew[1] = ColsPerChunk;

	dspcid = H5Screate_simple(2, dims, maxdims);
	cparms = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk( cparms, 2, dims);

	unsigned int fillval = HDF5_U32LE_DUMMYVALUE;
	status = H5Pset_fill_value(cparms, H5T_STD_U32LE, &fillval );
	dsetid = H5Dcreate2(fileid, dsetnm.c_str(), H5T_STD_U32LE, dspcid, H5P_DEFAULT, cparms, H5P_DEFAULT);

	//fill read/write buffer with dummy values which the preceeding write steps can parse
	for(size_t i = 0; i < RowsColsPerChunk; ++i) {
		u32le_buf[i] = HDF5_U32LE_DUMMYVALUE;
	}
	//and use in combination with the position offsets to always complete first individual chunks,
	//thereafter write them to file and refresh the read/write buffer with dummy values to avoid
	//having to read the file content for any preceeding thread whose data portion interval falls arbitrarily
	//on a chunk boundary
	TotalWriteHere = 0;
	TotalWritePortion = RowsColsPerChunk;
	TotalWriteAllUsed = nrows*ncols;
	TotalWriteAllChunk = NChunks*RowsColsPerChunk; //so many implicitly coded unsigned int to write in total (including potential buffer of trailing dummies)

	//do not close any resource handlers or the file, in what follows we will add data portions per thread
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_chunked_matrix_u32le( const string h5fn, const string dsetnm, vector<unsigned int> const & buffer )
{
	//write a local data portion using opened resources use previous allocated read/write buffer
	//and evaluating Total* interval positioning values to identify whether previous chunk is complete or not
	/*
	cout << "start/end/dims[0]/buffer.size()\t\t" << start << ";" << end << ";" << dims[0] << ";" << buffer.size() << endl;
	//MK::[start, end) mark positions on the initialized current dataset were the buffer data should be written to
	//check that [start,end) is within initialized bounds
	if ( start >= dims[0] || end > dims[0] )
		return WRAPPED_HDF5_OUTOFBOUNDS;
	//check that [start,end) length is the same as the buffer
	if ( (end-start) != buffer.size() )
		return WRAPPED_HDF5_ARGINCONSISTENT;
	*/

	//challenge is that interval [start,end) is not necessarily aligned with the implicitly defined chunk boundaries
	//hence find first chunk which contains array index position start
	size_t LocalCompleted = 0; //how much of feed buffer already processed
	size_t LocalTotal = buffer.size(); //total size of the feed buffer

	size_t BufferWriteHere = 0; //where to start filling in h5 write buffer u32le_buf
	size_t BufferWriteNow = 0; //how many remaining fillable before u32le_buf is full and needs a flush

	do { //loop for as many chunks as fillable with feed buffer content
		//on the fly filling and writing back of completed chunk
		//prior to filling the very first value from feed buffer TotalWriteHere == 0
		//right after filling n-th chunk TotalWriteHere % RowsCols == 0 so buffer was full but already emptied
		//in both these cases the buffer is in a state scientifically completely unfilled
		//MK::so if we feed the results from the threadlocal buffers in the strict ascending order
		//of their threadIDs, i.e. begin with MASTER, the h5 write buffer was not yet flushed to file

		BufferWriteHere = TotalWriteHere % RowsColsPerChunk;
		BufferWriteNow = RowsColsPerChunk;
		//correct write interval end because potential feed buffer not large enough
		//correct write interval end additionally based on whether still sufficient place on buffer available
		if ( (BufferWriteHere+(LocalTotal-LocalCompleted)) < RowsColsPerChunk )
			BufferWriteNow = BufferWriteHere + LocalTotal-LocalCompleted;

cout << "1::BufferWriteHere;Now;LocalCompleted;TotalWriteHere\t\t" << BufferWriteHere << ";" << BufferWriteNow << ";" << LocalCompleted << ";" << TotalWriteHere << endl;
		for( size_t w = BufferWriteHere; w < BufferWriteNow; w++ ) {
			u32le_buf[w] = buffer.at(LocalCompleted);
			LocalCompleted++;
			TotalWriteHere++;
		}
cout << "2::BufferWriteHere;Now;LocalCompleted;TotalWriteHere\t\t" << BufferWriteHere << ";" << BufferWriteNow << ";" << LocalCompleted << ";" << TotalWriteHere << endl;

		//potentially write to H5 file
		if ( TotalWriteHere % RowsColsPerChunk == 0 ) {
cout << "3::Writing a full chunk" << endl;
			dimsnew[0] = dimsnew[0] + RowsPerChunk;
			//second dimension of dimsnew needs to remain as is!
			//the chunk is always written completely, regardless how many values physically relevant or not!
			status = H5Dset_extent(dsetid, dimsnew);
			fspcid = H5Dget_space(dsetid);
			status = H5Sselect_hyperslab(fspcid, H5S_SELECT_SET, offset, NULL, dims, NULL);
			status = H5Dwrite(dsetid, H5T_STD_U32LE, dspcid, fspcid, H5P_DEFAULT, u32le_buf);
			offset[0] = offset[0] + RowsPerChunk;
			//second dimension of offset needs to remain as is!

			//refresh buffer with dummies completely
			for ( size_t i = 0; i < RowsColsPerChunk; ++i )
				u32le_buf[i] = HDF5_U32LE_DUMMYVALUE;
		}

	} while ( LocalCompleted < LocalTotal );
	//processing of feed buffer done

cout << "--->Local done\t\t" << TotalWriteHere << "\t\t" << TotalWriteAllUsed << "\t\t" << TotalWriteAllChunk << endl;

	//dont forget to pass potentially remaining part of the h5 write buffer to file
	if ( TotalWriteHere >= TotalWriteAllUsed ) { //everything which was planned in buffer
		//write h5 buffer out last time if not already happened
		if ( TotalWriteHere % RowsColsPerChunk != 0 ) {
cout << "4::Writing last chunk" << endl;

			dimsnew[0] = dimsnew[0] + RowsPerChunk;
			//second dimension of dimsnew needs to remain as is!
			//the chunk is always written completely, regardless how many values physically relevant or not!
			status = H5Dset_extent(dsetid, dimsnew);
			fspcid = H5Dget_space(dsetid);
			status = H5Sselect_hyperslab(fspcid, H5S_SELECT_SET, offset, NULL, dims, NULL);
			status = H5Dwrite(dsetid, H5T_STD_U32LE, dspcid, fspcid, H5P_DEFAULT, u32le_buf);
			offset[0] = offset[0] + RowsPerChunk;
			//second dimension of offset needs to remain as is!
		}
	}

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::reset_chunked_matrix_u32le_aftercompletion()
{
cout << "Resetting u32le" << endl;

	if ( u32le_buf != NULL ) {
		delete [] u32le_buf;
		u32le_buf = NULL;
	}

	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Sclose(fspcid);
	status = H5Pclose(cparms);
	status = H5Fclose(fileid);

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_contiguous_matrix_u8le_atonce( const string h5fn, const string dsetnm,
				const size_t nr, const size_t nc, vector<unsigned char> const & buffer )
{
	//open existing file and place data in there
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	dims[0] = nr; 		dims[1] = nc;

    dspcid = H5Screate_simple( 2, dims, NULL );
    dsetid = H5Dcreate( fileid, dsetnm.c_str(), H5T_STD_U8LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dsetid, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );
    status = H5Dclose (dsetid);
    status = H5Sclose (dspcid);
    status = H5Fclose (fileid);

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_contiguous_matrix_u32le_atonce( const string h5fn, const string dsetnm,
	const size_t nr, const size_t nc, vector<unsigned int> const & buffer )
{
	//open existing file and place data in there
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	dims[0] = nr; 		dims[1] = nc;

    dspcid = H5Screate_simple( 2, dims, NULL );
    dsetid = H5Dcreate( fileid, dsetnm.c_str(), H5T_STD_U32LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dsetid, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );
    status = H5Dclose (dsetid);
    status = H5Sclose (dspcid);
    status = H5Fclose (fileid);

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_contiguous_matrix_u64le_atonce( const string h5fn, const string dsetnm,
	const size_t nr, const size_t nc, vector<size_t> const & buffer )
{
	//open existing file and place data in there
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	dims[0] = nr; 		dims[1] = nc;

    dspcid = H5Screate_simple( 2, dims, NULL );
    dsetid = H5Dcreate( fileid, dsetnm.c_str(), H5T_STD_U64LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dsetid, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );
    status = H5Dclose (dsetid);
    status = H5Sclose (dspcid);
    status = H5Fclose (fileid);

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_contiguous_matrix_f32le_atonce( const string h5fn, const string dsetnm,
	const size_t nr, const size_t nc, vector<float> const & buffer )
{
	//open existing file and place data in there
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	dims[0] = nr; 		dims[1] = nc;

    dspcid = H5Screate_simple( 2, dims, NULL );
    dsetid = H5Dcreate( fileid, dsetnm.c_str(), H5T_IEEE_F32LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dsetid, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );
    status = H5Dclose (dsetid);
    status = H5Sclose (dspcid);
    status = H5Fclose (fileid);

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_contiguous_matrix_f64le_atonce( const string h5fn, const string dsetnm,
	const size_t nr, const size_t nc, vector<double> const & buffer )
{
	//open existing file and place data in there
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	dims[0] = nr; 		dims[1] = nc;

    dspcid = H5Screate_simple( 2, dims, NULL );
    dsetid = H5Dcreate( fileid, dsetnm.c_str(), H5T_IEEE_F64LE, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite( dsetid, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );
    status = H5Dclose (dsetid);
    status = H5Sclose (dspcid);
    status = H5Fclose (fileid);

	return WRAPPED_HDF5_SUCCESS;
}



int h5Hdl::write_contiguous_matrix_u8le_hyperslab( h5iometa const & h5info,
		h5offsets const & offsetinfo, vector<unsigned char> const & buffer )
{
	//buffer carries implicitly stored matrices row blocks glued together along columns
	//subset dimensions planned in relation to entire dataset
	//global context of the dataset a portion of it should now be written
	hsize_t mxdims[2] = { h5info.nr, h5info.nc };

	//do the write offsets remain within mxdims?
	if ( offsetinfo.nr1 < offsetinfo.nr0 || offsetinfo.nc1 < offsetinfo.nc0 ||
			offsetinfo.nr1 > mxdims[0] || offsetinfo.nc1 > mxdims[1] ) {
		return WRAPPED_HDF5_INCORRECTOFFSETS;
	}

	hsize_t dimssub[2] = { offsetinfo.nr1 - offsetinfo.nr0, offsetinfo.nc1 - offsetinfo.nc0 };


	//bounds check and consistency buffer has the same length than planed
	if ( buffer.size() == dimssub[0] * dimssub[1] ) {
		//define offset
		hsize_t offset[2] = { offsetinfo.nr0, offsetinfo.nc0 };
		hsize_t count[2] = { dimssub[0], dimssub[1] };
		hsize_t stride[2] = { 1, 1};
		hsize_t block[2] = { 1, 1};
		fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		dsetid = H5Dopen2( fileid, h5info.dsetnm.c_str(), H5P_DEFAULT );

		mspcid = H5Screate_simple(2, dimssub, mxdims );
		dspcid = H5Dget_space( dsetid );
		status = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block);
		status = H5Dwrite( dsetid, H5T_STD_U8LE, mspcid, dspcid, H5P_DEFAULT, buffer.data() );

	    status = H5Sclose(mspcid);
	    status = H5Sclose(dspcid);
	    status = H5Dclose(dsetid);
	    status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	else {
		return WRAPPED_HDF5_INCORRECTDIMENSIONS;
	}
}


int h5Hdl::write_contiguous_matrix_u32le_hyperslab( h5iometa const & h5info,
		h5offsets const & offsetinfo, vector<unsigned int> const & buffer )
{
	//buffer carries implicitly stored matrices row blocks glued together along columns
	//subset dimensions planned in relation to entire dataset
	//global context of the dataset a portion of it should now be written
	hsize_t mxdims[2] = { h5info.nr, h5info.nc };

	//do the write offsets remain within mxdims?
	if ( offsetinfo.nr1 < offsetinfo.nr0 || offsetinfo.nc1 < offsetinfo.nc0 ||
			offsetinfo.nr1 > mxdims[0] || offsetinfo.nc1 > mxdims[1] ) {
		return WRAPPED_HDF5_INCORRECTOFFSETS;
	}

	hsize_t dimssub[2] = { offsetinfo.nr1 - offsetinfo.nr0, offsetinfo.nc1 - offsetinfo.nc0 };


	//bounds check and consistency buffer has the same length than planed
	if ( buffer.size() == dimssub[0] * dimssub[1] ) {
		//define offset
		hsize_t offset[2] = { offsetinfo.nr0, offsetinfo.nc0 };
		hsize_t count[2] = { dimssub[0], dimssub[1] };
		hsize_t stride[2] = { 1, 1};
		hsize_t block[2] = { 1, 1};
		fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		dsetid = H5Dopen2( fileid, h5info.dsetnm.c_str(), H5P_DEFAULT );

		mspcid = H5Screate_simple(2, dimssub, mxdims );
		dspcid = H5Dget_space( dsetid );
		status = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block);
		status = H5Dwrite( dsetid, H5T_STD_U32LE, mspcid, dspcid, H5P_DEFAULT, buffer.data() );

	    status = H5Sclose(mspcid);
	    status = H5Sclose(dspcid);
	    status = H5Dclose(dsetid);
	    status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	else {
		return WRAPPED_HDF5_INCORRECTDIMENSIONS;
	}
}


int h5Hdl::write_contiguous_matrix_u64le_hyperslab( h5iometa const & h5info,
			h5offsets const & offsetinfo, vector<size_t> const & buffer )
{
	//buffer carries implicitly stored matrices row blocks glued together along columns
	//subset dimensions planned in relation to entire dataset
	//global context of the dataset a portion of it should now be written
	hsize_t mxdims[2] = { h5info.nr, h5info.nc };

	//do the write offsets remain within mxdims?
	if ( offsetinfo.nr1 < offsetinfo.nr0 || offsetinfo.nc1 < offsetinfo.nc0 ||
			offsetinfo.nr1 > mxdims[0] || offsetinfo.nc1 > mxdims[1] ) {
		return WRAPPED_HDF5_INCORRECTOFFSETS;
	}

	//the ni1 offsets (nr1,nc1) always give exclusive bounds, i.e. index the next array index beyond the last to read/write
	hsize_t dimssub[2] = { offsetinfo.nr1 - offsetinfo.nr0, offsetinfo.nc1 - offsetinfo.nc0 };


	//bounds check and consistency buffer has the same length than planed
	if ( buffer.size() == dimssub[0] * dimssub[1] ) {
		//define offset
		hsize_t offset[2] = { offsetinfo.nr0, offsetinfo.nc0 };
		hsize_t count[2] = { dimssub[0], dimssub[1] };
		hsize_t stride[2] = { 1, 1};
		hsize_t block[2] = { 1, 1};
		fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		dsetid = H5Dopen2( fileid, h5info.dsetnm.c_str(), H5P_DEFAULT );

		mspcid = H5Screate_simple(2, dimssub, mxdims );
		dspcid = H5Dget_space( dsetid );
		status = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block);
		status = H5Dwrite( dsetid, H5T_STD_U64LE, mspcid, dspcid, H5P_DEFAULT, buffer.data() );

	    status = H5Sclose(mspcid);
	    status = H5Sclose(dspcid);
	    status = H5Dclose(dsetid);
	    status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	else {
		return WRAPPED_HDF5_INCORRECTDIMENSIONS;
	}
}


int h5Hdl::write_contiguous_matrix_i16le_hyperslab( h5iometa const & h5info,
						h5offsets const & offsetinfo, vector<short> const & buffer )
{
	//buffer carries implicitly stored matrices row blocks glued together along columns
	//subset dimensions planned in relation to entire dataset
	//global context of the dataset a portion of it should now be written
	hsize_t mxdims[2] = { h5info.nr, h5info.nc };

	//do the write offsets remain within mxdims?
	if ( offsetinfo.nr1 < offsetinfo.nr0 || offsetinfo.nc1 < offsetinfo.nc0 ||
			offsetinfo.nr1 > mxdims[0] || offsetinfo.nc1 > mxdims[1] ) {
		return WRAPPED_HDF5_INCORRECTOFFSETS;
	}

	//the ni1 offsets (nr1,nc1) always give exclusive bounds, i.e. index the next array index beyond the last to read/write
	hsize_t dimssub[2] = { offsetinfo.nr1 - offsetinfo.nr0, offsetinfo.nc1 - offsetinfo.nc0 };


	//bounds check and consistency buffer has the same length than planed
	if ( buffer.size() == dimssub[0] * dimssub[1] ) {
		//define offset
		hsize_t offset[2] = { offsetinfo.nr0, offsetinfo.nc0 };
		hsize_t count[2] = { dimssub[0], dimssub[1] };
		hsize_t stride[2] = { 1, 1};
		hsize_t block[2] = { 1, 1};
		fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		dsetid = H5Dopen2( fileid, h5info.dsetnm.c_str(), H5P_DEFAULT );

		mspcid = H5Screate_simple(2, dimssub, mxdims );
		dspcid = H5Dget_space( dsetid );
		status = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block);
		status = H5Dwrite( dsetid, H5T_STD_I16LE, mspcid, dspcid, H5P_DEFAULT, buffer.data() );

	    status = H5Sclose(mspcid);
	    status = H5Sclose(dspcid);
	    status = H5Dclose(dsetid);
	    status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	else {
		return WRAPPED_HDF5_INCORRECTDIMENSIONS;
	}
}


int h5Hdl::write_contiguous_matrix_f32le_hyperslab( h5iometa const & h5info,
						h5offsets const & offsetinfo, vector<float> const & buffer )
{
	//buffer carries implicitly stored matrices row blocks glued together along columns
	//subset dimensions planned in relation to entire dataset
	//global context of the dataset a portion of it should now be written
	hsize_t mxdims[2] = { h5info.nr, h5info.nc };

	//do the write offsets remain within mxdims?
	if ( offsetinfo.nr1 < offsetinfo.nr0 || offsetinfo.nc1 < offsetinfo.nc0 ||
			offsetinfo.nr1 > mxdims[0] || offsetinfo.nc1 > mxdims[1] ) {
		return WRAPPED_HDF5_INCORRECTOFFSETS;
	}

	//the ni1 offsets (nr1,nc1) always give exclusive bounds, i.e. index the next array index beyond the last to read/write
	hsize_t dimssub[2] = { offsetinfo.nr1 - offsetinfo.nr0, offsetinfo.nc1 - offsetinfo.nc0 };


	//bounds check and consistency buffer has the same length than planed
	if ( buffer.size() == dimssub[0] * dimssub[1] ) {
		//define offset
		hsize_t offset[2] = { offsetinfo.nr0, offsetinfo.nc0 };
		hsize_t count[2] = { dimssub[0], dimssub[1] };
		hsize_t stride[2] = { 1, 1};
		hsize_t block[2] = { 1, 1};
		fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		dsetid = H5Dopen2( fileid, h5info.dsetnm.c_str(), H5P_DEFAULT );

		mspcid = H5Screate_simple(2, dimssub, mxdims );
		dspcid = H5Dget_space( dsetid );
		status = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block);
		status = H5Dwrite( dsetid, H5T_IEEE_F32LE, mspcid, dspcid, H5P_DEFAULT, buffer.data() );

	    status = H5Sclose(mspcid);
	    status = H5Sclose(dspcid);
	    status = H5Dclose(dsetid);
	    status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	else {
		return WRAPPED_HDF5_INCORRECTDIMENSIONS;
	}
}


int h5Hdl::write_contiguous_matrix_f64le_hyperslab( h5iometa const & h5info, h5offsets const & offsetinfo,
		vector<double> const & buffer )
{
	//buffer carries implicitly stored matrices row blocks glued together along columns
	//subset dimensions planned in relation to entire dataset
	//global context of the dataset a portion of it should now be written
	hsize_t mxdims[2] = { h5info.nr, h5info.nc };

	//do the write offsets remain within mxdims?
	if ( offsetinfo.nr1 < offsetinfo.nr0 || offsetinfo.nc1 < offsetinfo.nc0 ||
			offsetinfo.nr1 > mxdims[0] || offsetinfo.nc1 > mxdims[1] ) {
		return WRAPPED_HDF5_INCORRECTOFFSETS;
	}

	hsize_t dimssub[2] = { offsetinfo.nr1 - offsetinfo.nr0, offsetinfo.nc1 - offsetinfo.nc0 };


	//bounds check and consistency buffer has the same length than planed
	if ( buffer.size() == dimssub[0] * dimssub[1] ) {
		//define offset
		hsize_t offset[2] = { offsetinfo.nr0, offsetinfo.nc0 };
		hsize_t count[2] = { dimssub[0], dimssub[1] };
		hsize_t stride[2] = { 1, 1};
		hsize_t block[2] = { 1, 1};
		fileid = H5Fopen( h5resultsfn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
		dsetid = H5Dopen2( fileid, h5info.dsetnm.c_str(), H5P_DEFAULT );

		mspcid = H5Screate_simple(2, dimssub, mxdims );
		dspcid = H5Dget_space( dsetid );
		status = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block);
		status = H5Dwrite( dsetid, H5T_IEEE_F64LE, mspcid, dspcid, H5P_DEFAULT, buffer.data() );

	    status = H5Sclose(mspcid);
	    status = H5Sclose(dspcid);
	    status = H5Dclose(dsetid);
	    status = H5Fclose(fileid);
		return WRAPPED_HDF5_SUCCESS;
	}
	else {
		return WRAPPED_HDF5_INCORRECTDIMENSIONS;
	}
}


int h5Hdl::init_chunked_matrix_f64le( const string h5fn, const string dsetnm, const size_t nr, const size_t nc )
{
	//see additional comments for *_u32le version of these functions...
	nrows = nr; //real dataset size in number of elements on 1d implicit buffer of unsigned int is nr*nc
	ncols = nc;
	BytesPerChunk = static_cast<size_t>(ncols*1024*1024); //MK::use default chunk size 1MB, ##MK::performance study
	//##MK::value check needs to be integer multiple of sizeof(unsigned int)*ncols !
	RowsPerChunk = BytesPerChunk / (sizeof(double) * ncols);
	ColsPerChunk = ncols;
	RowsColsPerChunk = RowsPerChunk*ColsPerChunk;
	NChunks = static_cast<size_t>( ceil(static_cast<double>(nrows*ncols) / static_cast<double>(RowsColsPerChunk)) );

cout << "HDF5 RealDataRows/Cols_BytesPerChunk_Rows/Cols/Rows*ColsPerChunk_NChunksTotal\t\t" << nrows << ";" << ncols << ";" << BytesPerChunk << ";" << RowsPerChunk << ";" << ColsPerChunk << ";" << RowsColsPerChunk << ";" << NChunks << endl;

	//allocate a read/write buffer to pass data to the H5 library call functions
	if ( f64le_buf != NULL ) {
		delete [] f64le_buf;
		f64le_buf = NULL;
	}
	else {
		try {
			f64le_buf = new double[RowsPerChunk*ColsPerChunk];
		}
		catch (bad_alloc &h5croak) {
			stopping("Allocation error for write buffer in init_chunked_matrix_f64le HDF5");
			return WRAPPED_HDF5_ALLOCERR;
		}
	}

	//open an existent H5 file with read/write access
	fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	//initialize in this file a chunked data set in an existent group

	//start processing chunks into buffer successively and write to HDF5 file
	dims[0] = RowsPerChunk; 		dims[1] = ColsPerChunk;
	maxdims[0] = H5S_UNLIMITED;		maxdims[1] = H5S_UNLIMITED;
	offset[0] = 0;					offset[1] = 0;
	dimsnew[0] = 0;					dimsnew[1] = ColsPerChunk;

	dspcid = H5Screate_simple(2, dims, maxdims);
	cparms = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk( cparms, 2, dims);

	double fillval = HDF5_F64LE_DUMMYVALUE;
	status = H5Pset_fill_value(cparms, H5T_IEEE_F64LE, &fillval );
	dsetid = H5Dcreate2(fileid, dsetnm.c_str(), H5T_IEEE_F64LE, dspcid, H5P_DEFAULT, cparms, H5P_DEFAULT);

	//fill read/write buffer with dummy values which the preceeding write steps can parse
	for(size_t i = 0; i < RowsColsPerChunk; ++i) {
		f64le_buf[i] = HDF5_F64LE_DUMMYVALUE;
	}
	//and use in combination with the position offsets to always complete first individual chunks,
	//thereafter write them to file and refresh the read/write buffer with dummy values to avoid
	//having to read the file content for any preceeding thread whose data portion interval falls arbitrarily
	//on a chunk boundary
	TotalWriteHere = 0;
	TotalWritePortion = RowsColsPerChunk;
	TotalWriteAllUsed = nrows*ncols;
	TotalWriteAllChunk = NChunks*RowsColsPerChunk; //so many implicitly coded unsigned int to write in total (including potential buffer of trailing dummies)

	//do not close any resource handlers or the file, in what follows we will add data portions per thread
	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::write_chunked_matrix_f64le( const string h5fn, const string dsetnm, vector<double> const & buffer )
{
	//see additional comments for *_u32le version of this functions...
	size_t LocalCompleted = 0; //how much of feed buffer already processed
	size_t LocalTotal = buffer.size(); //total size of the feed buffer

	size_t BufferWriteHere = 0; //where to start filling in h5 write buffer u32le_buf
	size_t BufferWriteNow = 0; //how many remaining fillable before u32le_buf is full and needs a flush

	do { //loop for as many chunks as fillable with feed buffer content
		BufferWriteHere = TotalWriteHere % RowsColsPerChunk;
		BufferWriteNow = RowsColsPerChunk;

		if ( (BufferWriteHere+(LocalTotal-LocalCompleted)) < RowsColsPerChunk )
			BufferWriteNow = BufferWriteHere + LocalTotal-LocalCompleted;

		//fill buffer
cout << "1::BufferWriteHere;Now;LocalCompleted;TotalWriteHere\t\t" << BufferWriteHere << ";" << BufferWriteNow << ";" << LocalCompleted << ";" << TotalWriteHere << endl;
		for( size_t w = BufferWriteHere; w < BufferWriteNow; w++ ) {
			f64le_buf[w] = buffer.at(LocalCompleted);
			LocalCompleted++;
			TotalWriteHere++;
		}
cout << "2::BufferWriteHere;Now;LocalCompleted;TotalWriteHere\t\t" << BufferWriteHere << ";" << BufferWriteNow << ";" << LocalCompleted << ";" << TotalWriteHere << endl;

		//potentially flush buffer into H5 file
		if ( TotalWriteHere % RowsColsPerChunk == 0 ) {
			//##MK::write#####
cout << "3::Writing a full chunk" << endl;
			dimsnew[0] = dimsnew[0] + RowsPerChunk;
			//second dimension of dimsnew needs to remain as is!
			//the chunk is always written completely, regardless how many values physically relevant or not!
			status = H5Dset_extent(dsetid, dimsnew);
			fspcid = H5Dget_space(dsetid);
			status = H5Sselect_hyperslab(fspcid, H5S_SELECT_SET, offset, NULL, dims, NULL);
			status = H5Dwrite(dsetid, H5T_IEEE_F64LE, dspcid, fspcid, H5P_DEFAULT, f64le_buf);
			offset[0] = offset[0] + RowsPerChunk;
			//second dimension of offset needs to remain as is!

			//refresh buffer with dummies completely
			for ( size_t i = 0; i < RowsColsPerChunk; ++i )
				f64le_buf[i] = HDF5_F64LE_DUMMYVALUE;
		}

	} while ( LocalCompleted < LocalTotal );
	//processing of feed buffer done

cout << "--->Local done\t\t" << TotalWriteHere << "\t\t" << TotalWriteAllUsed << "\t\t" << TotalWriteAllChunk << endl;

	//dont forget to pass potentially remaining part of the h5 write buffer to file
	if ( TotalWriteHere >= TotalWriteAllUsed ) { //everything which was planned in buffer
		//write h5 buffer out last time if not already happened
		if ( TotalWriteHere % RowsColsPerChunk != 0 ) {
cout << "4::Writing last chunk" << endl;

			dimsnew[0] = dimsnew[0] + RowsPerChunk;
			//second dimension of dimsnew needs to remain as is!
			//the chunk is always written completely, regardless how many values physically relevant or not!
			status = H5Dset_extent(dsetid, dimsnew);
			fspcid = H5Dget_space(dsetid);
			status = H5Sselect_hyperslab(fspcid, H5S_SELECT_SET, offset, NULL, dims, NULL);
			status = H5Dwrite(dsetid, H5T_IEEE_F64LE, dspcid, fspcid, H5P_DEFAULT, f64le_buf);
			offset[0] = offset[0] + RowsPerChunk;
			//second dimension of offset needs to remain as is!
		}
	}

	return WRAPPED_HDF5_SUCCESS;
}


int h5Hdl::reset_chunked_matrix_f64le_aftercompletion()
{
cout << "Resetting f64le" << endl;

	if ( f64le_buf != NULL ) {
		delete [] f64le_buf;
		f64le_buf = NULL;
	}

	status = H5Dclose(dsetid);
	status = H5Sclose(dspcid);
	status = H5Sclose(fspcid);
	status = H5Pclose(cparms);
	status = H5Fclose(fileid);

	return WRAPPED_HDF5_SUCCESS;
}


//#endif


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
