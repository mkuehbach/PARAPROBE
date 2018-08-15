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


/*
	Implements functionalities to detect whether the microstructure has a percolating network of RX grains
	utilizes the Hoshen-Kopelman algorithm in 3D and implements a pathcompressed weighted union/find algorithm based
	on Kartik Kukreja's implementation http://kartikkukreja.wordpress.com
	Markus K\"uhbach, m.kuehbach (at) mpie.de

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "PARAPROBE_HoshenKopelman.h"


bool percAnalyzer::initialize( const bool* inputdata, const unsigned int nxx, const unsigned int nyy, const unsigned int nzz )
{
	setGridsize( nxx, nyy, nzz );

	cout << "Preconditioning for elimination of false negatives inside tip volume" << endl; //(bins too small for practical likelihood to contain ions)
	cout << "3D, SC site perc NX/Y/Z = " << NX << ";" << NY << ";" << NZ << endl;
	cout << "\t\tInitializing..." << endl;

	//get memory for the HK analysis and the union/find
	try { idd = new unsigned int[NXYZ]; }
	catch (bad_alloc &hkexc) {
		cout << "ERROR::HoshenKopelman unable to allocate label ID field!" << endl; return false;
	}
	try { finder = new UF(NXYZ); }
	catch (bad_alloc &hkexc) {
		cout << "ERROR::HoshenKopelman unable to allocate finder class instance" << endl; return false;
	}
	//##MK::consider optimize as this is a very large preallocation

	for (unsigned int c = 0; c < NXYZ; ++c ) { 
		//MK::BE CAREFUL HERE WE WOULD LIKE TO DETECT CLUSTER OF percolation FALSE'S indicating vacuum or connected bins inside tip volume
		//thus we invert original concept for RX find percolating RX cluster //id field 0 (false, no ions in bin), 1 (true, ions in bin)
		/*idd[c] = 0;
		if ( inputdata[c] == true ) 
			idd[c] = 1;
		*/
		if ( inputdata[c] == false )
			idd[c] = 1;
		else 
			idd[c] = 0;
	}
	cout << "\t\tHoshenKopelman::LabelingInitialized" << endl;
	return true;
}


void percAnalyzer::hk3d_core_nonperiodic( const unsigned int x, const unsigned int y, const unsigned int z )
{
	if ( idd[x+y*NX+z*NXY] == 1 ) { //only RXed (1) to analyze for percolation, periodic boundaries
		unsigned int front = ( z == 0 ? 0 : idd[x+y*NX+(z-1)*NXY] );
		unsigned int bottom = ( y == 0 ? 0 : idd[x+(y-1)*NX+z*NXY] );
		unsigned int left = ( x == 0 ? 0 : idd[(x-1)+y*NX+z*NXY] );

		bool bl = false;
		bool bb = false;
		bool bf = false;
		if ( left > 0 ) bl = true; //> 0, i.e. cell exists surplus has already (at least initially been labelled) once
		if ( bottom > 0 ) bb = true;
		if ( front > 0 ) bf = true;

		//tip covers most of the containers 
		//##MK::consider in further optimization that the likelihood for these cases is a function of the RX fraction to be able to test less ifs
		//##MK::utilize a checksum like bl*100+bb*10+bf*1 in a switch command maybe more efficient...
		//case 000 all neighbors are not-RXed
		if ( bl == false && bb == false && bf == false ) { //because UF assumes initially already NXYZ disjoint labels!
			idd[x+y*NX+z*NXY] = finder->initialAssgn( x+y*NX+z*NXY ); return;
		}

		if ( bl == true && bb == false && bf == false ) { //100
			idd[x+y*NX+z*NXY] = left; return;
		}

		if ( bl == false && bb == true && bf == false ) { //010
			idd[x+y*NX+z*NXY] = bottom; return;
		}

		if ( bl == false && bb == false && bf == true ) { //001
			idd[x+y*NX+z*NXY] = front; return;
		}

		if ( bl == true && bb == true && bf == false ) { //110
			idd[x+y*NX+z*NXY] = this->finder->merge( left, bottom );
			return;
		}
		if ( bl == true && bb == false && bf == true ) { //101
			idd[x+y*NX+z*NXY] = this->finder->merge( left, front );
			return;
		}
		if ( bl == false && bb == true && bf == true ) { //011
			idd[x+y*NX+z*NXY] = this->finder->merge( bottom, front );
			return;
		}
		if ( bl == true && bb == true && bf == true ) {
			//##MK::check if correct?
			//unsigned int cand1 = this->finder->merge( left, bottom );
			//unsigned int cand2 = this->finder->merge( left, front );
			unsigned int cand3 = this->finder->merge( bottom, front );
			idd[x+y*NX+z*NXY] = cand3; //because call for cand1 and cand2 modify tree representation
		}
	} //next cell
}


bool percAnalyzer::hoshen_kopelman( void )
{
	cout << "\t\tInitial identification of true's cluster..." << endl;
	for( unsigned int zz = 0; zz < NZ; zz++ ) {
		for( unsigned int yy = 0; yy < NY; yy++ ) {
			for( unsigned int xx = 0; xx < NX; xx++ ) {
				hk3d_core_nonperiodic( xx, yy, zz );
			}
		}
	}
	return true;
}

bool percAnalyzer::compactify( void )
{
	cout << "\t\tCompactifying..." << endl;

	map<unsigned int, unsigned int> old2new;
	map<unsigned int, unsigned int>::iterator which;
	unsigned int nnew = 0;
	unsigned int croot = 0;

	for ( unsigned int c = 0; c < NXYZ; c++ ) {
		if ( idd[c] > 0 ) { //MK::leaves in our case bins with ions untouched
			croot = this->finder->find_root( idd[c] );
			which = old2new.find( croot );
			if ( which != old2new.end() ) { //key found
				idd[c] = which->second;
			}
			else {
				nnew++;
				old2new[croot] = nnew;
				idd[c] = nnew;
			}
		}
	}

	results.nClusterTotal = nnew;
	//##MK::change format of verbosing
	cout << "\t\tTotal number of clusters = " << to_string(results.nClusterTotal) << endl;


	if ( Settings::IORAWHKClusterID == true ) {
		//##MK::BEGIN DEBUG output binary container with labels
		double tic, toc;
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".AdvPruningCmptfyCID." + to_string(NX) + "." + to_string(NY) + "." + to_string(NZ) + ".raw";
		cout << "Writing MPI I/O " << fn << " HK cluster id" << endl;

		MPI_File msFileHdl;
		MPI_Status msFileStatus;
		//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
		MPI_File_open( MPI_COMM_SELF, fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl );
		long long totalOffset = 0;
		MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

		unsigned int* rawdata = NULL;
		try { rawdata = new unsigned int[NXY]; }
		catch (bad_alloc &hkexc) {
			stopping("Unable to allocate memory for writing labels!");
			return false;
		}

		for( unsigned int z = 0; z < NZ; ++z ) {
			for ( unsigned int i = 0; i < NXY; ++i ) { rawdata[i] = 0; } //debug label
			for ( unsigned int y = 0; y < NY; ++y ) {
				for ( unsigned int x = 0; x < NX; ++x ) {
					rawdata[x+y*NX] = idd[x+y*NX+z*NXY];
				}
			}
			//xy layer at once
			MPI_File_write(msFileHdl, rawdata, NXY, MPI_UNSIGNED, &msFileStatus); //implicit advancement of fp
		} //next region z with regions on stacked top of one another in y
		//##MK::END DEBUG

		delete [] rawdata; rawdata = NULL;

		MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF

		toc = MPI_Wtime();
		cout << "AdvPruningHKClusterID written to MPIIO binary file in " << (toc-tic) << " seconds" << endl;
		//##MK::END DEBUG
	}


	return true;
}


bool percAnalyzer::checkLabeling( void )
{
	cout << "\t\tChecking the correct labeling with the compactified IDs..." << endl;
	//bool labelingerror = false;
	unsigned int LE,RI,FR,RE,BO,TO;
	for ( unsigned int z = 0; z < NZ; z++ ) {
		for (unsigned int y = 0; y < NY; y++ ) {
			for ( unsigned int x = 0; x < NX; x++ ) {
				if ( idd[x+y*NX+z*NXY] > 0 ) {
					//original Tobin Fricke i-->y, j --> x working on NX unit cube in R^3
					/*
					LE = ( x == 0 	? id[(N-1)+y*N+z*NN]	: id[(x-1)+y*N+z*NN] );
					RI = ( x == N-1	? id[0+y*N+z*NN] 		: id[(x+1)+y*N+z*NN] );
					FR = ( y == 0	? id[x+(N-1)*N+z*NN]	: id[x+(y-1)*N+z*NN] );
					RE = ( y == N-1	? id[x+0*N+z*NN]		: id[x+(y+1)*N+z*NN] );
					BO = ( z == 0	? id[x+y*N+(N-1)*NN]	: id[x+y*N+(z-1)*NN] );
					TO = ( z == N-1	? id[x+y*N+0*NN]		: id[x+y*N+(z+1)*NN] );
					*/
					//non-periodic domain
					LE = ( x == 0 		? 0	: idd[(x-1)+y*NX+z*NXY] );
					RI = ( x == NX-1	? 0	: idd[(x+1)+y*NX+z*NXY] );
					FR = ( y == 0		? 0	: idd[x+(y-1)*NX+z*NXY] );
					RE = ( y == NY-1	? 0	: idd[x+(y+1)*NX+z*NXY] );
					BO = ( z == 0		? 0	: idd[x+y*NX+(z-1)*NXY] );
					TO = ( z == NZ-1	? 0	: idd[x+y*NX+(z+1)*NXY] );
							
					unsigned int cand = idd[x+y*NX+z*NXY]; //von Neumann nearest neighbors must have the same label if they are not 0!
					//MK::PRODUCTION
					/*
					if ( LE != 0 && LE != cand ) return false; //labelingerror = true; //cout << "LE error " << x << ";" << y << ";" << z << "\t\t" << LE << "\t\t" << idd[x+y*NX+z*NXY] << endl; return false; }
					if ( RI != 0 && RI != cand ) return false; //labelingerror = true; //cout << "RI error " << x << ";" << y << ";" << z << "\t\t" << RI << "\t\t" << idd[x+y*NX+z*NXY] << endl; return false; }
					if ( FR != 0 && FR != cand ) return false; //labelingerror = true; //cout << "FR error " << x << ";" << y << ";" << z << "\t\t" << FR << "\t\t" << idd[x+y*NX+z*NXY] << endl; return false; }
					if ( RE != 0 && RE != cand ) return false; //labelingerror = true; //cout << "RE error " << x << ";" << y << ";" << z << "\t\t" << RE << "\t\t" << idd[x+y*NX+z*NXY] << endl; return false; }
					if ( BO != 0 && BO != cand ) return false; //labelingerror = true; //cout << "BO error " << x << ";" << y << ";" << z << "\t\t" << BO << "\t\t" << idd[x+y*NX+z*NXY] << endl; return false; }
					if ( TO != 0 && TO != cand ) return false; //labelingerror = true; //cout << "TO error " << x << ";" << y << ";" << z << "\t\t" << TO << "\t\t" << idd[x+y*NX+z*NXY] << endl; return false; }
					*/

					//##MK::DEBUG
					if ( LE != 0 && LE != cand ) { cout << "LE error " << x << ";" << y << ";" << z << "\t\t" << LE << "\t\t" << cand << endl; return false; }
					if ( RI != 0 && RI != cand ) { cout << "RI error " << x << ";" << y << ";" << z << "\t\t" << RI << "\t\t" << cand << endl; return false; }
					if ( FR != 0 && FR != cand ) { cout << "FR error " << x << ";" << y << ";" << z << "\t\t" << FR << "\t\t" << cand << endl; return false; }
					if ( RE != 0 && RE != cand ) { cout << "RE error " << x << ";" << y << ";" << z << "\t\t" << RE << "\t\t" << cand << endl; return false; }
					if ( BO != 0 && BO != cand ) { cout << "BO error " << x << ";" << y << ";" << z << "\t\t" << BO << "\t\t" << cand << endl; return false; }
					if ( TO != 0 && TO != cand ) { cout << "TO error " << x << ";" << y << ";" << z << "\t\t" << TO << "\t\t" << cand << endl; return false; }

				} //label consistency for cluster cell x,y,z checked
			}
		}
	}
	//not yet an inconsistency found?
	return true;
}


bool percAnalyzer::determine_clustersize_distr()
{
	//determine cluster size distribution first
	if ( results.nClusterTotal < 1 ) {
		cout << "ERROR::No cluster were detected so nothing to determine a size distribution from!" << endl;
		return false;
	}

	cout << "\t\tCharacterizing cluster size distribution..." << endl;

	//distribution of the cluster sizes
	vector<unsigned int> cdf;
	cdf.assign( results.nClusterTotal, 0 );

	//indices + 1 of cnt are the cluster names, the content cnt[i] their size
	unsigned int cand = 0;
	for ( unsigned int c = 0; c < NXYZ; ++c ) {
		cand = idd[c];
		if ( cand > 0 ) //largest cluster has ID 1
			cdf.at(cand-1) += 1;
	}

	//sort cluster by their size ascendingly and utilize this to find the largest cluster
	sort( cdf.begin(), cdf.end() );

	results.LargestClusterCnt = cdf.at(cdf.size()-1);
	results.LargestClusterID = cdf.at(cdf.size()-1) + 1;

	cout << "\t\tLargest cluster is cluster " << results.LargestClusterID << " with " << results.LargestClusterCnt << " cells." << endl;
	return true;
}


bool percAnalyzer::rebinarize( const unsigned int target, const unsigned int nsz, bool* bitmap )
{
	//get label at voxel with implicit coordinate target
	unsigned int vacuum_cluster_id = idd[target];
	cout << "\t\tRebinarizing against label " << vacuum_cluster_id << endl;

	for ( unsigned int b = 0; b < nsz; ++b ) {
		if ( idd[b] == vacuum_cluster_id ) //so either a voxel representing a bin with ions included or one of the flicker coincidentally too small and therefore voxels empty in ions deeply inside the tip
			bitmap[b] = true;
		else //really the tip
			bitmap[b] = false;
	}

	cout << "\t\tRebinarization successful!" << endl;


	if ( Settings::IORAWHKClusterID == true ) {
		//##MK::BEGIN DEBUG output binary container with labels
		double tic, toc;
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".AdvPruningRebinCID." + to_string(NX) + "." + to_string(NY) + "." + to_string(NZ) + ".raw";
		cout << "Writing MPI I/O " << fn << " HK cluster id" << endl;

		MPI_File msFileHdl;
		MPI_Status msFileStatus;
		//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
		MPI_File_open( MPI_COMM_SELF, fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl );
		long long totalOffset = 0;
		MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

		unsigned int* rawdata = NULL;
		try { rawdata = new unsigned int[NXY]; }
		catch (bad_alloc &hkexc) {
			stopping("Unable to allocate memory for writing labels!");
			return false;
		}

		for( unsigned int z = 0; z < NZ; ++z ) {
			for ( unsigned int i = 0; i < NXY; ++i ) { rawdata[i] = 0; } //debug label
			for ( unsigned int y = 0; y < NY; ++y ) {
				for ( unsigned int x = 0; x < NX; ++x ) {
					rawdata[x+y*NX] = idd[x+y*NX+z*NXY];
				}
			}
			//xy layer at once
			MPI_File_write(msFileHdl, rawdata, NXY, MPI_UNSIGNED, &msFileStatus); //implicit advancement of fp
		} //next region z with regions on stacked top of one another in y
		//##MK::END DEBUG

		delete [] rawdata; rawdata = NULL;

		MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF

		toc = MPI_Wtime();
		cout << "AdvPruningHKClusterID written to MPIIO binary file in " << (toc-tic) << " seconds" << endl;
		//##MK::END DEBUG
	}

	return true;
}
