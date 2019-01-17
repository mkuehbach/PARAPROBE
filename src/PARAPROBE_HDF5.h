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

#ifndef __PARAPROBE_HDFIO_H__
#define __PARAPROBE_HDFIO_H__

#include "PARAPROBE_VTKIO.h"

//comment out if no HDF5 support desired and switch then also EMPLOY_HDFSUPPORT OFF in the CMake file

#define UTILIZE_HDF5


#define PARAPROBE_RANGING						"/Ranging"
#define PARAPROBE_RANGING_IONTYPE_IDS			"/Ranging/IontypeID"
#define PARAPROBE_RANGING_IONTYPE_MQ			"/Ranging/IontypeMQ"

#define PARAPROBE_VOLRECON						"/VolumeRecon"
#define PARAPROBE_VOLRECON_XYZ					"/VolumeRecon/XYZ"
#define PARAPROBE_VOLRECON_MQ					"/VolumeRecon/MQ"
#define PARAPROBE_VOLRECON_IONTYPE_IDS			"/VolumeRecon/IontypeID"
#define PARAPROBE_VOLRECON_TOPO					"/VolumeRecon/Positions"
#define PARAPROBE_VOLRECON_SURFDISTSQR			"/VolumeRecon/SurfaceDistSQR"

#define PARAPROBE_SURFRECON						"/SurfaceRecon"
#define PARAPROBE_SURFRECON_ASHAPE				"/SurfaceRecon/AlphaShape"
#define PARAPROBE_SURFRECON_ASHAPE_HULL			"/SurfaceRecon/AlphaShape/TriangleHull"
#define PARAPROBE_SURFRECON_ASHAPE_INFO			"/SurfaceRecon/AlphaShape/DescrStats"
#define PARAPROBE_SURFRECON_ASHAPE_HULL_TOPO	"/SurfaceRecon/AlphaShape/TriangleHull/XDMFTopologyValues"
#define PARAPROBE_SURFRECON_ASHAPE_HULL_GEOM	"/SurfaceRecon/AlphaShape/TriangleHull/XDMFxyzValues"
#define PARAPROBE_SURFRECON_ASHAPE_ION2DIST		"/SurfaceRecon/AlphaShape/Ion2Distance"

#define PARAPROBE_CRYSTALLO						"/Crystallography"
#define PARAPROBE_CRYSTALLO_MATPOINTXYZ			"/Crystallography/SamplingPointCloudXYZ"
#define PARAPROBE_CRYSTALLO_MATPOINTMETA		"/Crystallography/SamplingPointCloudMeta"
#define PARAPROBE_CRYSTALLO_ELEVAZIMMETA		"/Crystallography/ElevAzimPointCloud"
#define PARAPROBE_CRYSTALLO_THREESTRONGEST		"/Crystallography/ThreeStrongest"

#define PARAPROBE_DESCRSTATS					"/DescrSpatStats"
#define PARAPROBE_DESCRSTATS_RDF				"/DescrSpatStats/RDF"
#define PARAPROBE_DESCRSTATS_RDF_PDF			"/DescrSpatStats/RDF/PDF"
#define PARAPROBE_DESCRSTATS_RDF_CDF			"/DescrSpatStats/RDF/CDF"
#define PARAPROBE_DESCRSTATS_KNN				"/DescrSpatStats/KNN"
#define PARAPROBE_DESCRSTATS_NCORR				"/DescrSpatStats/TwoPointStats"
#define PARAPROBE_DESCRSTATS_NCORR_CELLCENTER	"/DescrSpatStats/TwoPointStats/BinCenterXYZ"
#define PARAPROBE_DESCRSTATS_NCORR_BINNING		"/DescrSpatStats/TwoPointStats/BinNXNYNZ"

#define PARAPROBE_CLUST							"/Clustering"
#define PARAPROBE_CLUST_MAXSEP					"/Clustering/MaximumSeparation"
#define PARAPROBE_CLUST_MAXSEP_SZOUT			"/Clustering/MaximumSeparation/PrecSizeOut"
#define PARAPROBE_CLUST_MAXSEP_SZINN			"/Clustering/MaximumSeparation/PrecSizeInn"
#define PARAPROBE_CLUST_MAXSEP_XYZOUT			"/Clustering/MaximumSeparation/PrecXYZOut"
#define PARAPROBE_CLUST_MAXSEP_XYZINN			"/Clustering/MaximumSeparation/PrecXYZInn"
#define PARAPROBE_CLUST_MAXSEP_SZALL_CDF		"/Clustering/MaximumSeparation/PrecSizeAllCDF"
#define PARAPROBE_CLUST_MAXSEP_SZINN_CDF		"/Clustering/MaximumSeparation/PrecSizeInnCDF"


//##MK::so far Voronoi tessellation is generated in exclusive file as it might get extremely large
#define PARAPROBE_VOLTESS						"/VoronoiTess"
#define PARAPROBE_VOLTESS_DESCRSTATS			"/VoronoiTess/DescrStats"
#define PARAPROBE_VOLTESS_DESCRSTATS_NCELLS		"/VoronoiTess/DescrStats/NumberOfCells"
#define PARAPROBE_VOLTESS_DESCRSTATS_I2CELL		"/VoronoiTess/DescrStats/Ion2CellMappingValues"
#define PARAPROBE_VOLTESS_DESCRSTATS_I2TYPE		"/VoronoiTess/DescrStats/Ion2Iontype"
#define PARAPROBE_VOLTESS_DESCRSTATS_CELLPOS	"/VoronoiTess/DescrStats/CellPositions"
#define PARAPROBE_VOLTESS_DESCRSTATS_VOL		"/VoronoiTess/DescrStats/VolumeValues"
#define PARAPROBE_VOLTESS_DESCRSTATS_THREADID	"/VoronoiTess/DescrStats/ThreadID"
#define PARAPROBE_VOLTESS_DESCRSTATS_NFACES		"/VoronoiTess/DescrStats/NumberOfFacesValues"
#define PARAPROBE_VOLTESS_DESCRSTATS_BND		"/VoronoiTess/DescrStats/BoxContact"
#define PARAPROBE_VOLTESS_DESCRSTATS_NFTOTAL	"/VoronoiTess/DescrStats/NumberOfFacetsTotal"

#define PARAPROBE_VOLTESS_CELLS					"/VoronoiTess/CellGeometry"
#define PARAPROBE_VOLTESS_CELLS_TOPOLOGY		"/VoronoiTess/CellGeometry/XDMFTopologyValues"
#define PARAPROBE_VOLTESS_CELLS_GEOMETRY		"/VoronoiTess/CellGeometry/XDMFxyzValues"
#define PARAPROBE_VOLTESS_CELLS_THREADIDATTR	"/VoronoiTess/CellGeometry/XDMFThreadID"

//#ifdef UTILIZE_HDF5
	#include "hdf5.h"

	void debug_hdf5( void );

	//void reconstruction_read_hdf5( vector<vector<p3d>*> ppp, vector<vector<apt_real>* >)

	/*
	void write_pos_hdf5( vector<vector<pos>*> & in, const string h5_io_fn ); //##MK::add group names
	*/

	//dummy values for filling chunk buffers
	#define HDF5_U32LE_DUMMYVALUE				1 //UINT32MX //##MK::should be 1 if used for topology
	#define HDF5_F64LE_DUMMYVALUE				0.0

	//return codes for wrapped H5 functions
	#define WRAPPED_HDF5_SUCCESS				+1 //MK::following the HDF5 convention that error values are positiv in case of success or negative else
	#define WRAPPED_HDF5_ALLOCERR				-1
	#define WRAPPED_HDF5_OUTOFBOUNDS			-2
	#define WRAPPED_HDF5_ARGINCONSISTENT		-3
	#define WRAPPED_HDF5_EXECUTIONORDERISSUE	-4
	#define WRAPPED_HDF5_INCORRECTLOGIC			-5
	#define WRAPPED_HDF5_INCORRECTDIMENSIONS	-6
	#define WRAPPED_HDF5_INCORRECTOFFSETS		-7

	struct io_bounds
	{
		size_t n;
		size_t s;
		size_t e;
		io_bounds() : n(0), s(0), e(0) {}
		io_bounds( const size_t _n, const size_t _s, const size_t _e ) :
			n(_n), s(_s), e(_e) {}
	};

	struct voro_io_bounds
	{
		//assists coordinates writing of tessellation portions from individual threads
		//by defining the array bounds between which data from the individual threads are placed
		size_t cell_n;
		size_t cell_s; //defining [cell_s, cell_e)
		size_t cell_e;
		size_t topo_n;
		size_t topo_s;
		size_t topo_e;
		size_t geom_n;
		size_t geom_s;
		size_t geom_e;
		voro_io_bounds() : 	cell_n(0), cell_s(0), cell_e(0),
							topo_n(0), topo_s(0), topo_e(0),
							geom_n(0), geom_s(0), geom_e(0) {}
		voro_io_bounds( const size_t _cn, const size_t _cs, const size_t _ce,
						const size_t _tn, const size_t _ts, const size_t _te,
						const size_t _gn, const size_t _gs, const size_t _ge ) :
							cell_n(_cn), cell_s(_cs), cell_e(_ce),
							topo_n(_tn), topo_s(_ts), topo_e(_te),
							geom_n(_gn), geom_s(_gs), geom_e(_ge) {}
	};


	struct voro_io_info
	{
		size_t id;
		size_t nfacets;
		size_t ntopo;
		size_t ngeom;
		voro_io_info() : id(0), nfacets(0), ntopo(0), ngeom(0) {}
		voro_io_info(const size_t _id, const size_t _nf,
				const size_t _nt,	const size_t _ng) :
					id(_id), nfacets(_nf), ntopo(_nt), ngeom(_ng) {}
	};

	struct h5iometa
	{
		string h5fn;
		string dsetnm;
		size_t nr;
		size_t nc;
		h5iometa() : h5fn(""), dsetnm(""), nr(0), nc(0) {}
		h5iometa( const string _dn, const size_t _nr, const size_t _nc ) :
			h5fn(""), dsetnm(_dn), nr(_nr), nc(_nc) {}
		h5iometa( const string _fn, const string _dn, const size_t _nr, const size_t _nc ):
			h5fn(_fn), dsetnm(_dn), nr(_nr), nc(_nc) {}
	};

	ostream& operator<<(ostream& in, h5iometa const & val);


	struct h5dims
	{
		size_t nr;
		size_t nc;
		h5dims() : nr(0), nc(0) {}
		h5dims(const size_t _nr, const size_t _nc) : nr(_nr), nc(_nc) {}
	};


	struct h5offsets
	{
		size_t nr0;			//where to start reading writing off from row ID C indexing, i.e. inclusive bound
		size_t nr1;			//one past where to stop reading writing off from, i.e. exclusive bound
		size_t nc0;			//where to start reading writing off from col ID C indexing
		size_t nc1;
		size_t nrmax;		//maximum dimensions
		size_t ncmax;
		h5offsets() : nr0(-1), nr1(-1), nc0(-1), nc1(-1), nrmax(0), ncmax(0) {}
		h5offsets( const size_t _nr0, const size_t _nr1, const size_t _nc0, const size_t _nc1,
				const size_t _nrmx, const size_t _ncmx ) :
			nr0(_nr0), nr1(_nr1), nc0(_nc0), nc1(_nc1), nrmax(_nrmx), ncmax(_ncmx) {}
		bool is_within_bounds( const size_t _nr0, const size_t _nr1, const size_t _nc0, const size_t _nc1 );
		//size_t interval_length( const size_t _n0, const size_t _n1 );
	};

	ostream& operator<<(ostream& in, h5offsets const & val);


	class h5Hdl
	{
	//coordinating instance handling all (sequential) writing to HDF5 file wrapping HDF5 low-level C library calls
	public:
		h5Hdl();
		~h5Hdl();

		//file generation and closing
		void reinitialize();
		int create_file( const string h5fn );
		int create_paraprobe_results_file( const string h5fn );
		int create_paraprobe_clust_file( const string h5fn );
		int create_paraprobe_cryst_file( const string h5fn );
		int create_paraprobe_voronoi_file( const string h5fn );

		//group generation
		int create_group( const string h5fn, const string grpnm );

		//contiguous dataset generation and fill in for subsequent hyperslab writing
		int create_contiguous_matrix_u8le( h5iometa const & h5info );
		int create_contiguous_matrix_u32le( h5iometa const & h5info );
		int create_contiguous_matrix_u64le( h5iometa const & h5info );
		int create_contiguous_matrix_i16le( h5iometa const & h5info );
		int create_contiguous_matrix_f32le( h5iometa const & h5info );
		int create_contiguous_matrix_f64le( h5iometa const & h5info );

		//scalars
		int create_scalar_u64le( const string h5fn, const string grpnm, const size_t val );

		//nr x nc matrices
		int init_chunked_matrix_u32le( const string h5fn, const string dsetnm, const size_t nr, const size_t nc );
		int write_chunked_matrix_u32le(	const string h5fn, const string dsetnm, vector<unsigned int> const & buffer );
		int reset_chunked_matrix_u32le_aftercompletion();

		/*	//##MK::incremental writing of threadlocal results successively into contiguous dataset does not work
		int init_contiguous_matrix_u32le( const string h5fn, const string dsetnm, const size_t nr, const size_t nc );
		int write_contiguous_matrix_u32le(	const string h5fn, const string dsetnm, vector<unsigned int> const & buffer );
		int reset_contiguous_matrix_u32le_aftercompletion();
		*/

		int write_contiguous_matrix_u8le_atonce( const string h5fn, const string dsetnm,
				const size_t nr, const size_t nc, vector<unsigned char> const & buffer );
		int write_contiguous_matrix_u32le_atonce( const string h5fn, const string dsetnm,
				const size_t nr, const size_t nc, vector<unsigned int> const & buffer );
		int write_contiguous_matrix_u64le_atonce( const string h5fn, const string dsetnm,
						const size_t nr, const size_t nc, vector<size_t> const & buffer );
		int write_contiguous_matrix_f32le_atonce( const string h5fn, const string dsetnm,
						const size_t nr, const size_t nc, vector<float> const & buffer );
		int write_contiguous_matrix_f64le_atonce( const string h5fn, const string dsetnm,
				const size_t nr, const size_t nc, vector<double> const & buffer );

		int write_contiguous_matrix_u8le_hyperslab( h5iometa const & h5info,
				h5offsets const & offsetinfo, vector<unsigned char> const & buffer );
		int write_contiguous_matrix_u32le_hyperslab( h5iometa const & h5info,
				h5offsets const & offsetinfo, vector<unsigned int> const & buffer );
		int write_contiguous_matrix_u64le_hyperslab( h5iometa const & h5info,
				h5offsets const & offsetinfo, vector<size_t> const & buffer );
		int write_contiguous_matrix_i16le_hyperslab( h5iometa const & h5info,
				h5offsets const & offsetinfo, vector<short> const & buffer );
		int write_contiguous_matrix_f32le_hyperslab( h5iometa const & h5info,
				h5offsets const & offsetinfo, vector<float> const & buffer );
		int write_contiguous_matrix_f64le_hyperslab( h5iometa const & h5info,
				h5offsets const & offsetinfo, vector<double> const & buffer );

		int init_chunked_matrix_f64le( const string h5fn, const string dsetnm, const size_t nr, const size_t nc );
		int write_chunked_matrix_f64le( const string h5fn, const string dsetnm, vector<double> const & buffer );
		int reset_chunked_matrix_f64le_aftercompletion();

		//write data
		//we use buffer here that do not cover necessarily the entire dataset to
		//hide the entire writing of data to an HDF5 file from the higher level function call
		//how internally the data stripe [start, end) is stored chunked into only the HDF5 class object takes care of
		int write_chunked_matrix_f64le( const string dsetnm, const size_t nr, const size_t nc,
						const size_t start, const size_t end, vector<p3d> const & buffer );

		string h5resultsfn;

	private:
		//hold handler for the H5 file here
		herr_t status;
		hid_t fileid;
		hid_t groupid;
		hid_t mspcid;
		hid_t dsetid;
		hid_t dspcid;
		hid_t cparms;
		hid_t fspcid;

		hsize_t dims[2];
		hsize_t maxdims[2];
		hsize_t offset[2];
		hsize_t dimsnew[2];

		hsize_t nrows;
		hsize_t ncols;
		size_t BytesPerChunk;
		size_t RowsPerChunk;
		size_t ColsPerChunk;
		size_t RowsColsPerChunk;
		size_t NChunks;
		size_t TotalWriteHere; 			//specifies how many data elements in the currently processed dset were already written successfully
		size_t TotalWritePortion;		//specifies how many data elements are to be processed in an arbitrary portion as e.g. coming from threadlocal output
		size_t TotalWriteAllUsed;
		size_t TotalWriteAllChunk;			//specifies how many data elements ultimately have to become completed


		//temporary read/write buffer
		unsigned int* u32le_buf;
		double* f64le_buf;
	};

//#endif

#endif
