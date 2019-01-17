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


#ifndef __PARAPROBE_SOLVERHDL_H__
#define __PARAPROBE_SOLVERHDL_H__


#include "PARAPROBE_VoroXX.h"

//forward declaration classes
class solverHdl;
class reconstructor;
class threadmemory;
class decompositor;
class vxlizer;
class surfacer;
class rndlabeler;
class horderdist;
class clusterer;
class crystindexer;
class aptcrystHdl;
//##MK::class solver;


class solver
{
public:
	solver();
	~solver();

	solverHdl* owner;
	bool i_take_care( const size_t id, const unsigned long threadid, const unsigned int totalthreads );
	void volume_reconstruction();
	void spatial_decomposition();
	void volume_binning();
	void surface_triangulation();
	void surface_distancing2(); //##MK::faulty
	//##MK::deprecated VTK output void characterize_distances();
	void characterize_tip();
	void init_spatialindexing();
	void characterize_crystallography();
	void characterize_spatstat();
	void characterize_clustering();
	void tessellate_tipvolume();

	//computed location of ions in reconstructed space to accept during datamining spatially unorganized copy but reconstructed from solverHdl
	reconstructor* recon;
	decompositor* sp;
	vxlizer* binner;
	surfacer* surf;	//##MK::make independent from the binning
	rndlabeler* rndmizer;
	horderdist* hometrics;
	clusterer* cldetect;
	aptcrystHdl* aptcrystallo;
	tessHdl* tessellator;

	//MK::two structures to split the process of reconstruction from the rawdata_pos/_epos, sure we could pass the rawdata but then we
	//lack flexibility to instantiate multiple reconstructions by different algorithms hosted in the same solver

	class PeriodicTable mypse;				//periodic table instance to manage iontype identification

private:
	mathHdl mymath;							//internal mathematical subroutines
	bool healthy;							//is the reconstruction healthy or not
};


class threadmemory
{
public:
	threadmemory();
	~threadmemory();

	decompositor* owner;
	bool init_localmemory( p6d64 const & mybounds, p3d64 const & myhalo, const bool mlast );
	bool read_localions();
	inline size_t rectangular_transfer( const apt_xyz pos, const apt_xyz mi, const apt_xyz w, const apt_xyz mx, const size_t nvxl );
	void binning( sqb const & vxlgrid );
	void ion2surfdistance_bvh();
	void ion2surfdistance_init( const apt_xyz val );
	bool build_local_kdtree();

	inline apt_real get_zmi() const;
	inline apt_real get_zmx() const;
	inline bool me_last() const;

	vector<p3dm1> ionpp3;
	vector<size_t> ion2bin;
	vector<apt_xyz> ion2surf;

	//##MK::optimize
	vector<p3dm3> ionpp3_tess;
	vector<unsigned int> ionpp3_ionpp3idx;
	/*vector<unsigned int> distinfo_tess;
	//MK::SAME LENGTH THEN ionpp3_tess one-to-one correspondence
	//MK::stores an array index for ionpp3 which allows to identify which ion the respective entry in ionpp3_tess represents!
	//MK::avoids to store for every ion always an ID and thereby always to carry for every analysis task
	//larger array cache length than necessary, remind: only when the distancing info is to be used again
	//outside of the threadmemory scope, for instance in the Voro tess to identify how close an ion is
	//to the tip surface and use this information to eliminate cells with bias incorrect topology and geometry
	//MK::IN LATTER CASE ALLOWS ONLY to QUERY DISTANCE FOR VALIDZONE_IONS as their order of
	//populating ionpp3_tess during read_localions is instructed to proceed order preserving!*/

	//##MK::optimize, idea for instance, set up KDTree immediately would replace ionpp3 by ionpp3_kdtree and ion2surf by ion2surf_kdtree
	vector<p3dm1> ionpp3_kdtree;
	vector<apt_xyz> ion2surf_kdtree;

	kd_tree* threadtree;

	inline size_t get_memory_consumption();

private:
	p6d64 tessmii;
	p3d64 tesshalowidth;
	apt_real zmi;
	apt_real zmx;
	bool melast;
	mathHdl threadmath;
};


class decompositor
{
public:
	decompositor();
	~decompositor();

	solver* owner;
	void tip_aabb_get_extrema();
	void loadpartitioning();
	void reportpartitioning();

	vector<threadmemory*> db;				//threadlocal memory objects better suited for ccNUMA architectures
	vector<p6d64> spatialsplits;			//where is the point cloud split
	sqb rve;								//spatial bin geometry
	aabb3d tip;								//rigorous tip AABB bounds
	p3d64 halothickness;					//how much halo should be provided if tessellation desired
	bool kdtree_success;

private:
	bool healthy;
	mathHdl mymath;
};



class reconstructor
{
public:
	reconstructor();
	~reconstructor();

	solver* owner;

	bool reconstruction_accept_sequential( void );
	bool reconstruction_default_sequential( void );

	vector<vector<p3d>*> pp3;
	vector<vector<unsigned int>*> lbls; //##MK::redundant remove in the future!

private:
	bool healthy;							//is the reconstruction healthy or not
};


class vxlizer
{
public:
	vxlizer();
	~vxlizer();

	void characterize_binning();
	void identify_inside_bins();
	void identify_surface_adjacent_bins();
	bool identify_vacuum_bins();
	void binarization();
	bool allocate_binaries();
	sqb define_vxlgrid( const apt_xyz bwidth );
	void rectangular_binning();

	solver* owner;

	sqb vxlgrid;
	occupancy metadata;

	bool* IsVacuum;		//a bin with no ions laying in vacuum aka outside the tip
	bool* IsSurface;	//a bin with at least one ion in the Moore neighborhood of the vacuum
	bool* IsInside;		//a bin entirely enclosed with no Moore neighbor that is IsSurface == true
};


class surfacer
{
public:
	surfacer();
	~surfacer();

	solver* owner;

	bool read_vtk_io_tipsurface( const string vtk_io_fn );
	sqb define_binning_grid(const apt_xyz edgelen);
	bool new_pruning_mem();
	void del_pruning_mem();
	void seq_binarization();
	bool seq_hkbased_pruning();
	void seq_find_surface_adjacent_bins();
	void seq_filter_candidates( vector<p3d> & c );
	bool pruning_candidates( const apt_xyz db, vector<p3d> & cand );

	bool alphashape_core( vector<p3d> const & inp );
	bool alphashape();

	bool convexhull();
	bool marchingcubes();

	void report_tipsurface();

	bool build_rtree();
	void chop_rtree();

	vector<tri3d> tipsurface;
	Tree* bvh;

private:
	bool healthy;
	vector<p3d> candidates;

	bool* bitmap1; //helper to implement the pruning
	bool* bitmap2;
	aabb3d rve;
	sqb grid;

	string fn_alphashape;
	string fn_convexhull;
	string fn_marchcube;
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class rndlabeler
{
	//implements functionality to temporarily shuffle randomly the ion types to perform randomize descriptive spatial statistics
	//##MK::currently shuffling is done entirely deterministic by building global type array and shuffling it instead
	//building thread-local type arrays with thread-local PRNG to shuffle them
public:
	rndlabeler();
	~rndlabeler();
	void reset();

	void learn_original_types_global();
	void learn_original_types_cluster();
	void shuffle_types_mt19937();
	void apply_shuffled_types_global();
	void apply_shuffled_types_cluster();
	void reset_original_types_global();
	void reset_original_types_cluster();

	inline bool is_shuffled() const;
	inline bool is_applied() const;

	solver* owner;

	vector<unsigned int> initial_iontypid;
	vector<unsigned int> temporary_iontypid;

private:
	mt19937 urng;
	bool shuffled;
	bool applied;
};


class horderdist
{
public:
	horderdist();
	~horderdist();

	void initialize_descrstat_tasks();
	void compute_generic_spatstat2();   //for this reason here the thread team works successively through the region, with each thread team member getting a chunk of the region, downside here is thread fish frequently in memory of neighboring thread, but maximize load partitioning
	void compute_generic_spatstat3(); 	//removed at bounds check and repacks neighbor1dm1 results to avoid if ( type != desired ) continue branch mispredictions

	void report_apriori_descrstat2( const long tsktype, const string whichmetric, const string whichtarget, const string againstwhich, histogram & hist );
	void report_apriori_descrstat3( const string whichtarget, const string againstwhich, discrhistogram & hist );
	//##MK::single precision cumsum issues do not use! void report_apriori_descrstat4( const string whichtarget, const string againstwhich, npc3d & pcf );
	void report_npc3d_bincenters_hdf5( npc3d & pcf );
	void report_npc3d_histvalues_hdf5( const string whichtarget, const string againstwhich, npc3d & pcf );
	solver* owner;

private:
	bool healthy;

	vector<TypeSpecDescrStat> spatstat_tasks;
};


class clusterer
{
public:
	clusterer();
	~clusterer();

	void initialize_clustering_tasks( void );
	void maximum_separation_method( void );
	void maximum_separation_report( const string whichtarget,
			const string againstwhich, vector<dbscanres> const & results );

	solver* owner;

private:
	bool healthy;
	vector<TypeSpecDescrStat> clustering_tasks;
};


class clustertask
{
	//instance performing a single cluster analysis and subsequent spatial statistics with reduced ion cloud
public:

	clustertask();
	~clustertask();

	void initialize();
	dbscanres hpdbscan( const apt_real d, const size_t Nmin, const unsigned int runid );

	bool build_kdtree();
	void flag2exclude_guard( const apt_real dmx );
	void generic_spatstat( const unsigned int runid );
	void report_aposteriori_descrstat( const string whichmetric,
			const string whichtarget, const string againstwhich,
			const unsigned int rid, histogram & hist );
	void chop_kdtree();
	void unflag_all();
	void reset();

	clusterer* boss;
	TypeSpecDescrStat mission;
	size_t tskid;

	vector<p3dm1> ions_filtered;		//local copy of guys to work on, copy requires memory but allows to process more task characteristics specifically
	vector<apt_xyz> ion2surf_filtered;
	vector<p3dm1> ions_kdtree;
	vector<apt_xyz> ion2surf_kdtree;
	kd_tree* globaltree;
};


struct vicsmeta
{
	apt_real R;
	apt_real dR;
	apt_real binner;
	int NumberOfBins;
	vicsmeta() : R(0.f), dR(0.f), binner(0.f), NumberOfBins(0) {}
	vicsmeta( const apt_real _R, const apt_real _dR, const apt_real _binner, const int _nbins) :
		R(_R), dR(_dR), binner(_binner), NumberOfBins(_nbins) {}
};


inline bool SortAscElevAzimuthPairs( const pair<apt_real,apt_real> & a, const pair<apt_real,apt_real> & b)
{
	return a.second < b.second;
}


struct vicsresult
{
	apt_real elevation;
	apt_real azimuth;
	pair<float,float> max1; //bin, val
	pair<float,float> max2;
	pair<float,float> max3;
	/*bool fft_allc_success;
	bool fft_plac_success;
	bool fft_init_success;
	bool fft_comp_success;
	bool fft_free_success;
	bool pad1;
	bool pad2;
	bool pad3;*/
	vicsresult() : elevation(4*PI), azimuth(4*PI),
			max1(make_pair(0.f,0.f)),
			max2(make_pair(0.f,0.f)),
			max3(make_pair(0.f,0.f)) {}
	/*		fft_allc_success(false),
			fft_plac_success(false),
			fft_init_success(false),
			fft_comp_success(false),
			fft_free_success(false),
			pad1(false), pad2(false), pad3(false) {}*/
	vicsresult( const apt_real _el, const apt_real _az,
			pair<float,float> const & _m1,
			pair<float,float> const & _m2,
			pair<float,float> const & _m3 ) :
				elevation(_el), azimuth(_az),
				max1(_m1), max2(_m2), max3(_m3) {}
				/*fft_allc_success(false), fft_plac_success(false), fft_init_success(false),
				fft_comp_success(false), fft_free_success(false),
				pad1(false), pad2(false), pad3(false) {}*/
};


struct vicshistbin
{
	apt_real cnts;
	apt_real FFTr;
	apt_real FFTi;
	apt_real pw;
	vicshistbin() : cnts(0.f), FFTr(0.f), FFTi(0.f), pw(0.f) {}
	vicshistbin(const apt_real _c, const apt_real _fr, const apt_real _fi, const apt_real _pw) :
		cnts(_c), FFTr(_fr), FFTi(_fi), pw(_pw) {}
};




struct vicsfftsummary
{
	unsigned int threadID;
	unsigned int nFFTsSuccess;
	unsigned int nFFTsFailed;
	unsigned int nIons;
	vicsfftsummary() : threadID(0), nFFTsSuccess(0), nFFTsFailed(0), nIons(0) {}
	vicsfftsummary( const unsigned int _thr, const unsigned int _ffty,
			const unsigned int _fftn, const unsigned int _ni ) :
		threadID(_thr), nFFTsSuccess(_ffty), nFFTsFailed(_fftn), nIons(_ni) {}
};


class vics_materialpoint_result
{
public:
	vics_materialpoint_result();
	vics_materialpoint_result( const p3d here );
	~vics_materialpoint_result();

	p3d MatPointPos;
	vicsfftsummary FFTSummary;
	vector<vicsresult> ElevAzimTblSummary;
	//vector<vector<vicshistbin>> ElevAzimTblHistogram;
};


class crystindexer
{
public:
	crystindexer();
	~crystindexer();

	void configure();

	/*vicsresult ExecutePeaksFinding1D( vector<pair<double,double>> & CntsVsFreq );
	vicsfftstatus ExecuteSingleDiscreteFFT( vector<unsigned int> const & in, vector<pair<double,double>> & out );
	vicsresult SingleElevationAzimuth( const apt_real ee, const apt_real aa,
			vector<d3d> const & ions, vicsfftsummary & diary );
	void vicsmethod_core_incrementalfft( p3d const & target, vector<d3d> const & cand );*/
	void vicsmethod_core_incrementalfft2( p3d const & target, vector<d3d> const & cand );
	//void vicsmethod_core_batchfft_simd( p3d const & target, vector<d3d> const & cand );

	solver* owner;

	vicsmeta info;
	vector<float> window_coeff;
	vector<vics_materialpoint_result> res;
};


class aptcrystHdl
{
public:
	aptcrystHdl();
	~aptcrystHdl();

	void compute_crystallography_definegrid();
	void compute_crystallography_vicsmethod();
	void report_crystallography_results();
	void report_crystallography_results2();
	void delete_crystallography_results();

	solver* owner;
	cuboidgrid3d probehere;
	vector<p3d> samplingpoints;
	vector<p2d> elevazimpairs;
	vector<crystindexer*> workers;
};





class solverHdl
{
	//top-level construct implementing the worker instance at process rank within the MPI process level parallelism
	//there is only one solverHdl instance/object per process which may have different solver objects, one for every dataprocessing
	//reconstruction or mining) tasks. The solver again has internally a reconstructor and analyzer objects the latter two obey the solver object
	//the reconstructor reads the measurement raw data from its associated solverHdl, his owner and either feeds these through
	//into the x,y,z recon space values or performs a recon-algorithm to get these x,y,z recon values
	//the analyzer takes these x,y,z coordinates into a specific ##MK memory-locality aware data structure to perform analysis tasks, such as
	//tip surface reconstruction via CGAL and higher order metrics stuff

public:
	solverHdl();
	~solverHdl();

	bool generate_synthetic_tip();
	bool load_pos_sequentially( const string posfn );
	bool load_epos_sequentially( const string eposfn );
#ifdef UTILIZE_HDF5
	bool load_hdf5_sequentially( const string h5fn );
#endif
	bool compare_epos_pos( void );
	bool identify_ions( void );
	bool generate_hdf5_resultsfile();
	void delete_rawdata( void );

	/*
	/unsigned int plan_parameter_study_counts( void );
	unsigned int jobid2rank( const unsigned int jid );
	void plan_parameter_study( void );
	*/

	//setter, getter
	void set_mpidatatypes( void );
	inline int get_rank( void ) { return myRank; }
	inline int get_nranks( void ) { return nRanks; }
	void set_rank( const int rr ) { myRank = rr; }
	void set_nranks( const int nnrr ) { nRanks = nnrr; }

	//NUMA binding
	void initialize_thread_binding();

	spacebucket rawdata_alf;	//an artificial life form
	vector<vector<pos>*> rawdata_pos;
	vector<vector<epos>*> rawdata_epos;
	vector<vector<unsigned int>*> rawdata_iontype;
	unsigned int nevt;

	/*
	vector<jobreceipt> summary;
	vector<runparm> configuration;
	*/

	class PeriodicTable mypse;									//periodic table instance to manage iontype identification

	xdmfHdl xdmfmetaHdl;
	h5Hdl resultsh5Hdl;
	h5Hdl clusth5Hdl;
	h5Hdl crysth5Hdl;
	//h5Hdl voronoih5Hdl;
	profiler tictoc;

private:
	//MPI, OpenMP related
	int myRank;														//my MPI ID in the MPI_COMM_WORLD
	int nRanks;														//total number of MPI processes that work in the world
};

#endif
