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

#include "PARAPROBE_HPDBScan.h"
//#include "PARAPROBE_Profiling.h"

class solverHdl;
class reconstructor;
class threadmemory;
class decompositor;
class vxlizer;
class surfacer;
class rndlabeler;
class horderdist;
class clusterer;


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
	void characterize_distances();
	void init_spatialindexing();
	void characterize_spatstat();
	void characterize_clustering();

	//computed location of ions in reconstructed space to accept during datamining spatially unorganized copy but reconstructed from solverHdl
	reconstructor* recon;
	decompositor* sp;
	vxlizer* binner;
	surfacer* surf;	//##MK::make independent from the binning
	rndlabeler* rndmizer;
	horderdist* hometrics;
	clusterer* cldetect;

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
	bool init_localmemory( const apt_xyz zmin, const apt_xyz zmax, const bool mlast );
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

	//##MK::optimize, idea for instance, set up KDTree immediately would replace ionpp3 by ionpp3_kdtree and ion2surf by ion2surf_kdtree
	vector<p3dm1> ionpp3_kdtree;
	vector<apt_xyz> ion2surf_kdtree;

	kd_tree* threadtree;

	inline size_t get_memory_consumption();

private:
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

	sqb rve;								//spatial bin geometry
	aabb3d tip;								//rigorous tip AABB bounds
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

	//bool reconstruction_accept_synthetic( void );
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
	//void compute_generic_spatstat1();	//##MK::deprecated, threads work individually on only their own data, while this minimizes fishing in neighboring thread's memory, the downside is that for high thread team member counts individuals are passed too thin regions (close to bottom or top of tip) and Settings::SpatStatRadiusMax is large they do nothing and idle
	void compute_generic_spatstat2();   //for this reason here the thread team works successively through the region, with each thread team member getting a chunk of the region, downside here is thread fish frequently in memory of neighboring thread, but maximize load partitioning
	void report_apriori_descrstat1( const string whichmetric, const string whichtarget, const string againstwhich, histogram & hist );
	void report_apriori_descrstat2( const long tsktype, const string whichmetric, const string whichtarget, const string againstwhich, histogram & hist );

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
	void maximum_separation_report( const string whichtarget, const string againstwhich,
			vector<dbscanres> const & results );

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
	void flag2exclude_guard();
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

	profiler tictoc;

private:
	//MPI, OpenMP related
	int myRank;														//my MPI ID in the MPI_COMM_WORLD
	int nRanks;														//total number of MPI processes that work in the world
};

#endif
