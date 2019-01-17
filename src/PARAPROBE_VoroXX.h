/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __PARAPROBE_VOROXX_H__
#define __PARAPROBE_VOROXX_H__


#include "PARAPROBE_HPDBScan.h"



#include "thirdparty/VoroRycroft/voro++-0.4.6/src/voro++.hh"
using namespace voro;


//user defined testing parameter
#define GUARDZONE_ION			2
#define VALIDZONE_ION			1
#define DO_NOT_DEREFERENCE		UINT32MX
//##MK::#define LDRVEX					25.0		//nm
//##MK::#define LDRVEY					25.0
//##MK::#define LDRVEZ					250.0


//MK::given the specific case of Voronoi tessellations to atoms and provided the following holds:
//+no Voronoi cell has more than 255 unique coordinate values of its vertices
//+no Voronoi cell needs more than 255 individual values to encode its face topology
//we can use the unsigned char datatype instead of int to store metadata and topological pieces of info
//thereby compressing the dataset inplace by a factor of 4x
#define IMPLICIT_STORING_UCHAR			1
#define IMPLICIT_STORING_UINT32			2
#define IMPLICIT_STORING_UINT64			3

#define TESS_BOUNDARYCORR_UNNECESSARY	+2
#define TESS_BOUNDARYCORR_SUCCESS		+1
#define TESS_BOUNDARYCORR_NOTRIANGLES	-1
#define TESS_BOUNDARYCORR_NOCUTPLANE	-2
#define TESS_BOUNDARYCORR_CELLDELETED	-3
#define TESS_BOUNDARYCORR_FAILED		-4

struct cellimplicitkeys
{
	unsigned char ntopologyvalues;			//how many individual topological pieces of information
	unsigned char ncoordinatevalues;		//how many individual coordinate values
	//##MK::this has 6B trailing padding but for now lets accept it, to HDF5 will write only array of unsigned char pairs!
};

ostream& operator<<(ostream& in, cellimplicitkeys const & val);

struct cell_write_meta
{
	//precision demand
	size_t prec;			//which precision is necessary to store topology info of the cell
	size_t n_cells;			//how many cells we have?
	size_t n_topo;			//how many elements with given precision for topology required
	size_t n_geom;			//how many elements with given precision for geometry required
	cell_write_meta() : prec(0), n_cells(0), n_topo(0), n_geom(0) {}
	cell_write_meta(const size_t _pr, const size_t _nc, const size_t _nt, const size_t _ng) :
		prec(_pr), n_cells(_nc), n_topo(_nt), n_geom(_ng) {}
};

ostream& operator<<(ostream& in, cell_write_meta const & val);


struct buildstats
{
	size_t CellsInside;
	size_t CellsHalo;
	size_t CellsWallcontact;
	size_t CellsWallTruncSuccess;
	size_t CellsWallTruncFail;
	size_t CellsInconstructible;
	size_t CellsIncorrectGeometry;
	size_t CellsEroded;

	size_t NumberOfFacetsProcessed;
	size_t NumberOfVertsProcessed;

	buildstats() : CellsInside(0), CellsHalo(0), CellsWallcontact(0),
			CellsWallTruncSuccess(0), CellsWallTruncFail(0),
			CellsInconstructible(0), CellsIncorrectGeometry(0), CellsEroded(0),
			NumberOfFacetsProcessed(0), NumberOfVertsProcessed(0) {}
	buildstats(const size_t _vcin, const size_t _vchalo, const size_t _vcwall,
			const size_t _vctrcy, const size_t _vctrcn, const size_t _vcerr,
			const size_t _vigeo, const size_t _vcerode, const size_t _nfprc, const size_t _nvprc ) :
					CellsInside(_vcin), CellsHalo(_vchalo), CellsWallcontact(_vcwall),
					CellsWallTruncSuccess(_vctrcy), CellsWallTruncFail(_vctrcn),
					CellsInconstructible(_vcerr), CellsIncorrectGeometry(_vigeo), CellsEroded(_vcerode),
					NumberOfFacetsProcessed(_nfprc), NumberOfVertsProcessed(_nvprc) {}
};

ostream& operator<<(ostream& in, buildstats const & val);



struct wallstats
{
	//captures logical information which threadlocal domain walls cut the cell
	bool xmi_touch;
	bool xmx_touch;
	bool ymi_touch;
	bool ymx_touch;
	bool zmi_touch;
	bool zmx_touch;
	bool any_touch;
	//only evaluated for cells of VALIDZONE_ION we check in addition whether the haloregion was sufficient
	//MK::specifically whether or not enough halo guard was provided such that the cell makes
	//no contact with particular domain walls to which there are halos if so this is a clear
	//sign that the cell was cut by the domain wall instead of facets of cells in the halo
	//for all VALIDZONE_IONs the halo_insufficient needs to be false then the geometry of these cells
	//is correct also across the domain boundaries...
	bool halo_insufficient;		//if true indicates that halo was insufficient,
	//retessellation with larger halo guard width is recommended
	wallstats() : 	xmi_touch(false), xmx_touch(false),
					ymi_touch(false), ymx_touch(false),
					zmi_touch(false), zmx_touch(false),
					any_touch(false), halo_insufficient(false) {}
	wallstats(const bool _xmi, const bool _xmx, const bool _ymi, const bool _ymx,
			const bool _zmi, const bool _zmx, const bool _any, const bool _haloproblem ) :
				xmi_touch(_xmi), xmx_touch(_xmx),
				ymi_touch(_ymi), ymx_touch(_ymx),
				zmi_touch(_zmi), zmx_touch(_zmx),
				any_touch(_any), halo_insufficient(_haloproblem) {}
};


struct wallcharacter
{
	//carries information whether the threadlocal domain for which the tessellation is
	//computed is logically connected to neighboring domains
	//*_ishalo
	//-->if true means there is a neighbor domain wall is a halo boundary
	//-->else false means there is no neighbor domain, wall marks end of global nonperiodic dataset
	bool xmi_ishalo;
	bool xmx_ishalo;
	bool ymi_ishalo;
	bool ymx_ishalo;
	bool zmi_ishalo;
	bool zmx_ishalo;
	bool pad1;
	bool pad2;
	wallcharacter() : 	xmi_ishalo(false), xmx_ishalo(false),
						ymi_ishalo(false), ymx_ishalo(false),
						zmi_ishalo(false), zmx_ishalo(false),
						pad1(false), pad2(false) {}
	wallcharacter(const bool _xmi, const bool _xmx,
			const bool _ymi, const bool _ymx,
			const bool _zmi, const bool _zmx ) :
				xmi_ishalo(_xmi), xmx_ishalo(_xmx),
				ymi_ishalo(_ymi), ymx_ishalo(_ymx),
				zmi_ishalo(_zmi), zmx_ishalo(_zmx),
				pad1(false), pad2(false) {}
};


struct cellstats
{
	p3d pos;
	double volume;
	int nfaces;
	cellstats() : pos(p3d()), volume(0.f), nfaces(0) {}
	cellstats(const p3d _p, const double _v, const int _nf) : pos(_p), volume(_v), nfaces(_nf) {}
};


class tessHdl;
class Tree;

class tess
{
	//implementing a thread-local object keeping the heavy- and metadata of a thread-local tessellation
public:
	tess();
	~tess();

	void set_zbounds( const p6d64 these );
	void set_zbounds( zlim const & these );
	void set_halosize_p3d64( p3d64 const & these );
	void set_halosize_p3d( p3d const & these );
	void set_haloinfo();
	void pull_ionportion_inplace( vector<p3dm3> const & in, aabb3d const & globalbox);
	void pull_ionportion( vector<p3dm3> const & in, aabb3d const & globalbox );

	void build_tessportion0();

	vector<unsigned int> localID2worldID; //maps a thread-local Voronoi cell ID to a global ion ID
	vector<cellstats> cellmeta;
	vector<wallstats> wallmeta;
	vector<string> pvecmeta;
	vector<aabb3d> cellaabb;
	buildstats mystats;

	//heavy topology and geometry data of scientific interest
	vector<voro_io_info> io_info;
	vector<unsigned int> io_ion2cell;
	vector<float> io_cellposition;
	vector<float> io_vol;
	vector<unsigned int> io_nfaces; //##MK::could be made unsigned char...
	vector<unsigned char> io_wall;
	vector<unsigned char> io_ion2type;

#ifndef VALIDZONE_IONS_ONLY
	vector<short> io_halo;
#endif
	vector<size_t> io_topo;
	vector<float> io_geom;			//implicit 1d contiguous triplets of values define one 3d point coordinate

	cell_write_meta myprecisiondemand;

	bool i_store_tess_metainfo;
	bool i_store_tess_cellpos;
	bool i_store_tess_topogeom;
	tessHdl* owner;
	Tree* surfacehull;
	vector<tri3d>* trianglehull;
	//vector<unsigned int>* ion2pp3idx;
	vector<unsigned int>* surfdistance_ids;
	vector<apt_xyz>* surfdistance_vals;
	vector<apt_xyz>* ion2pp3surf;

private:
	/*void collect_cell_heavydata( unsigned int jd, cell_write_meta const & info, vector<int> const & nbors,
			vector<int> const & fvert, vector<double> const & verts, const double vol );*/
	cell_write_meta collect_cell_heavydata2( const short cellstatus, const unsigned int jd, vector<int> const & nbors,
			vector<int> const & fvert, vector<double> const & verts );
	p3i identify_blockpartitioning( const size_t p_total, const apt_real p_perblock_target, aabb3d const & roi );
	string identify_cell_pvector( vector<int> const & nbors, vector<int> const & fvert );
	int truncate_cell_wallcontact( const p3d vcenter, voronoicell_neighbor & vcell, aabb3d & vcbox );
	int truncate_cell_wallcontact2( const p3d pion, voronoicell_neighbor & vcell, aabb3d & vcbox );
	wallstats identify_cell_wallcontact( vector<int> const & nbors );
	aabb3d identify_cell_aabb3d( vector<double> const & verts );
	cell_write_meta identify_cell_storage_demands( vector<int> const & nbors,
			vector<int> const & fvert, vector<double> const & verts);
	void identify_minset_of_coordtriplets( vector<int> const & nbors,
			vector<int> const & fvert, vector<double> const & verts,
			vector<int> & fvert_reindexed, vector<p3d> & verts_triplets );

	vector<p3dm3> ionportion;
	aabb3d mywindow;
	zlim zbounds;
	p3d halosize;
	wallcharacter haloinfo;
	map<string,int> mypvectors;
};


class solver;

class tessHdl
{
	//coordinating instance handling the individual thread local portions of a large tessellation
	//and coordinates pulling pieces of information from the local sub-tessellations
public:
	tessHdl();
	~tessHdl();

	void configure( const bool store_metainfo, const bool store_cellpositions, const bool store_topogeom );
	int threaded_tessellating( vector<p3dm3> const & in, aabb3d const & tipbox,
			vector<zlim> const & zbounds, p3d const & halosize );
	//int report_tessellation_chunked();
	//##MK::deprecated int report_tessellation_contiguous_atonce();
	int report_tessellation_hyperslab_perthread();
	int report_tessellation_hyperslab_onlythethread();

	h5Hdl myh5;
	xdmfHdl myxdmf;
	solver* owner;

//private:
	vector<tess*> local_tessellations;
	vector<voro_io_bounds> iohelp;
	cell_write_meta allprecisiondemand;
	bool we_store_tess_metainfo; //if no metainfo no topogeom plotting possible
	bool we_store_tess_cellpos;
	bool we_store_tess_topogeom;
};

#endif
