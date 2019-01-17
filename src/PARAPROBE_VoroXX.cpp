/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#include "PARAPROBE_VoroXX.h"


ostream& operator<<(ostream& in, cellimplicitkeys const & val)
{
	in << "NCoordinatesVals/NTopologyVals\t\t" << val.ncoordinatevalues << ";" << val.ntopologyvalues << endl;
	return in;
}


ostream& operator<<(ostream& in, cell_write_meta const & val)
{
	in << "Precision/NCellsElementsTopoGeom\t\t" << val.prec << ";" << val.n_cells << ";" << val.n_topo << ";" << val.n_geom << endl;
	return in;
}



ostream& operator<<(ostream& in, buildstats const & val)
{
	in << "CellsInside\t\t" << val.CellsInside << "\n";
	in << "CellsInHalo\t\t" << val.CellsHalo << "\n";
	in << "CellWallcontact\t\t" << val.CellsWallcontact << "\n";
	in << "CellsInconstructible\t" << val.CellsInconstructible << "\n";
	in << "CellsIncorrectGeom\t" << val.CellsIncorrectGeometry << "\n";
	in << "CellsErodedAtSurface\t" << val.CellsEroded << "\n";
	in << "Truncation success/fail\t" << val.CellsWallTruncSuccess << ";" << val.CellsWallTruncFail << endl;
	return in;
}


tess::tess()
{
	i_store_tess_metainfo = false;
	i_store_tess_cellpos = false;
	i_store_tess_topogeom = false;

	owner = NULL;
	surfacehull = NULL;
	trianglehull = NULL;
	surfdistance_ids = NULL;
	surfdistance_vals = NULL;
	//ion2pp3surf = NULL;

	mywindow = aabb3d();
	zbounds = zlim();
	halosize = p3d();
	haloinfo = wallcharacter();

	mystats = buildstats();
	myprecisiondemand = cell_write_meta();
}


tess::~tess()
{
	//do not delete owner, surfacehull, trianglehull, ion2pp3idx, ion2pp3surf these are only a backreferences
}


void tess::set_zbounds( const p6d64 these )
{
	zbounds = zlim( these.zmi, these.zmx );
}


void tess::set_zbounds( zlim const & these )
{
	zbounds = these;
}


void tess::set_halosize_p3d64( p3d64 const & these )
{
	halosize = p3d( these.x, these.y, these.z );
}


void tess::set_halosize_p3d( p3d const & these )
{
	halosize = these;
}


void tess::set_haloinfo()
{
	//we use a z domain partitioning right now and partition a global nonperiodic tessellations
	//for this reason all x and y domain walls are ishalo = false
	//for omp_get_thread_num() MASTER upper z ishalo = true bottom false
	//for omp_get_threads_num()-1 lower z ishalo = true top false
	//for all other z regions top and bottom ishalo = true
	int mt = omp_get_thread_num(); //[0,omp_get_max_threads()-1]
	int nt = omp_get_num_threads();

	wallcharacter mysituation = wallcharacter();
	if ( nt > SINGLETHREADED ) {
		if ( mt > MASTER && mt < nt-1 ) {
			mysituation.zmi_ishalo = true;
			mysituation.zmx_ishalo = true;
		}
		else {
			if ( mt == MASTER ) {
				mysituation.zmi_ishalo = false;
				mysituation.zmx_ishalo = true;
			}
			else { //mt == nt-1
				mysituation.zmi_ishalo = true;
				mysituation.zmx_ishalo = false;
			}
		}
	}
	//nothing to reset in case of SINGLETHREADED as wallcharacter initialized by default to false

	haloinfo = mysituation;
}


void tess::pull_ionportion_inplace( vector<p3dm3> const & in, aabb3d const & globalbox)
{
	ionportion.clear();
	apt_real lzmii = zbounds.zmin;
	apt_real lzmxx = zbounds.zmax;
	apt_real zmii = zbounds.zmin - halosize.z;
	apt_real zmxx = zbounds.zmax + halosize.z;
	aabb3d bounds = aabb3d();
	int stats[3] = {0, 0, 0}; //out of range (should not occur), inside, in halo
	for( auto it = in.begin(); it != in.end(); ++it ) { // [zmii,zmx)
		//early reject
		if ( it->z < zmii || it->z >= zmxx ) {
			stats[0] +=1;
			continue;
		}
		//VALIDATOM or GUARDATOM
		if ( it->z >= lzmii && it->z < lzmxx ) {
			ionportion.push_back( p3dm3(it->x, it->y, it->z, it->id, VALIDZONE_ION, it->iontype) );
			stats[1] +=1 ;
		}
		else { //in guardzone
			ionportion.push_back( p3dm3(it->x, it->y, it->z, it->id, GUARDZONE_ION, it->iontype) );
			stats[2] += 1;
		}
		//MK::see that using such looping the implicit order of ion names from in is transferred to ionportion

		//##MK::old version make local box independent of global context
		/*if ( it->x <= bounds.xmi )		bounds.xmi = it->x;
		if ( it->x >= bounds.xmx )		bounds.xmx = it->x;
		if ( it->y <= bounds.ymi )		bounds.ymi = it->y;
		if ( it->y >= bounds.ymx )		bounds.ymx = it->y;*/

		if ( it->z <= bounds.zmi )		bounds.zmi = it->z;
		if ( it->z >= bounds.zmx )		bounds.zmx = it->z;
	}

	#pragma omp critical
	{
		cout << "Thread " << omp_get_thread_num() << " outofrange/valid/guard/total\t\t" << stats[0] << ";" << stats[1] << ";" << stats[2] << "\t\t" << (stats[0]+stats[1]+stats[2]) << endl;
 	}

	//make local box aware of global context the domain walls where we do not partition the domain
	bounds.xmi = globalbox.xmi;
	bounds.xmx = globalbox.xmx;
	bounds.ymi = globalbox.ymi;
	bounds.ymx = globalbox.ymx;

	mywindow = bounds;
	mywindow.add_epsilon_guard();
	mywindow.scale();
}



void tess::pull_ionportion( vector<p3dm3> const & in, aabb3d const & globalbox )
{
	ionportion.clear();
	apt_real zmii = zbounds.zmin - halosize.z;
	apt_real zmxx = zbounds.zmax + halosize.z;
	aabb3d bounds = aabb3d();
	int stats[2] = {0, 0};
 	for( auto it = in.begin(); it != in.end(); ++it ) { // [zmii,zmx)
		//early reject
		if ( it->z < zmii )
			continue;
		if ( it->z >= zmxx )
			continue;

		//VALIDATOM or GUARDATOM
		if ( it->z >= zbounds.zmin && it->z < zbounds.zmax ) {
			ionportion.push_back( p3dm3(it->x, it->y, it->z, it->id, VALIDZONE_ION, it->iontype) );
			stats[0] +=1 ;
		}
		else { //in guardzone
			ionportion.push_back( p3dm3(it->x, it->y, it->z, it->id, GUARDZONE_ION, it->iontype) );
			stats[1] += 1;
		}

		//##MK::old version make local box independent of global context
		/*if ( it->x <= bounds.xmi )		bounds.xmi = it->x;
		if ( it->x >= bounds.xmx )		bounds.xmx = it->x;
		if ( it->y <= bounds.ymi )		bounds.ymi = it->y;
		if ( it->y >= bounds.ymx )		bounds.ymx = it->y;*/

		if ( it->z <= bounds.zmi )		bounds.zmi = it->z;
		if ( it->z >= bounds.zmx )		bounds.zmx = it->z;
	}

	#pragma omp critical
	{
		cout << "Thread " << omp_get_thread_num() << " valid/guard/total\t\t" << stats[0] << ";" << stats[1] << "\t\t" << (stats[0]+stats[1]) << endl;
 	}

	//make local box aware of global context the domain walls where we do not partition the domain
	bounds.xmi = globalbox.xmi;
	bounds.xmx = globalbox.xmx;
	bounds.ymi = globalbox.ymi;
	bounds.ymx = globalbox.ymx;

	mywindow = bounds;
	mywindow.add_epsilon_guard();
	mywindow.scale();
}

/*
void tess::collect_cell_heavydata( unsigned int jd, cell_write_meta const & info, vector<int> const & nbors,
		vector<int> const & fvert, vector<double> const & verts, const double vol )
{
	//##MK::bounds check, narrowing contraction from size_t to unsigned int!
	io_info.push_back( voro_io_info(jd, info.n_topo, info.n_geom) );

	int j = 0;
	for( size_t i = 0; i < nbors.size(); ++i ) { //O(N)lg(N) identification of how many faces of certain type
		int nvertices_thisface = fvert[j];
		io_topo.push_back( nvertices_thisface );
		for( int k = 0; k < fvert[j]; k++) {
			//generating 3d points with implicit indices 0, 1, 2, .... on the double array of coordinate values belonging to the cell
			int l = 3*fvert[j+k+1];
			io_topo.push_back( l );
			io_topo.push_back( l+1 );
			io_topo.push_back( l+2 );
			//#####MK::accounting every vertex has three unsigned char indices!
		}
		j += fvert[j] + 1;
	}

	for( auto it = verts.begin(); it != verts.end(); ++it ) {
		io_geom.push_back( it->x );
		io_geom.push_back( it->y );
		io_geom.push_back( it->z );
	}

	io_vol.push_back( vol );
}
*/

cell_write_meta tess::collect_cell_heavydata2( const short cellstatus, const unsigned int jd,
		vector<int> const & nbors, vector<int> const & fvert, vector<double> const & verts )
{
	//when reporting the cell facet vertices in the global coordinate system the Voro++ library reports
	//a list of all the unique 1d coordinate values not the coordinate value triplets!
	//for visualization of the cell geometry with at least to my knowledge all tools
	//coordinate value triplets are required. If reported and stored naively this bloats up the dataset
	//unnecessarily, for instance take a cell with 36 unique float coordinate values and 108 vertices
	//when reporting the polygonal facet that build the cell, these 3d points are build from the
	//36 values, in effect multiple points will have the same float triplet
	//ideally one would store only the 36 floats and keep a datastructure of integer references which
	//instructs the renderer with which the vis software is instructed to generate the float triplets on the fly
	//neither with VTK nor XDMF/HDF5 this is necessarily efficient though
	//a such we seek to identify the unique float triplets and re-index the mapping of vertex to float triplets
	//this though requires to identify all unique combinations of three integer reference indices as
	//stored in fvert
	//finding unique float triplets is not necessarily numerically robust as floats have system-specific
	//realization and NaNs have to be handled properly
	//finding unique int, short or char triplets is robust however if implemented naive O(N^2) complex
	//instead we use here an approach based on the map container which finds existent entries in lg(N) time
	//resulting in O(NlgN) complexity, using an unordered_set has worst case complexity N, best constant time
	//so potentially ##MK::switching to unordered_set or unordered_map might give improvement

	//given that we seek to write the entire set of cells into an HDF5 file we cannot a priori
	//assume that the total number of vertex triplets stays below numeric_limits<unsigned int>::max()
	//at the same time we do not want to store the threadlocal tessellations individually but
	//write as if it is only one large tessellation, as such we strictly would have to define
	//the number of references to vertex coordinate triplets as size_t
	//for smaller tessellations typically encountered in most APT practical examples this is overkill
	//as only the lower likely 32bit of the 64bit integer will be occupied
	//##MK::hence we write the code as a first shot to use unsigned int+complain potentially via bounds
	//if really billion sized tessellations are to be visualized one would have to split
	//the population into cells for instance either based on their pvectors
	//or based on the number of cells
	//##MK::for both cases in particular the first implementations into HDF5 is tedious
	//as such we stick for now with this simpler version to have a visualization for the
	//first time within the community at all, if later there is a necessity for implementing either of the two
	//(I advise to use the pvector-based) one can add it
	cell_write_meta out = cell_write_meta();

	if ( i_store_tess_topogeom == false ) {
		io_info.push_back( voro_io_info(jd, nbors.size(), 0, 0) ); //##MK::in this case last two unused reporting n_geom how many TRIPLETS NOT SINGLE VALUES ! ##MK::bounds checks
		return out; //##MK::in this case no change on out necessary
	}
	else { //i_store_tess_topogeom == true )
		//first step: build list of unique p3d verts_triplets
		map<int,int> old2new;
		vector<int> fvert_reindexed;
		vector<p3d> vtriplets_unique;
		//because Voro++ reports vertices as a set of the unique 1d coordinate values
		//keys of the maps are build from the old verts indices where int indices u,v,w are converted into implicit 1d key
		//values of the map are the new indices referring to p3d verts_triplets specific to only the current cell
		//for subsequent I/O these packets of p3d triplets are fused together in a large array using an implicit position offset on the 1d array
		int ii = 0;
		//count number of unique float triplets
		int j = 0;
		for( size_t i = 0; i < nbors.size(); ++i ) { //O(N) processing
			for( int k = 0; k < fvert[j]; k++) { //generating 3d points with implicit indices 0, 1, 2, ....
				int l = 3*fvert[j+k+1];
				int implicitkey = l + (l+1)*1024 + (l+2)*1024*1024; //##MK::
				auto it = old2new.find( implicitkey );
				if ( it != old2new.end() ) { //O(lg(N)) found integer triplet exists already
					fvert_reindexed.push_back( it->second );
				}
				else { //integer triplet does not yet exist so create
					old2new.insert( make_pair( implicitkey, ii ) );
					fvert_reindexed.push_back( ii );
					vtriplets_unique.push_back( p3d(verts[l], verts[l+1], verts[l+2]) );
					ii++;
				}
			}
			j += fvert[j] + 1;
		}

		//second step: pass heavy data to threadlocal buffer
		//##MK::could at some point be fused with first step to save loop overhead...


		ii = 0;
		j = 0;
		for( size_t i = 0; i < nbors.size(); ++i ) { //O(N)lg(N) identification of how many faces of certain type
			unsigned int nvertices_thisface = fvert[j];
			io_topo.push_back( 3 ); //XDMF keyword to visualize an n-polygon
			io_topo.push_back( nvertices_thisface ); //promotion uint23 to size_t no problem
			//for( int k = 0; k < fvert[j]; k++) { //it is essential to sweep through this loop nest in the same order than above for carrying over the prescribed order on fvert_reindexed
			for( unsigned int k = 0; k < nvertices_thisface; k++, ii++ ) {
				io_topo.push_back( fvert_reindexed.at(ii) ); //promotion of int to size_t not a problem if int >= 0
			}
			j += fvert[j] + 1;
#ifndef VALIDZONE_IONS_ONLY
			io_halo.push_back( cellstatus );	//for every facet we need an own halo info because in XDMF/HDF5 the objects are the facets not the polyhedron! i.e as many facets as neighbors per polyhedron facet
			//this is a workaround which patches the faulty/nonexistent XDMF Paraview polyhedron support
#endif
		}

		for( auto it = vtriplets_unique.begin(); it != vtriplets_unique.end(); ++it ) {
			io_geom.push_back( it->x );
			io_geom.push_back( it->y );
			io_geom.push_back( it->z );
		}

		io_info.push_back( voro_io_info(jd, nbors.size(),
				fvert_reindexed.size() + 2*nbors.size(),
				vtriplets_unique.size()) ); //reporting n_geom how many TRIPLETS NOT SINGLE VALUES ! ##MK::bounds checks

		out.n_topo = fvert_reindexed.size() + 2*nbors.size();
		out.n_geom = vtriplets_unique.size() * 3;
		return out;
	}
}


p3i tess::identify_blockpartitioning( const size_t p_total, const apt_real p_perblock_target, aabb3d const & roi )
{
	p3i out = p3i( 1, 1, 1);
	if ( roi.xsz > EPSILON ) {
		apt_real yrel = roi.ysz / roi.xsz;
		apt_real zrel = roi.zsz / roi.xsz;
		apt_real ix = pow(
				(static_cast<apt_real>(p_total) / p_perblock_target) / (yrel * zrel),
				(1.f/3.f) );

		out.x = ( static_cast<int>(floor(ix)) > 0 ) ? static_cast<int>(floor(ix)) : 1;
		apt_real iy = yrel * ix;
		out.y = ( static_cast<int>(floor(iy)) > 0 ) ? static_cast<int>(floor(iy)) : 1;
		apt_real iz = zrel * ix;
		out.z = ( static_cast<int>(floor(iz)) > 0 ) ? static_cast<int>(floor(iz)) : 1;
	}
	return out;
}


string tess::identify_cell_pvector( vector<int> const & nbors, vector<int> const & fvert )
{
	//pvector see E. A. Lazar, J. K. Mason, R. D. MacPherson, D. J. Srolovitz
	//Complete topology of cells, grains, and bubbles in three-dimensional microstructures
	//definition of a pvector: i^th position tells how many faces the cell has with i-vertices
	//example a cube has six square faces hence 0006
	//a problem with the display of the pvector in Lazar's work though is what if the Voronoi cell
	//has more than 9 faces of a specific kind, to cure we add a minus sign i.e. 0006 becomes 0-0-0-6
	map<int,int> pvector;
	int j = 0;
	for( size_t i = 0; i < nbors.size(); ++i ) { //O(N)lg(N) identification of how many faces of certain type
		int nvertices_thisface = fvert[j];
		auto it = pvector.find(nvertices_thisface);
		if ( it != pvector.end() )
			it->second++;
		else
			pvector.insert( pair<int,int>(nvertices_thisface, 1) );
		j += fvert[j] + 1;
	}

	/*for ( auto it = pvector.begin(); it != pvector.end(); ++it ) {
		cout << "Which/Count\t\t" << it->first << "\t\t" << it->second << endl;
	}*/

	stringstream out;
	if( pvector.empty() == false ) {
		int highest_type = (--pvector.end())->first; //pvector contains only elements with counts > 0 !
		//cout << "Highest type " << highest_type << endl;
		for ( int i = 1; i < highest_type; i++ ) { //add highest type manually to avoid adding trailing separator character!
			auto it = pvector.find(i);
			if ( it != pvector.end() ) //element found so get counts of such face
				out << it->second << "-";
			else //element not found so no face of type i
				out << "0-";
		}
		out << (--pvector.end())->second;
		return out.str();
	}
	else {
		return "";
	}
}


int tess::truncate_cell_wallcontact( const p3d vcenter, voronoicell_neighbor & vcell, aabb3d & vcbox )
{
	//the cell topology, geometry and data object and also the vcbox aabb3d will be modified
	//there is no physically exact boundary correction procedure strictly speaking
	//because cell wraps an atom at the dataset boundary
	//even an ideal experiment of though APT measurement though can only recover the
	//discrete nature of the material, hence it is a matter of definition what one considers as dataset boundary
	//several strategies exist:
	//define plane3d which truncates the cell at some distance from ion position
	//--> this is arbitrary physically because how to justify that any particular distance is better than another
	//using the Voro++ functionality to do c.plane(ix,iy,iz) allows to cut the Voronoi cell by an
	//imaginary particle placement at ix,iy,iz
	//pros: for greedy truncation, just take mean of nearby triangle hull barycenter O(n) time!
	//cons: implicit inclination control only a better approach than the greedy is to
	//define additional local coordinate system to consistent outer unit normal to ion position
	//therewith compute where the particle has to be placed such that bisector cutplane of ion cell and particle fake cell
	//has exactly position and inclination desired, then however plane normal required so SVD or Eigendecomposition O(n^2.36) <= O(n^3) complexity...
	//use exact alpha shape triangle hull and cut with vcell polyhedron
	//--> numerically difficult and involved, need to identify closed loop of triangle/polygon facet lines + resulting volume
	//--> alphashape not necessarily watertight
	//--> not any better than first method as alpha shape was constructed by assuming the boundary connects the atoms, i.e.
	//the dataset boundary defining atoms are laying on the boundary
	//definition we use the first strategy

	//find triangles of the alpha shape dataset boundary whose AABB intrude into the cells current AABB
	//prune most triangles by utilizing bounded volume hierarchy (BVH) Rtree,
	//build axis-aligned bounding box about current Voronoi cell
	AABB e_aabb( trpl(vcbox.xmi, vcbox.ymi, vcbox.zmi), trpl(vcbox.xmx, vcbox.ymx, vcbox.zmx) );
	vector<unsigned int> candidates;
	candidates = surfacehull->query( e_aabb );
	//mathHdl mymath;

	if ( candidates.size() < 3 )
		return TESS_BOUNDARYCORR_NOTRIANGLES;
	else { //>= 3
		//extract barycenter of these triangles
		//vector<p3d> plane_support( candidates.size(), p3d() );
		p3d cutpoint = p3d(0.f, 0.f, 0.f);

		for( size_t i = 0; i < candidates.size(); ++i ) {
			unsigned int triidx = candidates.at(i);
			tri3d cand = trianglehull->at( triidx );
			p3d tribary = cand.barycenter(); //greedy way of doing it
			cutpoint.x = cutpoint.x + tribary.x;
			cutpoint.y = cutpoint.y + tribary.y;
			cutpoint.z = cutpoint.z + tribary.z;
			//plane_support.at(i) = cand.barycenter();
		}
		apt_xyz averager = static_cast<apt_xyz>(candidates.size());
		cutpoint.x /= averager;
		cutpoint.y /= averager;
		cutpoint.z /= averager;

		bool voroxx_status = vcell.plane( cutpoint.x, cutpoint.y, cutpoint.z );
		if ( voroxx_status == true) {
			return TESS_BOUNDARYCORR_SUCCESS;
		}
		else { //false, if the plane cut deleted the cell entirely
			return TESS_BOUNDARYCORR_CELLDELETED;
		}
	}
		/*
		//use collection of these points to fit plane3d
		plane3d cutplane = plane3d();
		if ( mymath.fit_leastsqr_plane3d( plane_support, cutplane ) == true ) {
			//use the plane to truncate the cell
			//Voro++ documentation http://math.lbl.gov/voro++/doc/refman/classvoro_1_1voronoicell__neighbor.html#ab71c4487bbca2d20b00e9af2fc6d1ed0
			//bool voro::voronoicell_neighbor::plane(	double x,double y,double z,	double rsq ) inline
			//This version of the plane routine just makes up the plane ID to be zero. It will only be referenced if
			//neighbor tracking is enabled.

			p3d64 pcutplane = p3d64(); //###########need a point on the cutplane
			double rsqrd = 1.0;  //need squared distance...

			bool voroxx_status = vcell.plane( pcutplane.x, pcutplane.y, pcutplane.z, rsqrd );
			if ( voroxx_status == true) {
				return TESS_BOUNDARYCORR_SUCCESS;
			}
			else { //false, if the plane cut deleted the cell entirely
				return TESS_BOUNDARYCORR_CELLDELETED;
			}
		}
		else {
			return TESS_BOUNDARYCORR_NOCUTPLANE;
		}
		 */
}


int tess::truncate_cell_wallcontact2( const p3d pion, voronoicell_neighbor & vcell, aabb3d & vcbox )
{
	//see further comments to truncate_cell_wallcontact
	//implements the following strategy:
	//i) find all ions whose fattened by R AABB, does intersect with at least three triangles whose
	//barycenter locations are inside the spherical region R about the ion
	//ii) use next these barycenter locations to fit a plane
	//iii) find generating point of a fake voronoi cell ptarget whose bisecting plane to the ion is inclined
	//as the unit normal of the fit plane, MK::assume that the ion is always behind or on the fit plane

	//find triangles of the alpha shape dataset boundary whose AABB intrude into the cells current AABB
	//prune most triangles by utilizing bounded volume hierarchy (BVH) Rtree,
	//build axis-aligned bounding box about current Voronoi cell
	apt_real R = 2.f; //2nm
	apt_real RSQR = SQR(2.f);

	//pruning ions with more than R + 0
	AABB e_aabb( trpl(pion.x-R, pion.y-R, pion.z-R), trpl(pion.x+R, pion.y+R, pion.z+R) );
	vector<unsigned int> first_cand;
	first_cand = surfacehull->query( e_aabb );

	//##MK::future optimization use, if existent, precomputed ion to surface distances to prune ions
	//this costs only a reading of a memory position, potentially a local one if value is delivered upon call
	//instead of doing a costly tree traversal
	if ( first_cand.size() < 3 ) { //MK::either because cell is deep in tip volume or
		//by chance ashape holes are large and thereby avoid to find >=3 triangles with barycenter inside sphere about pion
		return TESS_BOUNDARYCORR_UNNECESSARY;
		//##MK::STRICTLY SPEAKING THIS IS ONLY CORRECT IF ALPHASHAPE IS A MANIFOLD/HAS NOT HOLES!
		//####MK::which we have not tested so far...
		//return TESS_BOUNDARYCORR_NOTRIANGLES;
	}
	else { //>= 3
		//not only we need some three triangles but which have barycenter in sphere about pion
		vector<p3d> plane_support;
		for( size_t i = 0; i < first_cand.size(); ++i ) {
			unsigned int triidx = first_cand.at(i);
			tri3d cand = trianglehull->at( triidx );
			p3d tribary = cand.barycenter(); //##MK::optimize for speed
			if ( (SQR(tribary.x-pion.x)+SQR(tribary.y-pion.y)+SQR(tribary.z-pion.z)) <= RSQR )
				plane_support.push_back( tribary );
		}

		if ( plane_support.size() < 3 ) {
			return TESS_BOUNDARYCORR_NOTRIANGLES;
		}
		else { //>=3
			mathHdl mymath;
			plane3d bisectpln = plane3d(); //the bisecting plane with which we want to cut the Voronoi cell
			p3d bisectpoi = p3d(); //a test point on the bisecting plane
			if ( mymath.fit_leastsqr_plane3d( plane_support, bisectpoi, bisectpln ) == true ) {
				//is pion behind or on the plane?
				apt_real d = 	bisectpln.n0 * (bisectpoi.x - pion.x) +
								bisectpln.n1 * (bisectpoi.y - pion.y) +
								bisectpln.n2 * (bisectpoi.z - pion.z);
				if ( d > 0.0 ) { //sgn of dotproduct indicates is in front however
					//normal may be inconsistent so ASSUME, because we know the physical configuration,
					//that pion is behind, so flip normal
					bisectpln.n0 *= -1.f;
					bisectpln.n1 *= -1.f;
					bisectpln.n2 *= -1.f;
				}
				//now normal points away from pion to bisectpoi

				//now the fake generating point to build a bisecting plane as bisectpln to the cell of pion is
				//is to walk 2 times the distance from pion to the projected point whose normal
				//shoots through pion on the bisectpln in the direction of the normal n012
				p3d fakepoi = p3d();
				apt_real lambda = 2.f * abs(bisectpln.n0*pion.x+bisectpln.n1*pion.y+bisectpln.n2*pion.z);
				fakepoi.x = pion.x + lambda * (bisectpoi.x - pion.x);
				fakepoi.y = pion.y + lambda * (bisectpoi.y - pion.y);
				fakepoi.z = pion.z + lambda * (bisectpoi.z - pion.z);

				//use the plane to truncate the cell
				//Voro++ documentation http://math.lbl.gov/voro++/doc/refman/classvoro_1_1voronoicell__neighbor.html#ab71c4487bbca2d20b00e9af2fc6d1ed0
				//bool voro::voronoicell_neighbor::plane(	double x,double y,double z,	double rsq ) inline
				//This version of the plane routine just makes up the plane ID to be zero. It will only be referenced if
				//neighbor tracking is enabled.

				if ( vcell.plane( fakepoi.x, fakepoi.y, fakepoi.z ) == true )
					return TESS_BOUNDARYCORR_SUCCESS;
				//##MK::this cannot work in cases of an undulating or wavy local surface because
				//here the detailed cruvature of the patch defines were the barycenter build the plane
				//resulting potentially in pion in front of this face causing now a truncation of the
				//cell from the tip, like cutting it off with a scissor....

				else
					return TESS_BOUNDARYCORR_CELLDELETED; //false, if the plane cut deleted the cell entirely
			}
			else
				return TESS_BOUNDARYCORR_NOCUTPLANE;
		}
	}
}


wallstats tess::identify_cell_wallcontact( vector<int> const & nbors )
{
	//does cell make boundary contact or not ? (indicated by negative neighbor
	//numeral info of neighbor details which domain wall cuts the cell into a potential incorrect geometry
	wallstats out = wallstats();
	for( auto it = nbors.begin(); it != nbors.end(); ++it ) {
		if ( *it >= 0 ) {
			continue;
		}
		else { //a cell that has boundary contact
			int thisone = abs(*it); //-1 to 1, -2 to 2, and so forth, Voro++ domain wall neighbors are -1,-2,-3,-4,-5,-6
			switch (thisone)
			{
				case 1:
					out.xmi_touch = true;	out.any_touch = true;
					break;
				case 2:
					out.xmx_touch = true;	out.any_touch = true;
					break;
				case 3:
					out.ymi_touch = true;	out.any_touch = true;
					break;
				case 4:
					out.ymx_touch = true;	out.any_touch = true;
					break;
				//##MK::we use a partitioning into threadlocal z domains of a global nonperiodic tess
				//##MK::hence check additionally whether haloregion was sufficient for
				//truncating the cell only through the halo regions and not the domain walls
				//##MK::this consistency check identifies whether cells close to the
				//threadlocal domain boundaries have correct geometry
				case 5:
					out.zmi_touch = true;	out.any_touch = true;
					if ( haloinfo.zmi_ishalo == true && out.halo_insufficient == false ) {
						out.halo_insufficient = true;
					}
					break;
				case 6:
					out.zmx_touch = true;	out.any_touch = true;
					if ( haloinfo.zmx_ishalo == true && out.halo_insufficient == false ) {
						out.halo_insufficient = true;
					}
					break;
				default:
					break;
			}
		}
	} //all neighbors have to be tested to know for sure with which boundaries we make contact
	return out;
}


aabb3d tess::identify_cell_aabb3d( vector<double> const & verts )
{
	//find bounding box about Voronoi cell
	aabb3d out = aabb3d();
	for ( size_t i = 0; i < verts.size() / 3; ++i) {
		double xyz[3] = { verts.at(3*i+0), verts.at(3*i+1), verts.at(3*i+2) };
		if ( xyz[0] <= out.xmi )	out.xmi = xyz[0];
		if ( xyz[0] >= out.xmx )	out.xmx = xyz[0];
		if ( xyz[1] <= out.ymi )	out.ymi = xyz[1];
		if ( xyz[1] >= out.ymx )	out.ymx = xyz[1];
		if ( xyz[2] <= out.zmi )	out.zmi = xyz[2];
		if ( xyz[2] >= out.zmx )	out.zmx = xyz[2];
	}
	return out;
}


cell_write_meta tess::identify_cell_storage_demands( vector<int> const & nbors, vector<int> const & fvert, vector<double> const & verts)
{
	//number of edges for largest n polygon cell facet and total number of vertex coordinates value only
	//determines which minimum precision we need to store topological references...
	int demand = 0;
	int sum = 0; //sum up how many topological references there are in total
	int j = 0;
	for( size_t i = 0; i < nbors.size(); ++i ) {
		int nvertices_thisone = fvert[j];
		if ( nvertices_thisone >= demand ) { //number of edges per largest facet
			demand = nvertices_thisone;
		}
		sum = sum + (3*nvertices_thisone + 1);
		//*3 because every vertex 3d coordinate needs three references to coordinate values holded in verts
		//+1 because we need to store how many vertices current facet is build of
		j += fvert[j] + 1;
	}
	//MK::sum of (number of vertices per facet + one value to identify npolygon vertices count) for all facets
	//the sum itself does not affect precision demand because it defines only the length of the array but not its precision demand to allow
	//to hold a sufficiently discerning number of disjoint indices to topology or vertex coordinate values

	//number of unique vertex coordinate values defines maximum integer reference
	cell_write_meta out = cell_write_meta();
	if ( verts.size() > demand )
		out.prec = verts.size();
	else
		out.prec = demand;
	out.n_topo = sum; //MK::aim is to use only unsigned char keys to dereference faces, i.e. autocompression of topology data
	//yet we can store much more than 255 of such unsigned char values,
	//hence precision refers here only to the minimum precision to store the topology keys
	//in other words are there more than 256 coordinate values or facets with more than 255 edges we need to
	//use an unsigned short or even unsigned int to store the faces
	//clearly for smaller Voronoi cells as typically discussed in studies of interfaces one could
	//do quick-and-dirty and store always int to be sure it fits, however this scales poorer and creates
	//also for small datasets unnecessary overhead, sure one could use compression as bloating up the datatype
	//contribute low entropy information that efficiently compresses but why doing so in post?
	out.n_geom = verts.size();
	return out;
}


void tess::identify_minset_of_coordtriplets( vector<int> const & nbors,
			vector<int> const & fvert, vector<double> const & verts,
			vector<int> & fvert_reindexed, vector<p3d> & verts_triplets )
{
	//when reporting the cell facet vertices in the global coordinate system the Voro++ library reports
	//a list of all the unique 1d coordinate values not the coordinate value triplets!
	//for visualization of the cell geometry with at least to my knowledge all tools
	//coordinate value triplets are required. If reported and stored naively this bloats up the dataset
	//unnecessarily, for instance take a cell with 36 unique float coordinate values and 108 vertices
	//when reporting the polygonal facet that build the cell, these 3d points are build from the
	//36 values, in effect multiple points will have the same float triplet
	//ideally one would store only the 36 floats and keep a datastructure of integer references which
	//instructs the renderer with which the vis software is instructed to generate the float triplets on the fly
	//neither with VTK nor XDMF/HDF5 this is necessarily efficient though
	//a such we seek to identify the unique float triplets and re-index the mapping of vertex to float triplets
	//this though requires to identify all unique combinations of three integer reference indices as
	//stored in fvert
	//finding unique float triplets is not necessarily numerically robust as floats have system-specific
	//realization and NaNs have to be handled properly
	//finding unique int, short or char triplets is robust however if implemented naive O(N^2) complex
	//instead we use here an approach based on the map container which finds existent entries in lg(N) time
	//resulting in O(NlgN) complexity, using an unordered_set has worst case complexity N, best constant time
	//so potentially ##MK::switching to unordered_set or unordered_map might give improvement

	map<int,int> old2new;
	//keys are build from old verts indices where u,v,w is converted into implicit 1d key
	//values are new indices referring to p3d verts_triplets
	int ii = 0;
	//count number of unique float triplets
	int j = 0;
	int debug = 0;
	for( size_t i = 0; i < nbors.size(); ++i ) { //O(N) processing
		for( int k = 0; k < fvert[j]; k++) { //generating 3d points with implicit indices 0, 1, 2, ....
			int l = 3*fvert[j+k+1];
			int implicitkey = l + (l+1)*1024 + (l+2)*1024*1024; //##MK::
			auto it = old2new.find( implicitkey );
			if ( it != old2new.end() ) { //O(lg(N)) found integer triplet exists already
				fvert_reindexed.push_back( it->second );
			}
			else { //integer triplet does not yet exist so create
				old2new.insert( make_pair( implicitkey, ii ) );
				fvert_reindexed.push_back( ii );
				verts_triplets.push_back( p3d(verts[l], verts[l+1], verts[l+2]) );
				ii++;
			}
			debug++;
		}
		j += fvert[j] + 1;
	}

	//##MK::DEBUG
	#pragma omp critical
	{
		cout << "Number of reindexed vertices " << fvert_reindexed.size() << endl;
		cout << "Number of unique float triplets " << verts_triplets.size() << endl;
		cout << "Number of naive float triplets " << debug << endl;
	}
}


void tess::build_tessportion0()
{
	//5.0 as target value based on empirical result by Rycroft
	//http://math.lbl.gov/voro++/examples/timing_test/
	p3i blocks = identify_blockpartitioning( ionportion.size(), Settings::TessellationPointsPerBlock, mywindow );
	//blockpartitioning is the essential trick in Voro++ to reduce spatial querying costs
	//too few blocks too many points to test against
	//too many blocks too much memory overhead and cache trashing (L1,L2 and not to forget "page" TLB cache)

	bool periodicbox = false;
	//container con(	mywindow.xmi - AABBINCLUSION_EPSILON , mywindow.xmx + AABBINCLUSION_EPSILON,
	//					mywindow.ymi - AABBINCLUSION_EPSILON , mywindow.ymx + AABBINCLUSION_EPSILON,
	container con(	mywindow.xmi, mywindow.xmx,
					mywindow.ymi, mywindow.ymx,
					mywindow.zmi, mywindow.zmx,
					blocks.x, blocks.y, blocks.z,
					periodicbox, periodicbox, periodicbox, 1);

	//c_loop_all cl(con);

	int i = 0;
	//MK::need to start i at zero only then
	//when we can use the cl.pid() Voro++ internal cell id as an index on localID2worldID for getting corresponding cell from global dataset
	//and as an index on ionportion for the local portion respectively
	for ( auto it = ionportion.begin(); it != ionportion.end(); ++it ) { //##MK::is the container order-preserving ?
		con.put(i, it->x, it->y, it->z );
		//localID2worldID.push_back( it->id ); //want to allow mapping of local cell index into consistent and unique global cell index
		/*
		cellmeta.push_back( cellstats() );
		wallmeta.push_back( wallstats() );
		cellaabb.push_back( aabb3d() );
		pvecmeta.push_back( "" );
		*/
		i++;
	}

	#pragma omp critical
	{
		cout << "Thread " << omp_get_thread_num() << " blockpartitioning " << blocks.x << ";" << blocks.y << ";" << blocks.z << endl;
		cout << "Thread " << omp_get_thread_num() << " mywindow\t\t" << mywindow << " ionportion.size() " << ionportion.size() << endl;
		cout << "Thread " << omp_get_thread_num() << " number of particles in con " << con.total_particles() << endl;
	}

	//temporaries per cell, will be reinitialized during looping over cells
	vector<double> first_v;
	vector<int> first_neigh;
	vector<int> first_f_vert;

	c_loop_all cl(con);
	voronoicell_neighbor c;
	short thrid = 1 + static_cast<short>(omp_get_thread_num());
	//##MK::less threads +- interpreted than numeric_limits<short>::max() and ::lowest(), adding +1 because 0 has no sign,
	//i.e. threads reported in Fortran indexing style

	//size_t prec = 0; //storage precision demand
	if (cl.start()) { //this is an incremental computing of the tessellation which holds at no point the entire tessellation
		do {
			if ( con.compute_cell(c, cl) == true ) {
				int clpid = cl.pid();
				//MK::utilize that voroxx container stores order preserving
				unsigned int worldid = ionportion.at(clpid).id;
				unsigned char typid = ( ionportion.at(clpid).iontype > UCHARMX ) ? UCHARMX : ionportion.at(clpid).iontype;
				short cstatus = ( ionportion.at(clpid).sgn == VALIDZONE_ION ) ? +thrid : -thrid;

#ifdef VALIDZONE_IONS_ONLY
				if ( ionportion.at(clpid).sgn == VALIDZONE_ION ) {
#endif
					//extract only those Voronoi cells that are not formed by guard points in the halo region
					mystats.CellsInside++;

					//gather boundary contact information about the cell
					double x, y, z;
					cl.pos(x, y, z);
					c.neighbors(first_neigh);
					c.face_vertices(first_f_vert);
					c.vertices(x, y, z, first_v);

					mystats.NumberOfFacetsProcessed += first_neigh.size();
					mystats.NumberOfVertsProcessed += first_v.size();

					//find me the array position on ionpp3 and the equally implicitly name ordering coarray ion2surf
					//from which to extract the distance of the material point to the dataset boundary
					unsigned int next_validion2surf_idx = surfdistance_ids->at(clpid);
					apt_real distance = F32MX; //dummy value, everybody taken
					if ( next_validion2surf_idx != DO_NOT_DEREFERENCE ) { //a distance value per computed exists
						//MK::remember for guardzone ions no distance value exists!
						distance = surfdistance_vals->at(next_validion2surf_idx);
					}

					//cout << "Distance " << setprecision(18) << distance << "\t\t" << SQR(Settings::SurfaceCellsCarvingRadius) << endl;

					if ( distance >= SQR(Settings::SurfaceCellsCarvingRadius) ) {
						//need axis-aligned bounding box for wallcontact correction if cell has compute domain as boundary
						//identify boundary contact configuration
						wallstats rvewall = identify_cell_wallcontact( first_neigh );

						if ( rvewall.any_touch == true ) {
							mystats.CellsWallcontact++;
							if ( rvewall.halo_insufficient == true ) //made boundary contact but halo guard was still sufficient to let cell
								mystats.CellsIncorrectGeometry++; //can happen if guard not large enough
						}

						//gather information about the computed potentially truncated Voronoi cell
						/*
						aabb3d first_hull = identify_cell_aabb3d( first_v );
						cellaabb.at(clpid) = first_hull;
						wallmeta.at(clpid) = rvewall;
						cellmeta.at(clpid) = cellstats( p3d(x, y, z), c.volume(), c.number_of_faces() );
						*/

						/*
						//characterize cell topology using pvectors
						string pvec = "";
						if ( truncated == false )
							pvec = identify_cell_pvector( first_neigh, first_f_vert );
						else
							pvec = identify_cell_pvector( second_neigh, second_f_vert );

						pvecmeta.at(id) = pvec;

						auto it = mypvectors.find( pvec );
						if ( it != mypvectors.end() ) { //adding count on an existent cell topology
							it->second++;
						}
						else { //adding a no yet encountered cell topology
							mypvectors.insert( pair<string,int>( pvec, 1) );
						}
						*/

						if ( i_store_tess_metainfo == true ) {
							//##MK::potential for visualizing only the incorrect cells here
							//identify largest integer value required to store topology references and unique vertex coordinate values
#ifdef VALIDZONE_IONS_ONLY
							if ( rvewall.any_touch == false ) {
#endif
								//collect heavy data, remember: only VALIDZONE_IONs are stored
								io_ion2type.push_back( typid );
								io_ion2cell.push_back( worldid );
								if ( i_store_tess_cellpos == true ) {
									io_cellposition.push_back( x );
									io_cellposition.push_back( y );
									io_cellposition.push_back( z );
								}
								io_vol.push_back( c.volume() ); //##MK::precision loss from 64 to 32bit
								io_nfaces.push_back( first_neigh.size() );

								unsigned char has_wall_contact = rvewall.any_touch == false ? NO : YES;
								io_wall.push_back( has_wall_contact );
								myprecisiondemand.n_cells += 1;

								cell_write_meta demands = cell_write_meta();
								demands = collect_cell_heavydata2( cstatus, worldid, first_neigh, first_f_vert, first_v );

								//if ( demands.prec >= myprecisiondemand.prec ) { myprecisiondemand.prec = demands.prec; }
								myprecisiondemand.n_topo += demands.n_topo;
								myprecisiondemand.n_geom += demands.n_geom;
#ifdef VALIDZONE_IONS_ONLY
							}
#endif
						}
					} //a cell deep enough in tip volume
					else {
						mystats.CellsEroded++;
					}
#ifdef VALIDZONE_IONS_ONLY
				}
				else { //reporting an invalid Voronoi cell, i.e. one that is formed by a guardzone point and hence as a duplicate in another thread
					mystats.CellsHalo++;
				}
#endif
			} //done reporting a valid Voronoi cell
			else { //reporting an unconstructable Voronoi cell
				mystats.CellsInconstructible++;
			}
		} while (cl.inc());
	}

	#pragma omp critical
	{
		cout << "ThreadID " << omp_get_thread_num() << " portion build successfully " << "\n";
		cout << "PrecisionDemand " << myprecisiondemand << "\n";
		cout << mystats << "\n";
		cout << "NumberOfFacetsProcessed " << mystats.NumberOfFacetsProcessed << "\n";
		cout << "NumberOfVertsProcessed " << mystats.NumberOfVertsProcessed << "\n";
		cout << endl;
		/*
		//##MK::DEBUG::Report pvector of local tessellation
		cout << "Total number of disjoint topologies " << mypvectors.size() << endl;
		for( auto it = mypvectors.begin(); it != mypvectors.end(); ++it ) {
			cout << it->first << "\t\t\t\t" << it->second << endl;
		}
		*/
	}
}



tessHdl::tessHdl()
{
	allprecisiondemand = cell_write_meta();
	//MK::activate or deactivate if desired...
	we_store_tess_metainfo = true;
	we_store_tess_cellpos = false;
	we_store_tess_topogeom = false;
	owner = NULL;
}


tessHdl::~tessHdl()
{
	for( size_t i = 0; i < local_tessellations.size(); ++i ) {
		delete local_tessellations.at(i);
	}
	//do not delete owner only backreference
}


void tessHdl::configure( const bool store_metainfo, const bool store_cellpositions, const bool store_topogeom )
{
	we_store_tess_metainfo = store_metainfo;
	we_store_tess_cellpos = store_cellpositions;
	we_store_tess_topogeom = store_topogeom;
}


int tessHdl::threaded_tessellating( vector<p3dm3> const & in, aabb3d const & cuboid,
		vector<zlim> const & zbounds, p3d const & halosize )
{
	//MK::generate one such tessHdl object per running MPI process if want to split hybrid manner tessellation construction
	#pragma omp parallel //input parameter shared by default, storetessellation read only so no problem if shared as well
	{
		double mytic = omp_get_wtime();
		int mt = omp_get_thread_num();
		int nt = omp_get_num_threads();

		#pragma omp master
		{
			for( int thr = 0; thr < nt; ++thr )
				local_tessellations.push_back(NULL);
			if ( we_store_tess_metainfo == true ||
					we_store_tess_cellpos == true ||
						we_store_tess_topogeom == true ) {
				for( int thr = 0; thr < nt; ++thr )
					iohelp.push_back( voro_io_bounds() ); //we dont know the bounds yet...
			}
		}
		#pragma omp barrier //necessary omp master has no barrier and memory needs to be filled before correct address written to it

		//build thread-local instances of tessellation object, using thread-local memory and first-touch and write into it
		tess* thrlocaltess = NULL;
		thrlocaltess = new tess; //try/catch
		thrlocaltess->i_store_tess_metainfo = we_store_tess_metainfo;
		thrlocaltess->i_store_tess_cellpos = we_store_tess_cellpos;
		thrlocaltess->i_store_tess_topogeom = we_store_tess_topogeom;

		//pull atoms into local objects based on bounds, as this writes to the memory of the tess object it does the first-touch
		thrlocaltess->set_zbounds( zbounds.at(mt) );
		thrlocaltess->set_halosize_p3d( halosize );
		thrlocaltess->set_haloinfo();

		thrlocaltess->pull_ionportion( in, cuboid );

		thrlocaltess->build_tessportion0();
		//##MK::pull a memory snapshot of program here because now heavy data still exist in memory!

		//necessary all data have to be processed before cooperatively instructing further tasks
		double mytoc = omp_get_wtime();

		#pragma omp critical
		{
			local_tessellations.at(mt) = thrlocaltess; //linking address of thread-local object data to tessHdl necessary at the latest here before cooperative stuff

			if ( we_store_tess_metainfo == true ) {
				//communicate storage demands
				//if ( thrlocaltess->myprecisiondemand.prec > allprecisiondemand.prec ) { allprecisiondemand.prec = thrlocaltess->myprecisiondemand.prec; }
				allprecisiondemand.n_cells += thrlocaltess->myprecisiondemand.n_cells;
				iohelp.at(mt).cell_n = thrlocaltess->myprecisiondemand.n_cells;

				if ( we_store_tess_topogeom == true ) {
					allprecisiondemand.n_topo += thrlocaltess->myprecisiondemand.n_topo;
					allprecisiondemand.n_geom += thrlocaltess->myprecisiondemand.n_geom;
					//communicate also io implicit array bound demands for writing data
					iohelp.at(mt).topo_n = thrlocaltess->myprecisiondemand.n_topo;
					iohelp.at(mt).geom_n = thrlocaltess->myprecisiondemand.n_geom;
					//MK::latter values will be translated into absolute position intervals on implicit 1d data arrays during I/O operation
				}
			}
			cout << "Thread " << mt << " building local tessellation took " << (mytoc-mytic) << " seconds" << endl;
		}
	}
	//##MK::consider in the future to do multi-threaded I/O still
	return 0;
}


/*
//##MK::modify to take into account that voro_io_info uses now size_t !
int tessHdl::report_tessellation_chunked()
{
	if ( we_store_tess_metainfo == false ) { //##MK::CHUNKED VERSION NEEDS MODIFICATION
		return 0;
	}

	double mytic = omp_get_wtime();
	//create I/O plan, i.e. process cumulated writing positions on the implicit meta, topology, and geometry arrays
	int nthr = local_tessellations.size(); //##MK::bounds check
	for( int thr = MASTER; thr < nthr; thr++ ) {
		iohelp.at(thr).cell_s = 0;		iohelp.at(thr).topo_s = 0;		iohelp.at(thr).geom_s = 0;
		for( int trailing = 0; trailing < thr; trailing++ ) {
			iohelp.at(thr).cell_s += iohelp.at(trailing).cell_n;
			iohelp.at(thr).topo_s += iohelp.at(trailing).topo_n;
			iohelp.at(thr).geom_s += iohelp.at(trailing).geom_n;
		}
		iohelp.at(thr).cell_e = iohelp.at(thr).cell_s + iohelp.at(thr).cell_n;
		iohelp.at(thr).topo_e = iohelp.at(thr).topo_s + iohelp.at(thr).topo_n;
		iohelp.at(thr).geom_e = iohelp.at(thr).geom_s + iohelp.at(thr).geom_n;
		//##MK::DEBUG REPORT DETAILED I/O write position plan for heavy data of every thread
		cout << "I/O plan thread " << thr << "\tc\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << endl;
		cout << "I/O plan thread " << thr << "\tt\t" << iohelp.at(thr).topo_s << "\t\t" << iohelp.at(thr).topo_e << endl;
		cout << "I/O plan thread " << thr << "\tg\t" << iohelp.at(thr).geom_s << "\t\t" << iohelp.at(thr).geom_e << endl;
	}

	//once the I/O plan exists, the master thread can start initializing the actual HDF5 H5 output file
	int status = 0;
	const string thish5fn = "PVOROXX.DummyFile.h5";
	myh5.reinitialize();
	status = myh5.create_file( thish5fn );
	status = myh5.create_group( thish5fn, "/VoronoiTess" );


	//store mapping of ion to its Voronoi cell, whose ions Voronoi cell is it?
	status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/Ion2CellMappingCounts", allprecisiondemand.n_cells );
	status = myh5.init_chunked_matrix_u32le( thish5fn, "/VoronoiTess/Ion2CellMappingValues", allprecisiondemand.n_cells, 1);
	for( int thr = 0; thr < nthr; thr++ ) {
		status = myh5.write_chunked_matrix_u32le( thish5fn, "/VoronoiTess/Ion2CellMappingValues", local_tessellations.at(thr)->io_ion2cell );
		local_tessellations.at(thr)->io_ion2cell = vector<unsigned int>(); //MK::vector clearing hack to really delete memory and not just clear i.e. have it remaining allocated...
	}
	status = myh5.reset_chunked_matrix_u32le_aftercompletion();

	//report distribution of volume
	status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/VolumeCounts", allprecisiondemand.n_cells );
	status = myh5.init_chunked_matrix_f64le( thish5fn, "/VoronoiTess/VolumeValues", allprecisiondemand.n_cells, 1);
	for( int thr = 0; thr < nthr; thr++ ) {
		status = myh5.write_chunked_matrix_f64le( thish5fn, "VoronoiTess/VolumeValues", local_tessellations.at(thr)->io_vol );
		local_tessellations.at(thr)->io_vol = vector<double>();
	}
	status = myh5.reset_chunked_matrix_f64le_aftercompletion();

	//report distribution of number of faces
	status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/NumberOfFacesCounts", allprecisiondemand.n_cells );
	status = myh5.init_chunked_matrix_u32le( thish5fn, "/VoronoiTess/NumberOfFacesValues", allprecisiondemand.n_cells, 1);
	for( int thr = 0; thr < nthr; thr++ ) {
		status = myh5.write_chunked_matrix_u32le( thish5fn, "/VoronoiTess/NumberOfFacesValues", local_tessellations.at(thr)->io_nfaces );
		local_tessellations.at(thr)->io_nfaces = vector<unsigned int>();
	}
	status = myh5.reset_chunked_matrix_u32le_aftercompletion();

	//report collection of n-polygons and their topology build the Voronoi cells
	//re-label cell local topology to global topology knowing now the topology and geometry of every cell
	unsigned int VertexOffset = 0;
	for( int thr = 0; thr < nthr; thr++ ) {
		tess* thisone = local_tessellations.at(thr);
		size_t CurrentCellID = 0;
		vector<unsigned int>& thesefaces = thisone->io_topo;
		size_t ni = thisone->io_topo.size();
		for( size_t i = 1; i < ni;   ) { //first value tells XDMF topology typ index of first local cell ignore this value
			unsigned int CurrentCellNumberOfFaces = thisone->io_info.at(CurrentCellID).nfacets;
			for( unsigned int f = 0; f < CurrentCellNumberOfFaces; f++ ) {
				unsigned int CurrentFacetNumberOfIndices = thesefaces.at(i); //##MK::unsigned int to size_t no problem
				i++;
				for( unsigned int k = 0; k < CurrentFacetNumberOfIndices; k++ ) {
					thesefaces.at(i+k) = thesefaces.at(i+k) + VertexOffset;
				}
				i = i + CurrentFacetNumberOfIndices + 1; //skip updated values skip XDMF type key
			}
//cout << "CurrentCellNumberOfFaces--->" << CurrentCellNumberOfFaces << "\t\t\t" << VertexOffset << endl;
			VertexOffset = VertexOffset + thisone->io_info.at(CurrentCellID).ngeom;
			CurrentCellID++;
		}
	}
	cout << "Threadlocal and celllocal topology indexing relabeled to global indexing" << endl;

	status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/XDMFCellTopologyCounts", allprecisiondemand.n_topo );
	status = myh5.init_chunked_matrix_u32le( thish5fn, "/VoronoiTess/XDMFCellTopologyValues", allprecisiondemand.n_topo, 1);
	for( int thr = 0; thr < nthr; thr++ ) {
		status = myh5.write_chunked_matrix_u32le( thish5fn, "/VoronoiTess/XDMFCellTopologyValues",
				local_tessellations.at(thr)->io_topo );
			local_tessellations.at(thr)->io_topo = vector<unsigned int>();
	}
	status = myh5.reset_chunked_matrix_u32le_aftercompletion();

	//report collection of unique vertex triplets for all n-polygons
	status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/XDMFCellVerticesTripletCounts", (allprecisiondemand.n_geom / 3) );
	status = myh5.init_chunked_matrix_f64le( thish5fn, "/VoronoiTess/XDMFCellVerticesValues", (allprecisiondemand.n_geom / 3), 3);
	for( int thr = 0; thr < nthr; thr++ ) {
		status = myh5.write_chunked_matrix_f64le( thish5fn, "VoronoiTess/XDMFCellVerticesValues",
				local_tessellations.at(thr)->io_geom );
		local_tessellations.at(thr)->io_geom = vector<double>();
	}
	status = myh5.reset_chunked_matrix_f64le_aftercompletion();

	//write XDMF metadata file for visualizating with e.g. Paraview or VisIt
	//#######MK::status = myxdmf.open_emptyfile( "PVOROXX.DummyFile.xdmf" );
	//#######MK::status = myxdmf.write_voronoitess_meta();

	double mytoc = omp_get_wtime();
	cout << "Reporting of tessellation via XDMF/HDF5 file took " << (mytoc-mytic) << " seconds" << endl;
	return 0;
}
*/


/*
int tessHdl::report_tessellation_contiguous_atonce()
{
	//we may not need to write out topogeom but if we report at least we need metainfo
	if ( we_store_tess_metainfo == false )
		return 0;

	double mytic = omp_get_wtime();
	//create I/O plan, i.e. process cumulated writing positions on the implicit meta, topology, and geometry arrays
	size_t nthr = local_tessellations.size(); //##MK::bounds check
	for( size_t thr = MASTER; thr < nthr; thr++ ) {
		//if we store tess somehow we need at least these metadata
		iohelp.at(thr).cell_s = 0;
		for( int trailing = 0; trailing < thr; trailing++ )
			iohelp.at(thr).cell_s += iohelp.at(trailing).cell_n;
		iohelp.at(thr).cell_e = iohelp.at(thr).cell_s + iohelp.at(thr).cell_n;
//cout << "I/O plan thread " << thr << "\tc\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << endl;

		if ( we_store_tess_topogeom == true ) {
			iohelp.at(thr).topo_s = 0;		iohelp.at(thr).geom_s = 0;
			for( int trailing = 0; trailing < thr; trailing++ ) {
				iohelp.at(thr).cell_s += iohelp.at(trailing).cell_n;
				iohelp.at(thr).topo_s += iohelp.at(trailing).topo_n;
				iohelp.at(thr).geom_s += iohelp.at(trailing).geom_n;
			}
			iohelp.at(thr).topo_e = iohelp.at(thr).topo_s + iohelp.at(thr).topo_n;
			iohelp.at(thr).geom_e = iohelp.at(thr).geom_s + iohelp.at(thr).geom_n;
		}
//##MK::DEBUG REPORT DETAILED I/O write position plan for heavy data of every thread
//cout << "I/O plan thread " << thr << "\tt\t" << iohelp.at(thr).topo_s << "\t\t" << iohelp.at(thr).topo_e << endl;
//cout << "I/O plan thread " << thr << "\tg\t" << iohelp.at(thr).geom_s << "\t\t" << iohelp.at(thr).geom_e << endl;
	}
	reporting( "Voro output I/O plan created" );

	//////////////////////////////////
	//fuse all results on the master incrementally delete copies on threads ##MK::could be made even with finer granularity..
	//info
	vector<voro_io_info> & info_dump_master = local_tessellations.at(MASTER)->io_info;
	for ( size_t thr = MASTER+1; thr < nthr; thr++ ) {
		vector<voro_io_info> & info_pull_thr = local_tessellations.at(thr)->io_info;
		for( auto it = info_pull_thr.begin(); it != info_pull_thr.end(); ++it )
			info_dump_master.push_back( *it );
		info_pull_thr = vector<voro_io_info>();
	}
	reporting( "Voro_io_info pulled over" );
	//ion2cell mapping
	vector<unsigned int> & i2c_dump_master = local_tessellations.at(MASTER)->io_ion2cell;
	for( size_t thr = MASTER+1; thr < nthr; thr++ ) {
		vector<unsigned int> & i2c_pull_thr = local_tessellations.at(thr)->io_ion2cell;
		for( auto it = i2c_pull_thr.begin(); it != i2c_pull_thr.end(); ++it )
			i2c_dump_master.push_back( *it );
		i2c_pull_thr = vector<unsigned int>();
	}
	reporting( "Voro ion2cell pulled over" );
	//ion2type mapping
	vector<unsigned char> & i2t_dump_master = local_tessellations.at(MASTER)->io_ion2type;
	for( size_t thr = MASTER+1; thr < nthr; thr++ ) {
		vector<unsigned char> & i2t_pull_thr = local_tessellations.at(thr)->io_ion2type;
		for( auto it = i2t_pull_thr.begin(); it != i2t_pull_thr.end(); ++it )
			i2t_dump_master.push_back( *it );
		i2t_pull_thr = vector<unsigned char>();
	}
	reporting( "Voro ion2type pulled over" );
	//volume
	vector<float> & vol_dump_master = local_tessellations.at(MASTER)->io_vol;
	for( size_t thr = MASTER+1; thr < nthr; thr++ ) {
		vector<float> & vol_pull_thr = local_tessellations.at(thr)->io_vol;
		for( auto it = vol_pull_thr.begin(); it != vol_pull_thr.end(); ++it )
			vol_dump_master.push_back( *it );
		vol_pull_thr = vector<float>();
	}
	reporting( "Voro volume pulled over" );
	//number of faces
	vector<unsigned int> & nfaces_dump_master = local_tessellations.at(MASTER)->io_nfaces;
	for( size_t thr = MASTER+1; thr < nthr; thr++ ) {
		vector<unsigned int> & nfaces_pull_thr = local_tessellations.at(thr)->io_nfaces;
		for( auto it = nfaces_pull_thr.begin(); it != nfaces_pull_thr.end(); ++it )
			nfaces_dump_master.push_back( *it );
		nfaces_pull_thr = vector<unsigned int>();
	}
	reporting( "Voro number of faces pulled over" );
	//wallcontact
	vector<unsigned char> & wallinfo_dump_master = local_tessellations.at(MASTER)->io_wall;
	for( size_t thr = MASTER+1; thr < nthr; thr++ ) {
		vector<unsigned char> & wallinfo_pull_thr = local_tessellations.at(thr)->io_wall;
		for( auto it = wallinfo_pull_thr.begin(); it != wallinfo_pull_thr.end(); ++it )
			wallinfo_dump_master.push_back( *it );
		wallinfo_pull_thr = vector<unsigned char>();
	}
	reporting( "Voro wall info pulled over" );

	if( we_store_tess_topogeom == true ) {
		//topo
		vector<size_t> & topo_dump_master = local_tessellations.at(MASTER)->io_topo;
		for ( size_t thr = MASTER+1; thr < nthr; thr++ ) {
			vector<size_t> & topo_pull_thr = local_tessellations.at(thr)->io_topo;
			for( auto it = topo_pull_thr.begin(); it != topo_pull_thr.end(); ++it )
				topo_dump_master.push_back( *it );
			topo_pull_thr = vector<size_t>();
		}
		//re-label cell local topology to global topology knowing now the topology and geometry of every cell
		size_t VertexOffset = 0;
		size_t CurrentCellID = 0;
		vector<voro_io_info>& theseinfo = local_tessellations.at(MASTER)->io_info;
		vector<size_t>& thesefaces = local_tessellations.at(MASTER)->io_topo;
		size_t ni = thesefaces.size();
		for( size_t i = 1; i < ni;   ) { //first value tells XDMF topology typ index of first local cell ignore this value
			size_t CurrentCellNumberOfFaces = theseinfo.at(CurrentCellID).nfacets;
			for( size_t f = 0; f < CurrentCellNumberOfFaces; f++ ) {
				size_t CurrentFacetNumberOfIndices = thesefaces.at(i);
				i++;
				for( size_t k = 0; k < CurrentFacetNumberOfIndices; k++ ) {
					thesefaces.at(i+k) = thesefaces.at(i+k) + VertexOffset;
				}
				i = i + CurrentFacetNumberOfIndices + 1; //skip updated values skip XDMF type key
			}
	//cout << "CurrentCellNumberOfFaces--->" << CurrentCellNumberOfFaces << "\t\t\t" << VertexOffset << endl;
			VertexOffset = VertexOffset + theseinfo.at(CurrentCellID).ngeom;
			CurrentCellID++;
		}
		cout << "Threadlocal and celllocal topology indexing relabeled to global indexing" << endl;

		//cell geometry
		vector<float> & geom_dump_master = local_tessellations.at(MASTER)->io_geom;
		for ( size_t thr = MASTER+1; thr < nthr; thr++ ) {
			vector<float> & geom_pull_thr = local_tessellations.at(thr)->io_geom;
			for( auto it = geom_pull_thr.begin(); it != geom_pull_thr.end(); ++it )
				geom_dump_master.push_back( *it );
			geom_pull_thr = vector<float>();
		}
	}

	//once the I/O plan exists, the master thread can start initializing the actual HDF5 H5 output file
	int status = 0;
	const string thish5fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.h5";
	myh5.reinitialize();
	status = myh5.create_file( thish5fn );
	status = myh5.create_group( thish5fn, "/VoronoiTess" );
	status = myh5.create_group( thish5fn, "/VoronoiTess/DescrStats");

	if ( we_store_tess_topogeom == true ) {
		status = myh5.create_group( thish5fn, "/VoronoiTess/CellGeometry");
	}

	//store mapping of ion to its Voronoi cell, whose ions Voronoi cell is it?
	//ion2cell
	status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/DescrStats/NumberOfCells", allprecisiondemand.n_cells );
	//status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/DescrStats/Ion2CellMappingCounts", allprecisiondemand.n_cells );
	status = myh5.write_contiguous_matrix_u32le_atonce( thish5fn, "/VoronoiTess/DescrStats/Ion2CellMappingValues",
			allprecisiondemand.n_cells, 1, local_tessellations.at(MASTER)->io_ion2cell );
	//ion2type
	status = myh5.write_contiguous_matrix_u8le_atonce( thish5fn, "/VoronoiTess/DescrStats/Ion2Iontype",
			allprecisiondemand.n_cells, 1, local_tessellations.at(MASTER)->io_ion2type );
	//volume
	status = myh5.write_contiguous_matrix_f32le_atonce( thish5fn, "/VoronoiTess/DescrStats/VolumeValues",
			allprecisiondemand.n_cells, 1, local_tessellations.at(MASTER)->io_vol );
	//nfaces
	status = myh5.write_contiguous_matrix_u32le_atonce( thish5fn, "/VoronoiTess/DescrStats/NumberOfFacesValues",
			allprecisiondemand.n_cells, 1, local_tessellations.at(MASTER)->io_nfaces );
	//wallinfo
	status = myh5.write_contiguous_matrix_u8le_atonce( thish5fn, "/VoronoiTess/DescrStats/BoxContact",
			allprecisiondemand.n_cells, 1, local_tessellations.at(MASTER)->io_wall );
	//topo i.e. collection of n-polygons and their topology build the Voronoi cells
	size_t NumberOfFacets = 0;
	vector<unsigned int> & nfaces_pull_master = local_tessellations.at(MASTER)->io_nfaces;
	for( auto it = nfaces_pull_master.begin(); it != nfaces_pull_master.end(); ++it )
		NumberOfFacets += *it;
	cout << "Number of facets in reported in total " << NumberOfFacets << endl;
	status = myh5.create_scalar_u64le( thish5fn, "/VoronoiTess/DescrStats/NumberOfFacetsTotal", NumberOfFacets );

	if ( we_store_tess_topogeom == true ) {
		status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFTopologyNumberOfElements", NumberOfFacets );
		status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFTopologyDimensions", allprecisiondemand.n_topo );
		status = myh5.write_contiguous_matrix_u64le_atonce( thish5fn, "/VoronoiTess/CellGeometry/XDMFTopologyValues",
				allprecisiondemand.n_topo, 1, local_tessellations.at(MASTER)->io_topo );

		//geom i.e. collection of unique vertex triplets for all n-polygons
		status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFxyzTripletCounts", (allprecisiondemand.n_geom / 3) );
		status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFxyzDimensions", allprecisiondemand.n_geom );
		status = myh5.write_contiguous_matrix_f32le_atonce( thish5fn, "/VoronoiTess/CellGeometry/XDMFxyzValues",
				(allprecisiondemand.n_geom / 3), 3, local_tessellations.at(MASTER)->io_geom );

		//write XDMF metadata file for visualizating with e.g. Paraview or VisIt
		const string thisxmlfn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.xdmf";
		status = myxdmf.create_voronoicell_vis_file( thisxmlfn, NumberOfFacets,
				allprecisiondemand.n_topo, allprecisiondemand.n_geom, 0, thish5fn );
	}

	double mytoc = omp_get_wtime();
	cout << "Reporting of tessellation via XDMF/HDF5 file took " << (mytoc-mytic) << " seconds" << endl;
	return 0;
}
*/


int tessHdl::report_tessellation_hyperslab_perthread()
{
	//we may not need to write out topogeom but if we report at least we need metainfo
	if ( we_store_tess_metainfo == false )
		return 0;

	double mytic = omp_get_wtime();
	//create I/O plan, i.e. process cumulated writing positions on the implicit meta, topology, and geometry arrays
	size_t nthr = local_tessellations.size(); //##MK::bounds check
	for( size_t thr = MASTER; thr < nthr; thr++ ) {
		//if we store tess somehow we need at least these metadata
		iohelp.at(thr).cell_s = 0;
		for( size_t trailing = 0; trailing < thr; trailing++ ) {
			iohelp.at(thr).cell_s += iohelp.at(trailing).cell_n;
		}
		iohelp.at(thr).cell_e = iohelp.at(thr).cell_s + iohelp.at(thr).cell_n;
		if ( we_store_tess_topogeom == true ) {
			iohelp.at(thr).topo_s = 0;		iohelp.at(thr).geom_s = 0;
			for( int trailing = 0; trailing < thr; trailing++ ) {
				iohelp.at(thr).topo_s += iohelp.at(trailing).topo_n;
				iohelp.at(thr).geom_s += iohelp.at(trailing).geom_n;
			}
			iohelp.at(thr).topo_e = iohelp.at(thr).topo_s + iohelp.at(thr).topo_n;
			iohelp.at(thr).geom_e = iohelp.at(thr).geom_s + iohelp.at(thr).geom_n;
		}
//##MK::DEBUG REPORT DETAILED I/O write position plan for heavy data of every thread
//cout << "I/O plan thread " << thr << "\tc\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << endl;
//cout << "I/O plan thread " << thr << "\tt\t" << iohelp.at(thr).topo_s << "\t\t" << iohelp.at(thr).topo_e << endl;
//cout << "I/O plan thread " << thr << "\tg\t" << iohelp.at(thr).geom_s << "\t\t" << iohelp.at(thr).geom_e << endl;
	}
	reporting( "Voro output I/O plan created" );

	//once the I/O plan exists, the master thread can start initializing the actual HDF5 H5 output file and subgroups
	int status = 0;
	const string thish5fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.h5";
	myh5.reinitialize();
	status = myh5.create_file( thish5fn );
	status = myh5.create_group( thish5fn, PARAPROBE_VOLTESS );
	status = myh5.create_group( thish5fn, PARAPROBE_VOLTESS_DESCRSTATS);

	//MK::do not fuse results on the master
	//generate H5 datasets, fill incrementally with threadlocal buffer inplace subsequently
	//##MK::add PARAPROBE_VOLTESS_DESCRSTATS_NCELLS, "/VoronoiTess/DescrStats/NumberOfCells"
	h5iometa ifo = h5iometa();
	h5offsets offs = h5offsets();

	//ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_I2CELL, allprecisiondemand.n_cells, 1 );
	//status = myh5.create_contiguous_matrix_u32le( ifo );
	ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_I2TYPE, allprecisiondemand.n_cells, 1 );
	status = myh5.create_contiguous_matrix_u8le( ifo );
	if ( we_store_tess_cellpos == true ) {
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_CELLPOS, allprecisiondemand.n_cells, 3 );
		status = myh5.create_contiguous_matrix_f32le( ifo );
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_THREADID, allprecisiondemand.n_cells, 1 );
		status = myh5.create_contiguous_matrix_u8le( ifo );
	}
	ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_VOL, allprecisiondemand.n_cells, 1 );
	status = myh5.create_contiguous_matrix_f32le( ifo );
	ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_NFACES, allprecisiondemand.n_cells, 1 );
	status = myh5.create_contiguous_matrix_u32le( ifo );
	//ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_BND, allprecisiondemand.n_cells, 1 );
	//status = myh5.create_contiguous_matrix_u8le( ifo );
	//##MK::add PARAPROBE_VOLTESS_DESCRSTATS_NFTOTAL, "/VoronoiTess/DescrStats/NumberOfFacetsTotal"
	reporting( "Voro output I/O corresponding HDF5 groups and datasets initialized" );

	//in case of storing topology we need to know how many facets we have in total
	//topo i.e. collection of n-polygons and their topology build the Voronoi cells
	string mess = "";
	size_t NumberOfFacets = 0;

	//now that the contiguous datasets are initialize write threadlocal results knowing the offsets
	for(size_t thr = MASTER; thr < nthr; thr++ ) {
		/*
		vector<unsigned int> & u32thrbuffer = local_tessellations.at(thr)->io_ion2cell;
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_I2CELL, allprecisiondemand.n_cells, 1 );
		offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 1, u32thrbuffer.size(), 1);
		status = myh5.write_contiguous_matrix_u32le_hyperslab( ifo, offs, u32thrbuffer );
		u32thrbuffer = vector<unsigned int>();
		mess = "Thread " + to_string(thr) + " Voro ion2cell written"; reporting( mess );
		*/
		//ion2cell mapping
		local_tessellations.at(thr)->io_ion2cell = vector<unsigned int>();
		//ion2type mapping
		vector<unsigned char> & u8thrbuffer = local_tessellations.at(thr)->io_ion2type;
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_I2TYPE, allprecisiondemand.n_cells, 1 );
		offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 1, u8thrbuffer.size(), 1);
		cout << "cell_s/cell_e/cell_n/u8thrbuffer.size()\t\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << "\t\t" << iohelp.at(thr).cell_n << "\t\t" << u8thrbuffer.size() << endl;
		status = myh5.write_contiguous_matrix_u8le_hyperslab( ifo, offs, u8thrbuffer );
		cout << "Ion2type " << status << endl;
		u8thrbuffer = vector<unsigned char>();
		mess = "Thread " + to_string(thr) + " Voro ion2type written";
		reporting( mess );
		if ( we_store_tess_cellpos == true ) {
			//cell positions
			vector<float> & f32thrbuffer = local_tessellations.at(thr)->io_cellposition;
			ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_CELLPOS, allprecisiondemand.n_cells, 3);
			offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 3, (f32thrbuffer.size() / 3), 3);
			cout << "cell_s/cell_e/cell_n/f32thrbuffer.size()\t\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << "\t\t" << iohelp.at(thr).cell_n << "\t\t" << f32thrbuffer.size() << endl;
			status = myh5.write_contiguous_matrix_f32le_hyperslab( ifo, offs, f32thrbuffer );
			cout << "Cellpos " << status << endl;
			f32thrbuffer = vector<float>();
			mess = "Thread " + to_string(thr) + " Voro cell positions written";
			reporting( mess );
			//thread ID
			//size_t MaxThreadID = UCHARMX;
			unsigned char AttrVal = (thr < static_cast<size_t>(UCHARMX)) ? static_cast<unsigned char>(thr) : static_cast<unsigned char>(UCHARMX);
			vector<unsigned char> u8buf( iohelp.at(thr).cell_n, AttrVal );
			ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_THREADID, allprecisiondemand.n_cells, 1);
			offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 1, u8buf.size(), 1 );
			cout << "cell_s/cell_e/cell_n/u8buf.size()\t\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << "\t\t" << iohelp.at(thr).cell_n << "\t\t" << u8buf.size() << endl;
			status = myh5.write_contiguous_matrix_u8le_hyperslab( ifo, offs, u8buf );
			cout << "ThreadID " << status << endl;
			u8buf = vector<unsigned char>();
			mess = "Thread " + to_string(thr) + " Voro threads IDs written";
			reporting( mess );
		}
		//volume
		vector<float> & f32thrbuffer = local_tessellations.at(thr)->io_vol;
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_VOL, allprecisiondemand.n_cells, 1 );
		offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 1, f32thrbuffer.size(), 1);
		cout << "cell_s/cell_e/cell_n/f32thrbuffer.size()\t\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << "\t\t" << iohelp.at(thr).cell_n << "\t\t" << f32thrbuffer.size() << endl;
		status = myh5.write_contiguous_matrix_f32le_hyperslab( ifo, offs, f32thrbuffer );
		cout << "Vol " << status << endl;
		f32thrbuffer = vector<float>();
		mess = "Thread " + to_string(thr) + " Voro volume written";
		reporting( mess );
		//number of faces
		vector<unsigned int> & u32thrbuffer = local_tessellations.at(thr)->io_nfaces;
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_NFACES, allprecisiondemand.n_cells, 1 );
		offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 1, u32thrbuffer.size(), 1);
		cout << "cell_s/cell_e/cell_n/u32thrbuffer.size()\t\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << "\t\t" << iohelp.at(thr).cell_n << "\t\t" << u32thrbuffer.size() << endl;
		status = myh5.write_contiguous_matrix_u32le_hyperslab( ifo, offs, u32thrbuffer );
		cout << "Nfa " << status << endl;
		//count how many faces in total before deleting buffer content
		if( we_store_tess_topogeom == true ) {
			for( auto it = u32thrbuffer.begin(); it != u32thrbuffer.end(); ++it )
				NumberOfFacets += *it;
		}
		u32thrbuffer = vector<unsigned int>();
		mess = "Thread " + to_string(thr) + " Voro number of faces written";
		reporting( mess );
		//wall contact
		/*
		u8thrbuffer = local_tessellations.at(thr)->io_wall;
		ifo = h5iometa( PARAPROBE_VOLTESS_DESCRSTATS_BND, allprecisiondemand.n_cells, 1 );
		offs = h5offsets( iohelp.at(thr).cell_s, iohelp.at(thr).cell_e, 0, 1, u8thrbuffer.size(), 1);
		status = myh5.write_contiguous_matrix_u8le_hyperslab( ifo, offs, u8thrbuffer );
		u8thrbuffer = vector<unsigned char>();
		mess = "Thread " + to_string(thr) + " Voro wallcontact written";
		reporting( mess );
		*/
		local_tessellations.at(thr)->io_wall = vector<unsigned char>();
	}

	if ( we_store_tess_cellpos == true ) {
		const string thisxmlfn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroVolume.xdmf";
		status = myxdmf.create_voronoicell_vol_file( thisxmlfn, allprecisiondemand.n_cells, thish5fn );
	}

	if( we_store_tess_topogeom == true ) {
		cout << "Number of facets in reported in total " << NumberOfFacets << endl;

		status = myh5.create_group( thish5fn, PARAPROBE_VOLTESS_CELLS );
		ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_TOPOLOGY, allprecisiondemand.n_topo, 1 );
		status = myh5.create_contiguous_matrix_u64le( ifo );
		ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_GEOMETRY, (allprecisiondemand.n_geom / 3), 3 );
		status = myh5.create_contiguous_matrix_f32le( ifo );

		//MK::for topology knowing only the threadlocal entry positions and offsets
		//is insufficient we need to relabel threadlocal vertex IDs to process global vertex IDs for this pull first offsets counts per thread
		vector<voro_io_info> & info_dump_master = local_tessellations.at(MASTER)->io_info;
		for ( size_t thr = MASTER+1; thr < nthr; thr++ ) {
			vector<voro_io_info> & info_pull_thr = local_tessellations.at(thr)->io_info;
			for( auto it = info_pull_thr.begin(); it != info_pull_thr.end(); ++it )
				info_dump_master.push_back( *it );
			info_pull_thr = vector<voro_io_info>();
		}
		reporting( "Voro_io_info pulled from threads on master thread" );

		//re-label cell local topology to global topology knowing now the topology and geometry of every cell
		size_t VertexOffset = 0;
		size_t CurrentCellID = 0;
		vector<voro_io_info> & theseinfo = local_tessellations.at(MASTER)->io_info;
		tess* mypart = NULL;
		for( size_t thr = MASTER; thr < nthr; thr++ ) {
			//[MASTER,nthr) because IDs needing relabeling for all threads
			//modify threadlocal facet IDs knowing process global context
			//##MK::error
			mypart = local_tessellations.at(thr);
			vector<size_t>& thesefaces = mypart->io_topo;
			size_t ni = thesefaces.size();
			for( size_t i = 1; i < ni;   ) { //first value tells XDMF topology typ index of first local cell ignore this value
				size_t CurrentCellNumberOfFaces = theseinfo.at(CurrentCellID).nfacets;
				for( size_t f = 0; f < CurrentCellNumberOfFaces; f++ ) {
					size_t CurrentFacetNumberOfIndices = thesefaces.at(i);
					i++;
					for( size_t k = 0; k < CurrentFacetNumberOfIndices; k++ ) {
						thesefaces.at(i+k) = thesefaces.at(i+k) + VertexOffset;
					}
					i = i + CurrentFacetNumberOfIndices + 1; //skip updated values skip XDMF type key
				}
				//cout << "CurrentCellNumberOfFaces--->" << CurrentCellNumberOfFaces << "\t\t\t" << VertexOffset << endl;
				VertexOffset = VertexOffset + theseinfo.at(CurrentCellID).ngeom;
				CurrentCellID++;
			}
			mess = "Thread " + to_string(thr) + " threadlocal topology relabeled according to process global context";
			reporting( mess );

			//write topology vertex reference IDs to file
			ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_TOPOLOGY, allprecisiondemand.n_topo, 1 );
			offs = h5offsets( iohelp.at(thr).topo_s, iohelp.at(thr).topo_e, 0, 1, allprecisiondemand.n_topo, 1 );
			cout << "topo_s/topo_e/io_topo.size()\t\t" << iohelp.at(thr).topo_s << "\t\t" << iohelp.at(thr).topo_e << "\t\t" << mypart->io_topo.size() << endl;
			status = myh5.write_contiguous_matrix_u64le_hyperslab( ifo, offs, mypart->io_topo );
			cout << "Topo\t\t" << status << endl;
			mypart->io_topo = vector<size_t>();
			mess = "Thread " + to_string(thr) + " threadlocal topology written";
			reporting( mess );

			//cell geometry
			ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_GEOMETRY, (allprecisiondemand.n_geom / 3), 3);
			offs = h5offsets( (iohelp.at(thr).geom_s / 3), (iohelp.at(thr).geom_e / 3), 0, 3, (allprecisiondemand.n_geom / 3), 3 ); //##MK::-1,-1 indicating not used
			cout << "geom_s/geom_e/io_geom.size()\t\t" << iohelp.at(thr).geom_s << "\t\t" << iohelp.at(thr).geom_e << "\t\t" << mypart->io_geom.size() << endl;
			status = myh5.write_contiguous_matrix_f32le_hyperslab( ifo, offs, mypart->io_geom );
			cout << "Geom\t\t" << status << endl;
			//geom i.e. collection of unique vertex triplets for all n-polygons
			mypart->io_geom = vector<float>();
			mess = "Thread " + to_string(thr) + " threadlocal vertex positions written";
			reporting( mess );

			//##MK::exemplary attribute values, one value per polygon facets, not per cell
		} //write next thread's local results

		//finally write XDMF metadata file for assisting visualization of heavy data with e.g. Paraview or VisIt
		const string thisxmlfn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.xdmf";
		status = myxdmf.create_voronoicell_vis_file( thisxmlfn, NumberOfFacets,
				allprecisiondemand.n_topo, allprecisiondemand.n_geom, allprecisiondemand.n_cells, thish5fn );

		//##MK::so far not written into H5 file different than for *atonce version of this output routine
		//status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFTopologyNumberOfElements", NumberOfFacets );
		//status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFTopologyDimensions", allprecisiondemand.n_topo );
		//status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFxyzTripletCounts", (allprecisiondemand.n_geom / 3) );
		//status = myh5.create_scalar_u64le( thish5fn, "VoronoiTess/CellGeometry/XDMFxyzDimensions", allprecisiondemand.n_geom );
	}

	double mytoc = omp_get_wtime();
	cout << "Reporting of tessellation via XDMF/HDF5 file took " << (mytoc-mytic) << " seconds" << endl;
	return 0;
}


int tessHdl::report_tessellation_hyperslab_onlythethread()
{
	//we may not need to write out topogeom but if we report at least we need metainfo
	if ( we_store_tess_metainfo == false || we_store_tess_topogeom == false )
		return 0;

	double mytic = omp_get_wtime();
	//write individual tessellation of the threads into independent files to be able to visualize
	//results of them individually including the guardzones
	size_t nthr = local_tessellations.size();
	string mess = "";
	int status = 0;
	h5iometa ifo = h5iometa();
	h5offsets offs = h5offsets();
	for( size_t thr = MASTER; thr < nthr; thr++ ) { //if we store tess somehow we need at least I/O bounds but not accumulated because individually
		voro_io_bounds thrioinfo = voro_io_bounds(
			iohelp.at(thr).cell_n, iohelp.at(thr).cell_s, (iohelp.at(thr).cell_s + iohelp.at(thr).cell_n),
			iohelp.at(thr).topo_n, iohelp.at(thr).topo_s, (iohelp.at(thr).topo_s + iohelp.at(thr).topo_n),
			iohelp.at(thr).geom_n, iohelp.at(thr).geom_s, (iohelp.at(thr).geom_s + iohelp.at(thr).geom_n)    );
//##MK::DEBUG REPORT DETAILED I/O write position plan for heavy data of every thread
//cout << "I/O plan thread " << thr << "\tc\t" << iohelp.at(thr).cell_s << "\t\t" << iohelp.at(thr).cell_e << endl;
//cout << "I/O plan thread " << thr << "\tt\t" << iohelp.at(thr).topo_s << "\t\t" << iohelp.at(thr).topo_e << endl;
//cout << "I/O plan thread " << thr << "\tg\t" << iohelp.at(thr).geom_s << "\t\t" << iohelp.at(thr).geom_e << endl;
		mess = "Voro output I/O plan for thread " + to_string(thr) + " created";
		reporting(mess);

		tess* mypart = local_tessellations.at(thr);
		size_t NumberOfFacets = 0;
		vector<unsigned int> & u32thrbuffer = mypart->io_nfaces;
		for( auto it = u32thrbuffer.begin(); it != u32thrbuffer.end(); ++it )
			NumberOfFacets += *it;

		//##MK::relabeling of cell topology required because vertex indices are per cell
		size_t VertexOffset = 0;
		size_t CurrentCellID = 0;
		vector<voro_io_info> & theseinfo = mypart->io_info;
		vector<size_t> & thesefaces = mypart->io_topo;
		size_t ni = thesefaces.size();
		for( size_t i = 1; i < ni;   ) { //first value tells XDMF topology typ index of first local cell ignore this value
			size_t CurrentCellNumberOfFaces = theseinfo.at(CurrentCellID).nfacets;
			for( size_t f = 0; f < CurrentCellNumberOfFaces; f++ ) {
				size_t CurrentFacetNumberOfIndices = thesefaces.at(i);
				i++;
				for( size_t k = 0; k < CurrentFacetNumberOfIndices; k++ ) {
					thesefaces.at(i+k) = thesefaces.at(i+k) + VertexOffset;
				}
				i = i + CurrentFacetNumberOfIndices + 1; //skip updated values skip XDMF type key
			}
			//cout << "CurrentCellNumberOfFaces--->" << CurrentCellNumberOfFaces << "\t\t\t" << VertexOffset << endl;
			VertexOffset = VertexOffset + theseinfo.at(CurrentCellID).ngeom;
			CurrentCellID++;
		}

		mess = "Voro output I/O number of facets " + to_string(thr) + " is " + to_string(NumberOfFacets);
		reporting(mess);

		h5Hdl thrh5 = h5Hdl();
		const string thrh5fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.Thread." + to_string(thr) + ".h5";
		cout << thrh5fn << endl;
		status = thrh5.create_file( thrh5fn );
		cout << "Creation " << status << endl;
		status = thrh5.create_group( thrh5fn, PARAPROBE_VOLTESS );
		cout << "VoroTessGrp " << status << endl;
		status = thrh5.create_group( thrh5fn, PARAPROBE_VOLTESS_CELLS );
		cout << "VoroTessCellGrp " << status << endl;
		//##MK::so far no DescrStats only topology geometry and cell facet polygon attribute data namely the threadID

		//write topology vertex reference IDs to file
		ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_TOPOLOGY, thrioinfo.topo_n, 1 );
		status = thrh5.create_contiguous_matrix_u64le( ifo );
		cout << "Generate topo " << status << endl;
		offs = h5offsets( thrioinfo.topo_s, thrioinfo.topo_e, 0, 1, thrioinfo.topo_n, 1 );
		cout << "topo_s/topo_e/io_topo.size()\t\t" << thrioinfo.topo_s << "\t\t" << thrioinfo.topo_e << "\t\t" << thrioinfo.topo_n << "\t\t" << local_tessellations.at(thr)->io_topo.size() << endl;
		status = thrh5.write_contiguous_matrix_u64le_hyperslab( ifo, offs, local_tessellations.at(thr)->io_topo );
		cout << "Write topo " << status << endl;
		local_tessellations.at(thr)->io_topo = vector<size_t>();

		//cell geometry
		ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_GEOMETRY, (thrioinfo.geom_n / 3), 3 );
		status = thrh5.create_contiguous_matrix_f32le( ifo );
		cout << "Generate geom " << status << endl;
		offs = h5offsets( (thrioinfo.geom_s / 3), (thrioinfo.geom_e / 3), 0, 3, (thrioinfo.geom_n / 3), 3 );
		cout << "geom_s/geom_e/io_geom.size()\t\t" << thrioinfo.geom_s << "\t\t" << thrioinfo.geom_e << "\t\t" << thrioinfo.geom_n << "\t\t" << local_tessellations.at(thr)->io_geom.size() << endl;
		status = thrh5.write_contiguous_matrix_f32le_hyperslab( ifo, offs, local_tessellations.at(thr)->io_geom );
		cout << "Write geom " << status << endl;
		local_tessellations.at(thr)->io_geom = vector<float>();

		//save scalar halo attribute data per cell facet polygon positive thread IDs VALIDZONE_ION, negative thread IDs GUARDZONE_ION
#ifndef VALIDZONE_IONS_ONLY
		size_t halo_n = local_tessellations.at(thr)->io_halo.size();
		ifo = h5iometa( PARAPROBE_VOLTESS_CELLS_THREADIDATTR, halo_n, 1 );
		status = thrh5.create_contiguous_matrix_i16le( ifo );
		cout << "Generate attr " << status << endl;
		offs = h5offsets( 0, halo_n, 0, 1, halo_n, 1);
		cout << "halo_s/halo_e/halo_n/io_halo.size()\t\t" << 0 << "\t\t" << halo_n << "\t\t" << halo_n << "\t\t" << local_tessellations.at(thr)->io_halo.size() << endl;
		status = thrh5.write_contiguous_matrix_i16le_hyperslab( ifo, offs, local_tessellations.at(thr)->io_halo );
		cout << "Write halo " << status << endl;
		local_tessellations.at(thr)->io_halo = vector<short>();
#endif

		mess = "Thread " + to_string(thr) + " threadlocal data written";
		reporting( mess );

		//finally write XDMF metadata file for assisting visualization of heavy data with e.g. Paraview or VisIt
		const string thisxmlfn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.Thread." + to_string(thr) + ".xdmf";
		status = myxdmf.create_voronoicell_debug_file( thisxmlfn, NumberOfFacets,
				thrioinfo.topo_n, thrioinfo.geom_n, NumberOfFacets, thrh5fn );

	} //write next thread's local results

	double mytoc = omp_get_wtime();
	cout << "Reporting of tessellation for individual threads via XDMF/HDF5 file took " << (mytoc-mytic) << " seconds" << endl;
	return 0;
}
