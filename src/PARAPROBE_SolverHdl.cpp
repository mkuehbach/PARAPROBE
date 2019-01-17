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


#include "PARAPROBE_SolverHdl.h"


solver::solver()
{
	healthy = true;
	owner = NULL;

	recon = NULL;
	try {
		recon = new reconstructor;
		reporting( "Solver constructor allocated reconstructor");
	}
	catch (bad_alloc &croak) {
		stopping( "Solver constructor unable to allocate reconstructor!");
		healthy = false;
	}

	sp = NULL;
	try {
		sp = new decompositor;
		reporting( "Solver constructor allocated decompositor");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate decompositor!");
		healthy = false;
	}

	binner = NULL;
	try {
		binner = new vxlizer;
		reporting( "Solver constructor allocated vxlizer");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate vxlizer!");
		healthy = false;
	}

	surf = NULL;
	try {
		surf = new surfacer;
		reporting( "Solver constructor allocated surfacer");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate surfacer!");
		healthy = false;
	}

	rndmizer = NULL;
	try {
		rndmizer = new rndlabeler;
		reporting( "Solver constructor allocated rndlabeler");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate rndlabeler!");
		healthy = false;
	}

	hometrics = NULL;
	try {
		hometrics = new horderdist;
		reporting( "Solver constructor allocated spatstat");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate spatstat!");
		healthy = false;
	}

	cldetect = NULL;
	try {
		cldetect = new clusterer;
		reporting( "Solver constructor allocated clusterer");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate clusterer!");
		healthy = false;
	}

	aptcrystallo = NULL;
	try {
		aptcrystallo = new aptcrystHdl;
		reporting( "Solver constructor allocated aptcrystHdl");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate aptcrystHdl!");
		healthy = false;
	}

	tessellator = NULL;
	try {
		tessellator = new tessHdl;
		reporting( "Solver constructor allocated Voronoi tessellator");
	}
	catch (bad_alloc &croak) {
		stopping("Solver constructor unable to allocate Voronoi tessellator!");
		healthy = false;
	}

	//bvh = NULL;
}


solver::~solver()
{
	//MK::do not clear owner as it is only a backreference to my owner solverHdl

	//##MK::delete can be applied on NULL
	if ( recon != NULL ) {
		delete recon;
		recon = NULL;
	}
	if ( sp != NULL ) {
		delete sp;
		sp = NULL;
	}
	if ( binner != NULL ) {
		delete binner;
		binner = NULL;
	}
	if ( surf != NULL ) {
		delete surf;
		surf = NULL;
	}
	if ( rndmizer != NULL ) {
		delete rndmizer;
		rndmizer = NULL;
	}
	if ( hometrics != NULL ) {
		delete hometrics;
		hometrics = NULL;
	}
	if ( cldetect != NULL ) {
		delete cldetect;
		cldetect = NULL;
	}
	if ( aptcrystallo != NULL ) {
		delete aptcrystallo;
		aptcrystallo = NULL;
	}
	if ( tessellator != NULL ) {
		delete tessellator;
		tessellator = NULL;
	}

	//if ( bvh != NULL ) {
	//	delete bvh;
	//	bvh = NULL;
	//}
}


bool solver::i_take_care( const size_t id, const unsigned long threadid, const unsigned int totalthreads )
{
	if ( (1+id) % totalthreads != threadid ) //other guy takes care
		return false;
	return true;
}


void solver::volume_reconstruction()
{
	if ( healthy != true ) {
		stopping("Reconstruction::Solver is not healthy");
		healthy = false; return;
	}

	if ( Settings::SyntheticTips == true
		|| Settings::ReconstructionAlgo == E_RECON_READFROMSYNTH
		|| Settings::ReconstructionAlgo == E_RECON_READFROMPOS
		|| Settings::ReconstructionAlgo == E_RECON_READFROMEPOS ) {
		reporting( "Accepting reconstructed x,y,z rawdata as reconstruction coordinates");

		healthy = recon->reconstruction_accept_sequential();
	}
	else if ( Settings::ReconstructionAlgo == E_RECON_BAS_ETAL ) {

		healthy = recon->reconstruction_default_sequential();
	}
	else {
		stopping("Implementation reconstruction algorithm lacking");
		healthy = false;
	}
	//MK::now we have pp3 with reconstructed coordinates we can do better and organize these points spatially
	//in such a manner as to maximize parallel thread throughput by partitioning the ion locations to individidual threadlocal memory
	//but the effectiveness of such data organization will dependent on research question: for instance
	//higher order of only particular ion speci will be most efficient if ions are not only spatially organized but also by type
}


void solver::spatial_decomposition()
{
	//MK::partitioning the points on threadlocal memory to build thread local KDTrees for fast distance to surface computations, spatial range queries, and other pp3 intense operations
	sp->tip_aabb_get_extrema();

	sp->loadpartitioning();
}


void solver::volume_binning()
{
	binner->rectangular_binning();
}


void solver::surface_triangulation()
{
	//MK::why tip surface identification and triangulation ?
	//ions may lay too close to tip surface to allow unbiased computation of higher order neighbors
	//it is interesting to approximate the local tip surface curvature to enable computing a surface proximity based on which
	//the dataset can be segmented into ions living close to the surface and ions deeply embedded in the tip volume
	//ions deeply embedded have enough ions close to them such that unbiased statistics can be identified
	//for ions close to the surface spatial range query balls may protrude into chamber volume in which there is no ion,
	//thereby biasing the computation of spatial statistics
	//key idea we use the CGAL library for surface extraction has alpha shapes, hence also convex hulls is mature
	//we contribute advanced pruning scheme to reducing alpha shape construction costs

#ifdef UTILIZE_CGAL

	if ( Settings::SurfaceMeshingAlgo == E_NOSURFACE ) {
		reporting("Surface was not reconstructed higher order stuff may be biased as no distance information is available!");
	}
	if ( Settings::SurfaceMeshingAlgo == E_ALPHASHAPE_CGAL ) {
		if ( surf->alphashape() == true ) {
			string mess = "Surface was triangularized via AlphaShapes in " + to_string(surf->tipsurface.size()) + " triangles";
			reporting( owner->get_rank(), mess);
		}
	}
	else if ( Settings::SurfaceMeshingAlgo == E_CONVEXHULL_CGAL ) {
		if ( surf->convexhull() == false ) {
			stopping( "Convex hull currently not implemented!");
			//does not necessarily prevent all bias as ions can sit close to locally concave surface patches");
		}
	}
	else if ( Settings::SurfaceMeshingAlgo == E_MARCHINGCUBE_IGL ) {
		if ( surf->marchingcubes() == false ) {
			stopping( "Marching cubes currently not implemented!");
		}
	}
	else if ( Settings::SurfaceMeshingAlgo == E_LOADEXISTENTVTK ) {
		if ( surf->tipsurface.size() > 0 ) {
			stopping( "A triangle surface was already reconstructed!");
		}
		else {
			if ( surf->read_vtk_io_tipsurface( Settings::SurfaceFilenameIn ) == false ) {
				stopping( "Unable to load existent triangulation from file!");
			}
		}
	}
	else {
		stopping( "Unknown triangulation algorithm!");
	}


	if ( Settings::IOTriangulation == true && surf->tipsurface.size() > 0 ) {
		surf->report_tipsurface();
	}

#else

	if ( Settings::SurfaceMeshingAlgo == E_NOSURFACE ) {
		reporting("Surface was not reconstructed higher order stuff may be biased as no distance information is available!");
	}
	else if ( Settings::SurfaceMeshingAlgo == E_LOADEXISTENTVTK ) {
		if ( surf->tipsurface.size() > 0 ) {
			stopping( "A triangle surface was already reconstructed!");
		}
		else {
			if ( surf->read_vtk_io_tipsurface( Settings::SurfaceFilenameIn ) == false ) {
				stopping( "Unable to load existent triangulation from file!");
			}
		}
	}
	else {
		stopping( "Unknown triangulation algorithm!");
	}

#endif
}


void solver::surface_distancing2()
{
	//##MK::development version patching performance issues of surface_distancing1 load imbalance
	//because points are balanced across threads but triangles not, hence as tip narrows more triangle tests for same number of ions
	//##MK::key idea of this patch, global Rtree build, threadlocal filling with maximum probed distance
	//then sequentialized processing of regions coorperative distancing with all threads
	//threads collect results locally than update distances for close to the tip ions in critical region and proceed with next region

	//build an Rtree of triangles sequentially ##MK, then use in parallel for computing ion to closest triangle distance
	double tic = MPI_Wtime();

	//check if a triangle hull exists at all
	bool hullexists = ( surf->tipsurface.size() > 0 ) ? true : false;

	cout << "HullStatus ___" << hullexists << "___ tipsurface.size() " << surf->tipsurface.size() << endl;

	//##MK::DEBUG even when it does do  not compute distances e.g. for speeding up debugging
	if ( Settings::DebugComputeDistance == false ) {
		hullexists = false;
	}

	//anyway fill with default distances
	#pragma omp parallel shared(hullexists)
	{
		//unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		apt_xyz R = max(Settings::SpatStatRadiusMax, Settings::SurfaceCellsCarvingRadius);
		apt_xyz RSQR = SQR(R); //##MK::max of spat and clustering
		sp->db.at(mt)->ion2surfdistance_init( RSQR );

		#pragma omp critical
		{
			cout << "Thread " << mt << " initialized distance field with R " << R << " for " << sp->db.at(mt)->ionpp3.size() << " ions" << endl;
		}
	}

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->tictoc.prof_elpsdtime_and_mem( "Ion2SurfDistInit", APT_GEO, APT_IS_PAR, mm, tic, toc);

	cout << "TipSurface distance fields initialized in " << (toc-tic) << " seconds" << endl;

	if ( hullexists == false ) { //done and out
		return;
	}

	tic = MPI_Wtime();

	//computing distances desired and therefore building threadglobal Rtree
	surf->build_rtree();

	toc = MPI_Wtime();

	mm = owner->tictoc.get_memoryconsumption();
	owner->tictoc.prof_elpsdtime_and_mem( "SurfTriangleHullRTreeConstruct", APT_BVH, APT_IS_SEQ, mm, tic, toc);
	cout << "Building triangle tree BVH " << (toc-tic) << " seconds" << endl;
	cout << "Now executing surface distancing ..." << endl;

	//cooperative computation of distances region by region
	tic = MPI_Wtime();

	#pragma omp parallel
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		//threadlocal buffering of results to reduce time spent in pragma omp criticals
		vector<pdist> myres;

		//threadlocal copy of link to global BVH
		//##MK::potentially build problem specific smaller trees to profit from indirect memory locality improvement
		Tree* mylinkedbvh = surf->bvh;

		//distance values, SQR to avoid sqrt computations, threadmath.closestPointOnTriangle returns squared distance!
		const apt_xyz R = max(Settings::SpatStatRadiusMax,Settings::SurfaceCellsCarvingRadius); //##MK::change at some point to max(spatstat and clustering)
		const apt_xyz RSQR = SQR(R);

		//threadlocal math
		mathHdl genius;
		vector<tri3d> const & hull = surf->tipsurface;

		//profiling
		double mytic, mytoc;

		for(unsigned int thr = MASTER; thr < nt; thr++) { //everybody executes enforced one after another processing of data for thread thr by
			vector<p3dm1> & theseions = sp->db.at(thr)->ionpp3;
			size_t ni = theseions.size();
			vector<unsigned int> candidates;
			size_t myworkload = 0;
			size_t nochecks = 0;

			//parallel for takes active thread team to distribute loop iterations in an attempt to improve load balancing
			//caching threadlocal results to myres if intruding triangles were found
			myres.clear();
			mytic = omp_get_wtime();

			//MK::distancing2, use defaults that is schedule(static,1) --> caused substantial load imbalance
			#pragma omp for schedule(dynamic,1) nowait
			for(size_t i = 0; i < ni; ++i) {
				p3d me = p3d( theseions.at(i).x, theseions.at(i).y, theseions.at(i).z ); //reffering to memory owned by potentially different thread than me

				//prune most triangles by utilizing bounded volume hierarchy (BVH) Rtree
				AABB e_aabb( trpl(me.x-R, me.y-R, me.z-R), trpl(me.x+R, me.y+R, me.z+R) );

				//if there is no triangle inside spherical ball about me, I am definately inside
				candidates.clear();
				candidates = mylinkedbvh->query( e_aabb ); //MK::assure that this is read only ##MK future improve for better locality multiple Rtrees
				myworkload += candidates.size();

				//if axis-aligned bounding box about triangle does not intersect this probe volume triangle must be farther away than RCutoff
				if ( candidates.size() > 0 ) { //no triangle close nothing to do distance value was already set, most likely case
					apt_xyz CurrentBest = RSQR;
					apt_xyz CurrentValue = RSQR;
					for ( size_t c = 0; c < candidates.size(); ++c ) {
						size_t triangleid = candidates.at(c);
						CurrentValue = genius.closestPointOnTriangle( hull.at(triangleid), me ); //##MK::far field fetch but much fewer candidates in total
						if ( CurrentValue > CurrentBest ) { //most likely
							continue;
						}
						else {
							CurrentBest = CurrentValue;
			//cout << "i/c/nc/best " << i << "\t\t" << c << "\t\t" << nc << "\t\t" << CurrentBest << endl;
						}
					}
					//MK::CurrentBest <= RSQR
					myres.push_back( pdist(i, CurrentBest) );
				} //all candidate triangles probed
				else {
					++nochecks;
				}
			} //continue to process ions in parallel across all threads

			//MK::implicit barrier at pragma omp for unless nowait clause

			#pragma omp critical
			{
				//MK::updating values in memory of thr! by other threads and thr
				size_t nj = myres.size();
				vector<apt_xyz> & there = sp->db.at(thr)->ion2surf;
				for(size_t j = 0; j < nj; ++j) {
					pdist thatone = myres.at(j);
					there.at(thatone.pid) = thatone.dist;
				}
				mytoc = omp_get_wtime();
				cout << "Thread " << mt << " processed " << myres.size() << "/" << nochecks << " with/without checks @" << myworkload << " triangles for region " << thr << " took " << (mytoc-mytic) << " seconds" << endl;
			}

			myres.clear();

			mytoc = omp_get_wtime();

			//threads work cooperatively on all ions of region
			//##MK::assuming regions are large and potentially billions of triangle distance to be computed load
			//should be much better distributed, making it at least for now admissible to spend some time at barrier
			#pragma omp barrier

			#pragma omp master
			{
				cout << "Cooperative distancing for ions on thread " << thr << " successful took " << (mytoc-mytic) << "/" << (omp_get_wtime()-mytic) << " without/with barrier" << endl;
			} //no barrier at the end of omp master
		} //distance ions of next region
	}

	//technically BVH no longer required
	//surf->chop_rtree();

	toc = MPI_Wtime();
	mm = owner->tictoc.get_memoryconsumption();
	owner->tictoc.prof_elpsdtime_and_mem( "Ion2SurfThreadparallelDistancing", APT_GEO, APT_IS_PAR, mm, tic, toc);
	cout << "Computing parallelized distance to surface in " << (toc-tic) << " seconds" << endl;
}


/*
void solver::characterize_distances()
{
	//MK::do something with distance information
	//##MK::simplest as of now report in vtk file
	size_t nevents = 0;
	for(size_t mt = MASTER; mt < sp->db.size(); mt++ ) { //deterministic output for given loadpartitioning and a priori fixed number of threads
		if ( sp->db.at(mt)->ion2surf.size() != sp->db.at(mt)->ionpp3.size() ) {
			string mess = "Thread " + to_string(mt) + " dissimilar array size ionpp3, ion2surf";
			complaining(mess);
			return;
		}
		nevents += sp->db.at(mt)->ion2surf.size();
	}

	if ( Settings::IOIonTipSurfDists == true ) {
		double tic = MPI_Wtime();

		//MK::write VTK file showing the positions of all ions in reconstructed space, type and distance to surface
		string vtk_io_fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".IonSurfDistance.vtk";
		string mess = "VTKIO writing ion distance to surface to " + vtk_io_fn;
		reporting( mess );

		ofstream vtk;
		vtk.open( vtk_io_fn.c_str() );
		if ( vtk.is_open() == true ) {
			vtk << "# vtk DataFile Version 2.0\n"; //header
			vtk << "PARAPROBE AlphaShapeTriHullDistance\n";
			vtk << "ASCII\n";
			vtk << "DATASET POLYDATA\n";
			vtk << "POINTS " << nevents << " double\n";

			//sequential outputting thread by thread
			for(size_t mt = MASTER; mt < sp->db.size(); mt++ ) { //deterministic output for given loadpartitioning and a priori fixed number of threads
				vector<p3dm1>& these = sp->db.at(mt)->ionpp3;
				size_t ni = these.size();
				for(size_t i = 0; i < ni; ++i) {
					vtk << these.at(i).x << " " << these.at(i).y << " " << these.at(i).z << "\n";
				}
			}
			vtk << endl;
			vtk << "VERTICES " << nevents << " " << 2*nevents << "\n";
			for ( size_t e = 0; e < nevents; ++e ) {
				vtk << 1 << " " << e << "\n";
			}

			//MK::ranged ion type as field data for coloring in Paraview
			vtk << "POINT_DATA " << nevents << "\n";
			vtk << "FIELD FieldData 2\n";
			vtk << "IonType 1 " << nevents << " float\n"; //do not make int because numeric_limits<unsigned int>::max() is not readable then
			for(size_t mt = MASTER; mt < sp->db.size(); mt++ ) { //deterministic output for given loadpartitioning and a priori fixed number of threads
				vector<p3dm1>& these = sp->db.at(mt)->ionpp3;
				size_t ni = these.size();
				for(size_t i = 0; i < ni; ++i) {
					vtk << these.at(i).m << "\n";
				}
			}
			vtk << "Dist2Surf 1 " << nevents << " float\n";
			for(size_t mt = MASTER; mt < sp->db.size(); mt++ ) { //deterministic output for given loadpartitioning and a priori fixed number of threads
				vector<apt_xyz>& these = sp->db.at(mt)->ion2surf;
				size_t ni = these.size();
				for(size_t i = 0; i < ni; ++i) {
					vtk << these.at(i) << "\n";
				}
			}
			vtk << endl;

			vtk.flush();
			vtk.close();
		}

		double toc = MPI_Wtime();
		memsnapshot mm = memsnapshot();
		owner->tictoc.prof_elpsdtime_and_mem( "VTKFileOutputIon2Surf", APT_IO, APT_IS_SEQ, mm, tic, toc);
	}

	cout << "Ion distance to surface distribution characterized " << endl;
}
*/


void solver::characterize_tip()
{
	if ( Settings::IOReconstruction == true || Settings::IOIonTipSurfDists == true ) {
		double tic = MPI_Wtime();

		size_t NumberOfIons = 0;
		size_t nthr = sp->db.size();
		vector<io_bounds> iohelp( nthr, io_bounds() );
		vector<io_bounds> topohelp( nthr, io_bounds() );
		for( size_t thr = MASTER; thr < nthr; thr++ ) {
			if ( sp->db.at(thr)->ionpp3.size() == sp->db.at(thr)->ion2surf.size() ) {
				NumberOfIons += sp->db.at(thr)->ionpp3.size();
				iohelp.at(thr).n = sp->db.at(thr)->ionpp3.size();
				iohelp.at(thr).s = 0;
				topohelp.at(thr).n = sp->db.at(thr)->ionpp3.size() * 3; //MK::triplet of Polyvertex keyword, number of vertices, and global vertex ID
				topohelp.at(thr).s = 0;
				for ( size_t trailing = 0; trailing < thr; ++trailing ) {
					iohelp.at(thr).s += iohelp.at(trailing).n;
					topohelp.at(thr).s += topohelp.at(trailing).n;
				}
				iohelp.at(thr).e = iohelp.at(thr).s + iohelp.at(thr).n;
				topohelp.at(thr).e = topohelp.at(thr).s + topohelp.at(thr).n;
			}
			else {
				complaining( "Container length for ionpp3 and ion2surf inconsistent" ); return;
			}
		}

		if ( NumberOfIons >= UINT32MX ) {
			complaining( "This I/O function cannot handle more than UINT32MX ions change uint32 topology type to size_t" );
			return;
		}

		int status = 0;
		h5iometa ifo = h5iometa();
		h5offsets offs = h5offsets();

		//##MK::add topology pieces of information allows to filter points in Paraview and VisIt
		ifo = h5iometa( PARAPROBE_VOLRECON_TOPO, 3*NumberOfIons, 1 );
		status = owner->resultsh5Hdl.create_contiguous_matrix_u32le( ifo );
		unsigned int GlobalVertexID = 0; //##MK::change to size_t if more than UINT32MX ions

		ifo = h5iometa( PARAPROBE_VOLRECON_XYZ, NumberOfIons, 3 );
		status = owner->resultsh5Hdl.create_contiguous_matrix_f32le( ifo );

		ifo = h5iometa( PARAPROBE_VOLRECON_IONTYPE_IDS, NumberOfIons, 1 );
		status = owner->resultsh5Hdl.create_contiguous_matrix_u8le( ifo );

		if ( Settings::IOIonTipSurfDists == true ) {
			ifo = h5iometa( PARAPROBE_VOLRECON_SURFDISTSQR, NumberOfIons, 1 );
			status = owner->resultsh5Hdl.create_contiguous_matrix_f32le( ifo );
			reporting( "Ion distance I/O corresponding HDF5 groups and datasets initialized" );
		}

		for( size_t thr = MASTER; thr < nthr; thr++ ) {
			vector<float> f32thrbuf;
			vector<unsigned char> u8thrbuf;
			f32thrbuf.reserve( 3*iohelp.at(thr).n );
			u8thrbuf.reserve( 1*iohelp.at(thr).n );
			for( auto it = sp->db.at(thr)->ionpp3.begin(); it != sp->db.at(thr)->ionpp3.end(); ++it ) {
				f32thrbuf.push_back( it->x );
				f32thrbuf.push_back( it->y );
				f32thrbuf.push_back( it->z );
				unsigned char itypid = (it->m >= static_cast<unsigned int>(UCHARMX)) ? UCHARMX : static_cast<unsigned char>(it->m);
				//##MK::replace by unsigned short if more than 255 types required
				u8thrbuf.push_back( itypid );
			}

			ifo = h5iometa( PARAPROBE_VOLRECON_XYZ, NumberOfIons, 3 );
			offs = h5offsets( iohelp.at(thr).s, iohelp.at(thr).e, 0, 3, f32thrbuf.size() / 3, 3);
			cout << "Thread " << thr << " xyz reconstructed coordinates" << endl;
			status = owner->resultsh5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, f32thrbuf );
			cout << status << endl;
			f32thrbuf = vector<float>();

			ifo = h5iometa( PARAPROBE_VOLRECON_IONTYPE_IDS, NumberOfIons, 1 );
			offs = h5offsets( iohelp.at(thr).s, iohelp.at(thr).e, 0, 1, u8thrbuf.size(), 1);
			cout << "Thread " << thr << " iontype reconstructed ions" << endl;
			status = owner->resultsh5Hdl.write_contiguous_matrix_u8le_hyperslab( ifo, offs, u8thrbuf );
			cout << status << endl;
			u8thrbuf = vector<unsigned char>();

			//##MK::we add XDMF Polyvertex topology support
			//therewith we can filter ions based on attribute data within Paraview
			vector<unsigned int> u32topobuf( 3*iohelp.at(thr).n, 1 ); //MK::we need to store a triplet of XDMF topology type count and vertex ID
			size_t ni = sp->db.at(thr)->ionpp3.size();
			for ( size_t i = 0; i < ni; ++i ) { //two loop instead of one to avoid having to allocate even larger intermediate buffer
				u32topobuf.at((3*i)+2) = GlobalVertexID;
				GlobalVertexID++; //consistent numbering of vertices across threads
			}

			ifo = h5iometa( PARAPROBE_VOLRECON_TOPO, 3*NumberOfIons, 1 );
			offs = h5offsets( topohelp.at(thr).s, topohelp.at(thr).e, 0, 1, u32topobuf.size(), 1);
			cout << "Thread " << thr << " topo information" << endl;
			status = owner->resultsh5Hdl.write_contiguous_matrix_u32le_hyperslab( ifo, offs, u32topobuf );
			cout << status << endl;
			u32topobuf = vector<unsigned int>();

			if ( Settings::IOIonTipSurfDists == true ) {
				ifo = h5iometa( PARAPROBE_VOLRECON_SURFDISTSQR, NumberOfIons, 1);
				offs = h5offsets( iohelp.at(thr).s, iohelp.at(thr).e, 0, 1, iohelp.at(thr).n, 1);
				cout << "Thread " << thr << " ion2surf distance squared" << endl;
				status = owner->resultsh5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, sp->db.at(thr)->ion2surf );
				cout << status << endl;
				//MK::do not delete ion2surf values, will need them still!
			}
		}

		string xdmffn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VolReconstruction.xdmf";
		//owner->xdmfmetaHdl.create_iondistance_file( xdmffn, NumberOfIons, owner->resultsh5Hdl.h5resultsfn );
		owner->xdmfmetaHdl.create_volrecon_file( xdmffn, NumberOfIons, owner->resultsh5Hdl.h5resultsfn );

		double toc = MPI_Wtime();
		memsnapshot mm = memsnapshot();
		owner->tictoc.prof_elpsdtime_and_mem( "IOHDF5Ion2SurfDistance", APT_IO, APT_IS_SEQ, mm, tic, toc);
		cout << "Ion distance to surface distribution characterized " << endl;
	}
}


void solver::init_spatialindexing()
{
	double tic = MPI_Wtime();

	bool AllKDTreesValid = true; //attempt to disprove

	#pragma omp parallel shared(AllKDTreesValid)
	{
		//unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		bool mystatus = sp->db.at(mt)->build_local_kdtree();

		#pragma omp critical
		{
			if ( mystatus == false )
				AllKDTreesValid = false;
		}
	} //end of parallel region explicit barrier

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->tictoc.prof_elpsdtime_and_mem( "ThreadparallelKDTreeConstruct", APT_BVH, APT_IS_PAR, mm, tic, toc);

	tic = MPI_Wtime();

	if ( AllKDTreesValid == true ) {
		sp->kdtree_success = true;
		reporting( owner->get_rank(), "KDTree construction was globally successful!");

		if ( Settings::IOKDTreePartitioning == true ) {
			sp->reportpartitioning();
		}
	}
	else {
		sp->kdtree_success = false;
		stopping( owner->get_rank(), "KDTree construction was not successful globally!");
	}

	toc = MPI_Wtime();
	mm = owner->tictoc.get_memoryconsumption();
	owner->tictoc.prof_elpsdtime_and_mem( "ThreadparallelKDTreeReportResults", APT_IO, APT_IS_SEQ, mm, tic, toc);
	cout << "Threadlocal KDTrees constructed in parallel in " << (toc-tic) << " seconds" << endl;
}



void solver::characterize_crystallography()
{
	cout << "Starting with Vicente Araullo-Peters crystallographic feature identification..." << endl;
	aptcrystallo->compute_crystallography_vicsmethod();

	cout << "Reporting results of crystallographic feature extraction..." << endl;
	aptcrystallo->report_crystallography_results2();
	aptcrystallo->delete_crystallography_results();
}


void solver::characterize_spatstat()
{
	//MK::provided that all trees are true
	if ( sp->kdtree_success == true ) {
		cout << "Starting with spatial distribution function characterization mode " << Settings::SpatialDistributionTask << endl;
		cout << "Descriptive spatial statistics mode code ___" << Settings::DescrStatTasksCode << "___" << endl;

		//identify which specific iontype pairing to analyze
		hometrics->initialize_descrstat_tasks();

		//perform analysis for all tasks using initial assignment of the labels
		//hometrics->compute_generic_spatstat2(); //stable used for Paper14
		hometrics->compute_generic_spatstat3(); //developmental with ordering

		if ( Settings::SpatStatAddLabelRnd == true ) {
			//MK::now randomize iontype assignment over entire tip to get randomize spatial statistics
			rndmizer->reset();
			rndmizer->learn_original_types_global();
			rndmizer->shuffle_types_mt19937();
			rndmizer->apply_shuffled_types_global();

			//so characterize again but this time with randomization
			//hometrics->compute_generic_spatstat2();
			hometrics->compute_generic_spatstat3();

			rndmizer->reset_original_types_global(); //MK::necessary otherwise subsequent tasks use incorrect physical ion labeling!
		}
	}
	else {
		stopping( owner->get_rank(), "KDTree construction was previously not successful globally!");
	}
}


void solver::characterize_clustering()
{
	if ( sp->kdtree_success == true ) {
		cout << "Starting with clustering characterization mode " << Settings::ClusteringTask << endl;
		cout << "Clustering characterization mode code ___" << Settings::ClusteringTasksCode << "___" << endl;

		cldetect->initialize_clustering_tasks();

		switch(Settings::ClusteringTask)
		{
			case E_DBSCAN:
				return; //NOT IMPLEMENTED
			case E_MAXSEPARATION:
				cldetect->maximum_separation_method();
				break;
			case E_ISOSURF:
				return; //NOT IMPLEMENTED

			default:
				return;
		}
	}
	else {
		stopping( owner->get_rank(), "KDTree construction was previously not successful globally!");
	}
}


void solver::tessellate_tipvolume()
{
	double tic = MPI_Wtime();
	memsnapshot mm1 = memsnapshot();

	tessellator->configure( Settings::IOVoronoiDescrStats,
			Settings::IOVoronoiCellPositions, Settings::IOVoronoiTopoGeom ); //set lean and heavy storage in tessHdl

	//MK::generate one such tessHdl object per running MPI process if want to split hybrid manner tessellation construction
	#pragma omp parallel shared(mm1)//input parameter shared by default, storetessellation read only so no problem if shared as well
	{
		double mytic = omp_get_wtime();
		int mt = omp_get_thread_num();
		int nt = omp_get_num_threads();

		#pragma omp master
		{
			//master allocates space for pointer to individual threadlocal workers BUT not FOR THEIR HEAVY DATA!
			for( int thr = 0; thr < nt; ++thr )
				tessellator->local_tessellations.push_back(NULL);

			if ( tessellator->we_store_tess_metainfo == true ||
					tessellator->we_store_tess_cellpos == true ||
						tessellator->we_store_tess_topogeom == true ) {
				for( int thr = 0; thr < nt; ++thr )
					tessellator->iohelp.push_back( voro_io_bounds() ); //we dont know the bounds yet...
			}
		}
		#pragma omp barrier //necessary omp master has no barrier and memory needs to be filled before correct address written to it

		//build thread-local instances of tessellation object
		//the thread-local workers allocate their own memory
		//given that we are in parallel region they will pull memory preferentially and more frequently
		//IF AND ONLY IF sufficient local memory still available and OS instructed to have thread pulling from
		//local memory (this is why we pin the thread by overwriting the thread affinity mask and why it is
		//useful in addition to overwrite the standard STL malloc3 allocator by e.g. jemalloc
		//IMPORTANT IS further that we not only allocate the memory (as this will only allocate virtually)
		//but we have to write to it, this is known as first touch police, so allocate locally and write to it
		tess* ltess = NULL;
		ltess = new tess; //##MK::exception handling try/catch
		ltess->i_store_tess_metainfo = tessellator->we_store_tess_metainfo; //first touch
		ltess->i_store_tess_cellpos = tessellator->we_store_tess_cellpos;
		ltess->i_store_tess_topogeom = tessellator->we_store_tess_topogeom;
		ltess->owner = tessellator;
		ltess->surfacehull = surf->bvh;
		ltess->trianglehull = &surf->tipsurface;
		ltess->surfdistance_ids = &sp->db.at(mt)->ionpp3_ionpp3idx;
		ltess->surfdistance_vals = &sp->db.at(mt)->ion2surf;
		//ltess->ion2pp3surf = &sp->db.at(mt)->ion2surf;

		//pull atoms into local objects based on bounds, as this writes to the memory of the tess object it does the first-touch
		ltess->set_zbounds( sp->spatialsplits.at(mt) );
		ltess->set_halosize_p3d64( sp->halothickness );
		ltess->set_haloinfo();

		ltess->pull_ionportion_inplace( sp->db.at(mt)->ionpp3_tess, sp->tip );

		//##MK::clear copy in databasesp->db.at(mt)->ionpp3_tess.clear();
		sp->db.at(mt)->ionpp3_tess = vector<p3dm3>();
		//MK::do not clear ionpp3_ionpp3idx, we need this during the build phase to identify
		//the distance of a point in ionportion to the triangulated surface to potentially avoid
		//counting Voronoi cell to the surface of the dataset whose topology and geometry in affect by boundary
		//! MK::this is necessary currently because we do not have a exact triangle manifold and manifold cell intersection algorithm

		ltess->build_tessportion0(); //MK::does not make any correction attempt for cells at boundary just detects them based on distance and does not report them subsequently
		//! MK::this allows to erode the surface Voronoi cells with incomplete or faulty geometry from the dataset
		//! MK::without recomputing the exact distance of a point to the surface regardless whether hull is a manifold or not

		//ltess->build_tessportion1(); //MK::does not get cells at the surface correct, ackward spaceship-formed-elongated cells even for deeper embedded ions how numerically never make wall contact
		//ltess->build_tessportion2(); //MK::not tested but not expected to work because for cases of concave/convex local surface patch configurations fit plane cannot SO FAR simplistically assured to be always in front or on the ion position hence cut procedure by bisecting with fake voronoi cell generator point to use voroxx plane function does not work consistently
		//ltess->build_tessportion3(); //MK::excludes cell if respective ion is closer to the boundary than R <= Settings::SpatialStatsRadiusMax
		//##MK::pull a memory snapshot of program here because now heavy data still exist in memory!

		//now we do not need any longer the ionpp3_ionpp3idx helper array
		sp->db.at(mt)->ionpp3_ionpp3idx = vector<unsigned int>();

		#pragma omp master
		{
			mm1 = owner->tictoc.get_memoryconsumption();
		}
		//necessary all data have to be processed before cooperatively instructing further tasks
		double mytoc = omp_get_wtime();

		#pragma omp critical
		{
			tessellator->local_tessellations.at(mt) = ltess; //linking address of thread-local object data to tessHdl necessary at the latest here before cooperative stuff

			if ( tessellator->we_store_tess_metainfo == true ) {
				//communicate storage demands
				//if ( ltess->myprecisiondemand.prec > allprecisiondemand.prec ) { allprecisiondemand.prec = ltess->myprecisiondemand.prec; }
				tessellator->allprecisiondemand.n_cells += ltess->myprecisiondemand.n_cells;
				tessellator->iohelp.at(mt).cell_n = ltess->myprecisiondemand.n_cells;

				if (  tessellator->we_store_tess_topogeom == true ) {
					tessellator->allprecisiondemand.n_topo += ltess->myprecisiondemand.n_topo;
					tessellator->allprecisiondemand.n_geom += ltess->myprecisiondemand.n_geom;
					//communicate also io implicit array bound demands for writing data
					tessellator->iohelp.at(mt).topo_n = ltess->myprecisiondemand.n_topo;
					tessellator->iohelp.at(mt).geom_n = ltess->myprecisiondemand.n_geom;
					//MK::latter values will be translated into absolute position intervals on implicit 1d data arrays during I/O operation
				}
			}
			cout << "Thread " << mt << " building local tessellation took " << (mytoc-mytic) << " seconds" << endl;
		}
	} //end of parallel region

	double toc = MPI_Wtime();
	owner->tictoc.prof_elpsdtime_and_mem( "ComputeTessellation", APT_TES, APT_IS_PAR, mm1, tic, toc);
	cout << "Voro++ powered threaded Voronoi tessellation computation completed took " << (toc-tic) << " seconds" << endl;

	tic = MPI_Wtime();

	//##MK::consider in the future to do multi-threaded I/O still
	//tessellator->report_tessellation_contiguous_atonce(); //##MK::in the future change
#ifdef VALIDZONE_IONS_ONLY
	tessellator->report_tessellation_hyperslab_perthread();
#else
	//##MK::use the following function only in combination with undefined VALIDZONE_IONS_ONLY
	//to report per thread based results with halo regions
	tessellator->report_tessellation_hyperslab_onlythethread();
#endif

	toc = MPI_Wtime();
	memsnapshot mm2 = owner->tictoc.get_memoryconsumption();
	owner->tictoc.prof_elpsdtime_and_mem( "HDF5XDMFReportVoronoiTess", APT_IO, APT_IS_SEQ, mm2, tic, toc);
	cout << "HDF5/XDMF based reporting of tessellation completed took " << (toc-tic) << " seconds" << endl;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

reconstructor::reconstructor()
{
	owner = NULL;
	healthy = true;
}


reconstructor::~reconstructor()
{
	//MK::do not delete owner, it is only a backreference to the solver owing the reconstructor

	size_t nb = pp3.size();
	for (size_t b = 0; b < nb; ++b) {
		delete pp3.at(b); pp3.at(b) = NULL;
	}

	nb = lbls.size();
	for (size_t b = 0; b < nb; ++b) {
		delete lbls.at(b); lbls.at(b) = NULL;
	}
}



bool reconstructor::reconstruction_accept_sequential( void )
{
	double tic = MPI_Wtime();

	//feed through rawdata_pos/_epos, where are they?
	solverHdl* here = owner->owner;
	size_t nb = 0;
	if ( Settings::SyntheticTips == true )
		nb = here->rawdata_alf.buckets.size();
	else {
		if ( Settings::InputFileformat == E_IO_EPOS )
			nb = here->rawdata_epos.size();
		else if ( Settings::InputFileformat == E_IO_POS || Settings::InputFileformat == E_IO_HDF5 )
			nb = here->rawdata_pos.size();
		else //nothing
			return false;
	}

	//MK::if process pragma omp parallel here potential to comply with first touch policy that OpenMP threads get portions of the rawdata in threadlocal memory
	//MK::HOWEVER, at this stage the rawdata in the solverHdl object are not necessarily reconstructed yet, therefore we have only a very approximate
	//MK>>idea of a good spatial partitioning of these, in particular if we seek to reconstruct not all ions, therefore, we better build the OpenMP
	//thread/local stuff after the reconstruction stage...
	for ( size_t b = 0; b < nb; ++b) {
		size_t n = 0;
		if ( Settings::SyntheticTips == true ) {
			if ( here->rawdata_alf.buckets.at(b) != NULL )
				n = here->rawdata_alf.buckets.at(b)->size();
		}
		else {
			if ( Settings::InputFileformat == E_IO_EPOS ) {
				if ( here->rawdata_epos.at(b) != NULL )
					n = here->rawdata_epos.at(b)->size();
			}
			else { //if ( Settings::InputFileformat == E_IO_POS || Settings::InputFileformat == E_IO_HDF5 ) {
				if ( here->rawdata_pos.at(b) != NULL )
					n = here->rawdata_pos.at(b)->size();
			}
		}

//cout << "Reconstruction accept bucket b/n number of buckets " << b << ";" << n << endl;

		//a bucket can at least be empty
		pp3.push_back(NULL);
		lbls.push_back(NULL);

		if ( n > 0) { //if data exist...
			vector<p3d>* wpbucket = NULL;
			vector<unsigned int>* wlbucket = NULL;
			try {
				wpbucket = new vector<p3d>;
				wlbucket = new vector<unsigned int>;
			}
			catch (bad_alloc &reconexc) {
				complaining( "Unable to allocate memory in reconstructor for storing reconstructed ion locations");
				return false; //already allocated stuff class destructor will take care of
			}

			if ( Settings::SyntheticTips == true ) {
				vector<pos>* rp = here->rawdata_alf.buckets.at(b);
				vector<unsigned int>* rl = here->rawdata_iontype.at(b);
				for(size_t i = 0; i < n; ++i) {
					wpbucket->push_back( p3d(rp->at(i).x, rp->at(i).y, rp->at(i).z) ); //##use fill routine
					wlbucket->push_back( rl->at(i) );
				}
			}
			else {
				if ( Settings::InputFileformat == E_IO_EPOS ) {
					vector<epos>* rp = here->rawdata_epos.at(b);
					vector<unsigned int>* rl = here->rawdata_iontype.at(b);
					for(size_t i = 0; i < n; ++i) {
						wpbucket->push_back( p3d(rp->at(i).x, rp->at(i).y, rp->at(i).z) );
						wlbucket->push_back( rl->at(i) );
					}
				}
				else { //Settings::InputFileformat == E_IO_POS
					vector<pos>* rp = here->rawdata_pos.at(b);
					vector<unsigned int>* rl = here->rawdata_iontype.at(b);
					for(size_t i = 0; i < n; ++i) {
						wpbucket->push_back( p3d(rp->at(i).x, rp->at(i).y, rp->at(i).z) );
						wlbucket->push_back( rl->at(i) );
					}
				}
			}
			//register atom in pp3 and lbls
			//it gets filled
			pp3.at(b) = wpbucket;
			lbls.at(b) = wlbucket;

//cout << "Added databucket pp3.size()/lbls.size() " << pp3.back()->size() << ";" << lbls.back()->size() << endl;
		}
	} //process next bin

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "AcceptReconstruction", APT_REC, APT_IS_SEQ, mm, tic, toc);
	reporting( "Reconstruction accepted in " + to_string(toc-tic) + " seconds");

	/*if ( Settings::IOReconstruction == true ) {
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".AcceptedRecon.vtk";
		reconstruction_vtk( pp3, lbls, runparm( 0.f, 0.f, 0.f, 0, 0, 0), fn );

		toc = MPI_Wtime();
		memsnapshot mm = memsnapshot();
		owner->owner->tictoc.prof_elpsdtime_and_mem( "VTKFileOutputReconstruction", APT_IO, APT_IS_SEQ, mm, tic, toc);
	}*/

	return true;
}


bool reconstructor::reconstruction_default_sequential()
{
	//##MK::Bas et al. reconstruction algorithm

	if ( Settings::InputFileformat != E_IO_EPOS ) { //voltage data are not in *.pos file therefore not reconstruction possible
		return false;
	}

	double tic = MPI_Wtime();

	//feed from rawdata_epos, where are they?
	solverHdl* here = owner->owner;
	size_t nb = here->rawdata_epos.size();

	vector<epos>* rbucket = NULL;
	vector<p3d>* wpbucket = NULL;
	for( size_t b = 0; b < nb; ++b ) {
		rbucket = here->rawdata_epos.at(b);
		wpbucket = NULL;
		pp3.push_back(NULL);
		if ( rbucket != NULL ) {
			size_t n = rbucket->size();
			try {
				wpbucket = new vector<p3d>;
				wpbucket->reserve(n);
			}
			catch (bad_alloc &reconexc) {
				complaining( "Unable to allocate memory in reconstructor for storing reconstructed ion locations");
				return false;
			}
			for(size_t i = 0; i < n; ++i)
				wpbucket->push_back( p3d() );

			pp3.back() = wpbucket;
		}
	}

	reporting( here->get_rank(), "Bas et al. reconstruction running...");

	size_t nevents = static_cast<size_t>(here->nevt);
	apt_real ETA = (Settings::DetEffMin + Settings::DetEffMax) * 0.5;
	apt_real KF = (Settings::KFMin + Settings::KFMax) * 0.5;
	apt_real ICF = (Settings::ICFMin + Settings::ICFMax) * 0.5;

	//##MK::as of now process via one buffer only

	//allocate recon temporary memory
	//##MK::should also become leaner
	apt_real* ang = NULL;
	apt_real* thetaP = NULL;
	apt_real* theta = NULL;
	apt_real* vdcpu = NULL;
	try { ang = new apt_real[nevents]; }
	catch (bad_alloc&reconexc) {
		stopping( here->get_rank(), "Unable to allocate temporary reconstruction array ang in Bas reconstruction!"); return false;
	}
	try { thetaP = new apt_real[nevents]; }
	catch (bad_alloc&reconexc) {
		stopping( here->get_rank(), "Unable to allocate temporary reconstruction array thetaP in Bas reconstruction!"); return false;
	}
	try { theta = new apt_real[nevents]; }
	catch (bad_alloc&reconexc) {
		stopping( here->get_rank(), "Unable to allocate temporary reconstruction array theta in Barr reconstruction!"); return false;
	}
	try { vdcpu = new apt_real[nevents]; }
	catch (bad_alloc&reconexc) {
		stopping( here->get_rank(), "Unable to allocate temporary reconstruction array vdcpu in Barr reconstruction!"); return false;
	}

	reporting( here->get_rank(), "Temporary reconstruction arrays in Barr allocated!");

	//begin reconstruction
	size_t e = 0;

	//Cartesian detector coordinates into polar coordinates
	//in-place effective detector area
	//compressed angle
	//prep for radius evolution from voltage curve (in nm)
	apt_real maxrad = 0.f;
	apt_real _flen = 1.0 / Settings::FlightLength;
	apt_real _kfFevap = 1.0 / (KF * Settings::EvaporationField);

	for( size_t b = 0; b < nb; ++b ) {
		rbucket = here->rawdata_epos.at(b);
		size_t n = rbucket->size();
		for(size_t i = 0; i < n; ++i) {
			ang[e] = atan2( rbucket->at(i).dety, rbucket->at(i).detx );

			apt_real rad = sqrt( SQR(rbucket->at(i).detx) + SQR(rbucket->at(i).dety) );

			thetaP[e] = atan( rad * _flen );
			vdcpu[e] = rbucket->at(i).vdc + rbucket->at(i).vpu; //dont forget later to * _kfFevap;

//cout << ang[e] << "\t\t" << rad << "\t\t" << thetaP[e] << "\t\t" << vdcpu[e] << endl;

			e++;
			if ( rad < maxrad )
				continue;
			else
				maxrad = rad; //>=
		}
	}

	//radius evolution from voltage curve (in nm) in-place
	//launch angle
	apt_real ICF1 = ICF - 1.0;
	for (    e = 0; e < nevents; ++e ) {
		//double sinthP = ICF1 * sin(thetaP[e]);
		theta[e] = thetaP[e] + asin(ICF1 * sin(thetaP[e]));
	}

	//distance from axis and z shift of each hit


	////##MK see the equivalence
	////for(unsigned int e = 0; e < nevents; ++e) {
	////	zP[e] = rspec[e]*cos(theta[e]);
	////	d[e] = rspec[e]*sin(theta[e]);
	////	x[e] = d[e]*cos(ang[e]);
	////	y[e] = d[e]*sin(ang[e]);
	////	zP[e] = rspec[e]-zP[e]; //rspec-rspec*cos(theta)=rspec*(1-cos(theta))
	////}

	////double omegaFactor = 1.0 / Settings::AtomicDensity;
	////omegaFactor = omegaFactor * ( SQR(Settings::FlightLength) * SQR(KF) * SQR(Settings::EvaporationField) / (ETA * Adet * SQR(ICF)) );
	apt_real Adet = SQR(maxrad) * PI;
	apt_real ADENS = Settings::AtomicDensity;
	apt_real FLEN = Settings::FlightLength;
	apt_real FEVAP = Settings::EvaporationField;
	apt_real omegaFactor = 1.0 / ADENS * SQR(FLEN) * SQR(KF) * SQR(FEVAP) / (ETA * Adet * SQR(ICF));

	apt_real cumsumdz = 0.f;
	e = 0;
	for( size_t b = 0; b < nb; ++b ) {
		rbucket = here->rawdata_epos.at(b);
		if ( rbucket != NULL ) {
			size_t n = rbucket->size();
			vector<p3d>* wpbucket = pp3.at(b);

			assert( rbucket->size() == wpbucket->size() );

			for(size_t i = 0; i < n; ++i) {
			//for (   e = 0; e < nevents; ++e) {
				//zP[e] = rspec[e]*cos(theta[e]);
				//the z shift with respect to the top of the cap is Rspec - zP
				//but can be fused as it is done here because it is an elementwise operation on already known value rspec[e]
				//zP[e] = rspec[e] - rspec[e]*cos(theta[e]);

				//as zP is not summated
				//##MK::zP[e] = (vdcpu[e] * _kfFevap) * (1.0 - cos(theta[e]));

				apt_real d = (vdcpu[e] * _kfFevap) * sin(theta[e]);


				cumsumdz = cumsumdz + (omegaFactor / SQR(vdcpu[e])); //##MK::cumsumdz + dz[e]; //##MK::dz[e] = omegaFactor / SQR(vdcpu[e]);

				apt_real xrecon = d * cos(ang[e]);
				apt_real yrecon = d * sin(ang[e]);
				apt_real zrecon = cumsumdz + ((vdcpu[e] * _kfFevap) * (1.0 - cos(theta[e]))); //##MK::cumsumdz + zP[e];
			//}
				//reset dummy value with reconstruction space coordinates
				wpbucket->at(i) = p3d( xrecon, yrecon, zrecon );

//cout << setprecision(32) << wpbucket->at(i).x << "\t\t" << wpbucket->at(i).y << "\t\t" << wpbucket->at(i).z << endl;
				e++;
			}
		}
	}

	memsnapshot mm1 = owner->owner->tictoc.get_memoryconsumption();
	delete [] ang;
	delete [] thetaP;
	delete [] theta;
	delete [] vdcpu;

	//##MK::improve in future
	vector<unsigned int>* lbucket = NULL;
	vector<unsigned int>* wlbucket = NULL;
	for ( size_t b = 0; b < nb; ++b) {
		lbucket = here->rawdata_iontype.at(b);
		wlbucket = NULL;
		lbls.push_back(NULL);
		if ( lbucket != NULL ) {
			size_t n = lbucket->size();
			try {
				wlbucket = new vector<unsigned int>;
				wlbucket->reserve(n);
			}
			catch (bad_alloc &reconexc) {
				complaining( "Unable to allocate memory in reconstructor for storing reconstructed ion labels");
				return false;
				//already allocated stuff class destructor will take care of
			}
			for ( size_t i = 0; i < n; ++i)
				wlbucket->push_back( lbucket->at(i) );

			lbls.back() = wlbucket;
		}
	}


	double toc = MPI_Wtime();
	memsnapshot mm2 = owner->owner->tictoc.get_memoryconsumption();
	memsnapshot mm = (mm2.residentmem < mm1.residentmem) ? mm1 : mm2;
	owner->owner->tictoc.prof_elpsdtime_and_mem( "BasEtAlReconstruction", APT_REC, APT_IS_SEQ, mm, tic, toc);
	string mess = "Reconstruction including memory setup took " +  to_string(toc - tic) + " seconds!";
	reporting( here->get_rank(), mess);

	/*if ( Settings::IOReconstruction == true ) {
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".BasEtAlRecon.vtk";
		reconstruction_vtk( pp3, lbls, runparm( ETA, KF, ICF, 0, 0, 0), fn );

		toc = MPI_Wtime();
		owner->owner->tictoc.prof_elpsdtime_and_mem( "VTKFileOutputReconstruction", APT_REC, APT_IS_SEQ, mm, tic, toc);
	}*/

	return true;
}


threadmemory::threadmemory()
{
	owner = NULL;
	threadtree = NULL;
	tessmii = p6d64();
	tesshalowidth = p3d64();
	zmi = F32MX;
	zmx = F32MI;
	melast = false;
}


threadmemory::~threadmemory()
{
	//MK::do not delete owner only backreference to decompositor who owns me!
	if ( threadtree != NULL ) {
		delete threadtree; threadtree = NULL;
	}
}


bool threadmemory::init_localmemory( p6d64 const & mybounds, p3d64 const & myhalo, const bool mlast )
{
	//MK::CALLED FROM WITHIN PARALLEL REGION
	tessmii = mybounds;
	tesshalowidth = myhalo;
	zmi = mybounds.zmi;
	zmx = mybounds.zmx;
	melast = mlast;
	return true;
}


bool threadmemory::read_localions()
{
	//MK::CALLED FROM WITHIN PARALLEL REGION

	//picks from reconstructed volume all positions inside me into local structure
	reconstructor* pointshere = owner->owner->recon;
	size_t nbp = pointshere->pp3.size();
	reconstructor* labelshere = owner->owner->recon;
	size_t nbl = labelshere->lbls.size();

	if ( nbp != nbl ) {
		stopping("Dissimilar length of geometry and label container");
		return false;
	}

	//##MK::consider to better pack pp3 in p3dm1 to have labels directly
	size_t nb = nbp;
	//cout << "--->Read localions nb/nbp/nbl\t\t" << nb << "\t\t" << nbl << "\t\t" << nbl << endl;

	bool TessellationYesOrNo = (Settings::VolumeTessellation != E_NOTESS) ? true : false;
	//size_t IonLeakage = 0;
	unsigned int worldid = 0; //MK::consecutive ID to remain capable of storing in the measured ion evaporation sequence
	for(size_t b = 0; b < nb; ++b) {
		vector<p3d>* thesepoints = pointshere->pp3.at(b);
		vector<unsigned int>* theselabels = labelshere->lbls.at(b);
		if ( thesepoints != NULL && theselabels != NULL ) {
			size_t pn = thesepoints->size();
			size_t ln = theselabels->size();
			if ( pn == ln ) {
				//two case, if a tessellation at some point is desired we store ions and halo
				//##MK::leaner would be to store only halo additionally
				apt_real local_zmi = tessmii.zmi; //by construction the same as zmi
				apt_real local_zmx = tessmii.zmx; //by construction the same as zmx
				apt_real voro_zmi = tessmii.zmi - tesshalowidth.z;
				apt_real voro_zmx = tessmii.zmx + tesshalowidth.z;
				for( size_t i = 0; i < pn; ++i, worldid++ ) { //pruning most likely not even in halo
					p3dm1 pp = p3dm1(thesepoints->at(i).x, thesepoints->at(i).y, thesepoints->at(i).z, theselabels->at(i));
					if ( TessellationYesOrNo == true ) {
						if ( pp.z >= local_zmi && pp.z < local_zmx && pp.z >= zmi && pp.z < zmx) { //construction of this->zmi and this->zmx assures they are pairwise the same
							//in the region for a valid zone ion + the threadlocal subdomain
							ionpp3.push_back( pp );
							//the length of ionpp3_tess is potentially longer than of ionpp3 because for tessellation we need a haloregion
							//while for normal tasks we dont
							ionpp3_tess.push_back( p3dm3(pp.x, pp.y, pp.z, worldid, VALIDZONE_ION, pp.m) );
							//derived threadlocal quantities like ion2surf are ordered as in ionpp3

							//we need a helper array on ionpp3_tess which identifies whether distance information
							//is available (it is only for VALIDZONE_IONs) and an array index were we can find it on the threadlocal ion2surf
							ionpp3_ionpp3idx.push_back( ionpp3.size()-1 );
						}
						else if ( pp.z >= voro_zmi && pp.z < voro_zmx ) { //MK::important to use else if cause domain volume needs to be excluded!
							ionpp3_tess.push_back( p3dm3(pp.x, pp.y, pp.z, worldid, GUARDZONE_ION, pp.m) );
							ionpp3_ionpp3idx.push_back( DO_NOT_DEREFERENCE );
						}
						else {
							/*#pragma omp critical
							{
								cout << "Leakage1 " << setprecision(32) << pp.x << ";" << pp.y << ";" << pp.z << endl;
								cout << "localzmi/zmx " << local_zmi << ";" << local_zmx << ";" << " zmi/zmx " << zmi << ";" << zmx << endl;
							}
							IonLeakage++;*/
							continue;
						}
					}
					else {
						if ( pp.z >= zmi && pp.z < zmx) {
							//in the region for a valid zone ion + the threadlocal subdomain
							ionpp3.push_back( pp );
						}
						else {
							/*#pragma omp critical
							{
								cout << "Leakage2 " << setprecision(32) << pp.x << ";" << pp.y << ";" << pp.z << endl;
								cout << "zmi/zmx " << zmi << ";" << zmx << endl;
							}
							IonLeakage++;*/
							continue;
						}
					}
					//formulating it like this is a little bit less efficient in terms of branch prediction
					//but avoids running two loop solution as shown below
				}
/*
				//MK::the ionpp3_ionpp3idx array allows to trace back the implicit ID/name of the ion
				//if a value in this array reads UINT32MX it indicates that this ion is a guardzone ion an hence
				//has no corresponding reference to an ion_pp3 ion because guardzone ions are copies outside the
				//region the thread is responsible for and guardzone ions are only used to assure a consistent
				//computation of the tessellation!

				//old two loop solution
				if ( Settings::VolumeTessellation != E_NOTESS ) { //tessellation eventually done
					apt_real local_zmi = tessmii.zmi;
					apt_real local_zmx = tessmii.zmx;
					apt_real voro_zmi = tessmii.zmi - tesshalowidth.z;
					apt_real voro_zmx = tessmii.zmx + tesshalowidth.z;
					for( size_t i = 0; i < pn; ++i ) { //pruning most likely not even in halo
						if ( thesepoints->at(i).z < voro_zmi )
							continue;
						if ( thesepoints->at(i).z >= voro_zmx )
							continue;
						//one of mine, either in local region or halo i.e. VALIDATOM or GUARDATOM
						if ( thesepoints->at(i).z >= local_zmi && thesepoints->at(i).z < local_zmx ) { //##MK::handle global ion order
							ionpp3_tess.push_back( p3dm3(thesepoints->at(i).x, thesepoints->at(i).y,
									thesepoints->at(i).z, 0, VALIDZONE_ION, theselabels->at(i)) );
						}
						else { //in guardzone
							ionpp3_tess.push_back( p3dm3(thesepoints->at(i).x, thesepoints->at(i).y,
									thesepoints->at(i).z, 0, GUARDZONE_ION, theselabels->at(i)) );
						}
					}
				}
				for(size_t i = 0; i < pn; ++i) { //z < zmii ? not included : (z < zmxx) ? included : not included
					if ( thesepoints->at(i).z < zmi )
						continue;
					if ( thesepoints->at(i).z < zmx )
						ionpp3.push_back( p3dm1(thesepoints->at(i).x,
								thesepoints->at(i).y, thesepoints->at(i).z, theselabels->at(i)) );
					else
						continue;
				}
*/
			}
		}
		else if ( thesepoints == NULL && theselabels == NULL ) {
			//the bucket is empty but that should not mean that we have a problem
			//remember that the thesepoints/and synchronized theselabel buckets index points spatially
		}
		else {
			stopping("Status of thesepoints is inconsistent with theselabel bucket");
			return false;
		}
	}

	#pragma omp critical
	{
		cout << "Thread " << omp_get_thread_num() << " read_localions ionpp3.size() << " << ionpp3.size() << " ionpp3_tess.size() " << ionpp3_tess.size() << " worldid " << worldid << endl; //"ion leakage is " << IonLeakage << endl;
	}
	return true;
}


inline size_t threadmemory::rectangular_transfer( const apt_xyz pos, const apt_xyz mi, const apt_xyz w, const apt_xyz mx, const size_t nvxl )
{
	//maps from position in continuous space to position in discrete space
	apt_xyz tmpd = (pos - mi) / w;
	size_t tmpsz = static_cast<size_t>(tmpd);
	size_t loc = (tmpsz != nvxl-1) ? 1+tmpsz : nvxl-1;
	//because total bin grid width is nx include left guard zone 1+ .......  +1 right guard zone at nx-1 last addressible bin ix nx-2
	return loc;
}


void threadmemory::binning( sqb const & vxlgrid )
{
	//MK::CALLED FROM WITHIN PARALLEL REGION

	//bins all local ions according to instructed vxlgrid
	ion2bin.clear();

	//##MK::potentially also local reallocate beneficial ?

	size_t ni = ionpp3.size();
	aabb3d limits = vxlgrid.box;
	apt_xyz binwidth = vxlgrid.width;

	for(size_t i = 0; i < ni; ++i) {
		size_t thisbin = -1;

		size_t x = rectangular_transfer( ionpp3.at(i).x, limits.xmi, binwidth, limits.xmx, vxlgrid.nx );
		size_t y = rectangular_transfer( ionpp3.at(i).y, limits.ymi, binwidth, limits.ymx, vxlgrid.ny );
		size_t z = rectangular_transfer( ionpp3.at(i).z, limits.zmi, binwidth, limits.zmx, vxlgrid.nz );

		thisbin = x + y*vxlgrid.nx + z*vxlgrid.nxy;

//cout << i << "\t\t" << thisbin << endl;

		ion2bin.push_back( thisbin );
	}
}


void threadmemory::ion2surfdistance_bvh()
{
	//MK::CALLED FROM WITHIN PARALLEL REGION
	//MK::hint on boost::dynamic_bitset
	//http://www.boost.org/doc/libs/1_64_0/libs/dynamic_bitset/
	// An example of setting and reading some bits. Note that operator[] goes from the least-significant bit at 0 to the most significant
	// bit at size()-1. The operator<< for dynamic_bitset prints the bitset from most-significant to least-significant, since that is
	// the format most people are used to reading. flag all ions close to the boundary
	double tic, toc;
	tic = omp_get_wtime();

	apt_real R = 0.f; //##MK::was Settings::MaximumKernelRadius = 2.0 nm;
	apt_real FarFromSurface = owner->owner->sp->tip.diag() * 2.0;
	apt_real SQRFarFromSurface = SQR(FarFromSurface);
	//SQR to avoid sqrt computations, threadmath.closestPointOnTriangle returns squared distance!

	ion2surf.clear();

	size_t ni = ionpp3.size();
	Tree* thisbvh = owner->owner->surf->bvh;
	vector<tri3d>* rbucket = &owner->owner->surf->tipsurface;
	vector<unsigned int> candidates;
	size_t myworkload = 0;
	for(size_t i = 0; i < ni; ++i) {
		//get my position
		p3d me = p3d( ionpp3.at(i).x, ionpp3.at(i).y, ionpp3.at(i).z );

		//instead of naive O(N^2) approach test for all ions against all triangles
		//prune most triangles by utilizing bounded volume hierarchy (BVH) Rtree
		AABB e_aabb( trpl(me.x-R, me.y-R, me.z-R), trpl(me.x+R, me.y+R, me.z+R) );

		candidates.clear();
		candidates = thisbvh->query( e_aabb ); //MK::assure that this is read only ##MK future improve for better locality multiple Rtrees

		myworkload += candidates.size();

		//if axis-aligned bounding box about triangle does not intersect this probe volume triangle must be farther away than RCutoff
		if ( candidates.size() < 1 ) { //no triangle close, most likely case
			ion2surf.push_back( FarFromSurface );
		}
		else {
			apt_real CurrentBest = SQRFarFromSurface;
			apt_real CurrentValue = SQRFarFromSurface;
			size_t nc = candidates.size();
			for ( size_t c = 0; c < nc; ++c ) {
				CurrentValue = threadmath.closestPointOnTriangle( rbucket->at(candidates.at(c)), me ); //##MK::far field fetch but much fewer candidates in total

				if ( CurrentValue > CurrentBest ) {
					continue; //most likely
				}
				else {
					CurrentBest = CurrentValue;
//cout << "i/c/nc/best " << i << "\t\t" << c << "\t\t" << nc << "\t\t" << CurrentBest << endl;
				}
			}
			//final decision on dSQR against minimum distance required to probe that sphere of radius at least dSQRCutoff does not protrude trihull
			//if dSQR > RSQR, deepinside the tip volume
			//else, ion lays too close to the tipsurface therefore attempting probing a sphere would protrude into vacuum and yield fewer neighbors i.e unbiased higher order statistics
			ion2surf.push_back( CurrentBest );
		}
	} //for each ion

	toc = omp_get_wtime();
	#pragma omp critical
	{
		cout << "Thread " << omp_get_thread_num() << " processed " << ionpp3.size() << " probed " << myworkload << " took " << (toc-tic) << " seconds" << endl;
	}
}


void threadmemory::ion2surfdistance_init( const apt_xyz val )
{
	//MK::CALLED FROM WITHIN PARALLEL REGION
	ion2surf.clear();
	size_t ni = ionpp3.size();
	for(size_t i = 0; i < ni; ++i) {
		ion2surf.push_back( val ); //mentality distance is at least val
	}
	//MK::see that using such looping the implicit order of ion names is the same on ionpp3 and ion2surf
}


bool threadmemory::build_local_kdtree()
{
	//MK::called from within parallel region

	//MK::there is a principal conflict, the input heavy data on owner->owner->owner->rawdata are stored only approximately in order along increasing depth, i.e. negative z direction
	//within an arbitrary plane however the sequence probes almost randomly over the detector space, hence the order in
	//owner->owner->sp->db[] is not expected in any reasonable way memory efficiently organized
	//for all so far conducted tasks this is irrelevant as processing operations such as ranging, finding bounding boxes are essential single time trivial parallel reads
	//spatial range queries though require an as best as possible balance grouping of ion data in logical groups (i.e. tree structures) to remain efficient.
	//MK::hence when building a tree the structure on ionpp3 will be reordered, and should afterwards be once also physically restructured to place ions in the same leaf of the tree
	//adjacent in memory to improve locality once probing the pruned search space
	//MK::we do not store an entire copy of the data in the KDTree but require a one time twice allocation to restructure however
	//if we then were to release the original ionpp3 all previous functions face now a different logical sequence of ions when machining off ionpp3
	//##MK:: as of now: build tree, restructure and keep two copies, future improvement, build tree already in spatial partitioning phase and implement leaner via for instance vector<vector*>

	//##MK::if sufficient threads are used also unsigned int will be sufficient

	double tic, toc;
	tic = omp_get_wtime();

	size_t ioncnt = ionpp3.size();
	vector<p3d> ioninput;
	vector<size_t> permutations;
	try {
		ioninput.reserve(ioncnt);
		permutations.reserve(ioncnt);
		ionpp3_kdtree.reserve(ioncnt);
		ion2surf_kdtree.reserve(ioncnt);
	}
	catch (bad_alloc &threadcroak) {
		string mess = "Allocation error during building threadlocal KDTree in " + to_string(omp_get_thread_num());
		stopping(mess);
		return false;
	}

	for(size_t i = 0; i < ioncnt; ++i) {
		ioninput.push_back( p3d( ionpp3.at(i).x, ionpp3.at(i).y, ionpp3.at(i).z) );
	}

	threadtree = new kd_tree; //this calls the default constructor

	//this really constructs the nodes of the tree
	threadtree->build( ioninput, permutations );

	#pragma omp critical
	{
		//cout << "Thread-local KDTree on " << to_string(omp_get_thread_num()) << " has [" << threadtree->min.x << ";" << threadtree->min.y << ";" << threadtree->min.z << "----" << threadtree->max.x << "] with " << threadtree->nodes.size() << " nodes in total"  << endl;
		//cout << "Thread-local KDTree on " << to_string(omp_get_thread_num()) << " has [" << threadtree->min << ";" << threadtree->max.x << "] with " << threadtree->nodes.size() << " nodes in total"  << endl;
		cout << "Thread-local KDTree on " << to_string(omp_get_thread_num()) << " has " << threadtree->nodes.size() << " nodes in total consuming at least " << threadtree->get_treememory_consumption() << endl;
	}

	threadtree->pack_p3dm1_dist( permutations, ionpp3, ion2surf, ionpp3_kdtree, ion2surf_kdtree );

	//permutation array no longer necessary
	//ioninput as well, both are temporaries will be freed automatically upon exiting subroutine

	if ( threadtree->verify( ionpp3_kdtree ) == false ) {
		string mess = "Threadlocal KDTree for " + to_string(omp_get_thread_num()) + " has overlapping indices!";
		stopping( mess );
		return false;
	}

	toc = omp_get_wtime();
	#pragma omp critical
	{
		cout << "Threadlocal KDTree build on " << to_string(omp_get_thread_num()) << " completed in " << (toc-tic) << " seconds" << endl;
	}

	return true;
}



inline apt_real threadmemory::get_zmi() const
{
	return this->zmi;
}

inline apt_real threadmemory::get_zmx() const
{
	return this->zmx;
}

inline bool threadmemory::me_last() const
{
	return this->melast;
}



inline size_t threadmemory::get_memory_consumption()
{
	size_t bytes = 0;
	bytes += ionpp3.size() * sizeof(p3dm1);
	bytes += ion2bin.size() * sizeof(size_t);
	bytes += ion2surf.size() * sizeof(apt_xyz);

	bytes += ionpp3_kdtree.size() * sizeof(p3dm1);
	bytes += ion2surf_kdtree.size() * sizeof(apt_xyz);

	bytes += threadtree->nodes.size()*sizeof(node);
	//do not account for elementar type locals

	return bytes;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////

decompositor::decompositor()
{
	halothickness = p3d64();
	kdtree_success = false;

	owner = NULL;
	healthy = true;
}

decompositor::~decompositor()
{
	//do not delete owner it is only a backreference
	for(size_t mt = MASTER; mt < db.size(); ++mt) {
		delete db.at(mt);
	}
}


void decompositor::tip_aabb_get_extrema()
{
	double tic = MPI_Wtime();
	//computes opposite extreme points axis-aligned bounding box ion locations

	//parallel version, no SIMD, global variables for reduction
	//aabb3d tipOMP;

/*
	apt_real gxmi = numeric_limits<apt_real>::max();
	apt_real gxmx = numeric_limits<apt_real>::lowest();
	apt_real gymi = numeric_limits<apt_real>::max();
	apt_real gymx = numeric_limits<apt_real>::lowest();
	apt_real gzmi = numeric_limits<apt_real>::max();
	apt_real gzmx = numeric_limits<apt_real>::lowest();

	//##MK::#pragma omp parallel reduction(min:gxmi), reduction(max:gxmx), reduction(min:gymi), reduction(max:gymx), reduction(min:gzmi), reduction(max:gzmx)
	{
		//unsigned long mytid = omp_get_thread_num();
		//unsigned long nthr = omp_get_num_threads();
		for ( size_t b = 0; b < owner->recon->pp3.size(); ++b ) {
			//if ( i_take_care(b, mytid, nthr ) == true ) {
			vector<struct p3d>* bucket = owner->recon->pp3.at(b);
			unsigned int n = bucket->size();
			for (unsigned int i = 0; i < n; ++i) {
				//##MK::change to SIMD at some point
				if ( bucket->at(i).x <= gxmi )	gxmi = bucket->at(i).x;
				if ( bucket->at(i).x >= gxmx )	gxmx = bucket->at(i).x;
				if ( bucket->at(i).y <= gymi )	gymi = bucket->at(i).y;
				if ( bucket->at(i).y >= gymx )	gymx = bucket->at(i).y;
				if ( bucket->at(i).z <= gzmi )	gzmi = bucket->at(i).z;
				if ( bucket->at(i).z >= gzmx )	gzmx = bucket->at(i).z;
			}
		}
		/////} //thread distribute buckets among them
		////#pragma omp critical
		////{
		////	if ( gxmi <= tipOMP.xmi )	tipOMP.xmi = gxmi;
		////	if ( gxmx >= tipOMP.xmx )	tipOMP.xmx = gxmx;
		////	if ( gymi <= tipOMP.ymi )	tipOMP.ymi = gymi;
		////	if ( gymx >= tipOMP.ymx )	tipOMP.ymx = gymx;
		////	if ( gzmi <= tipOMP.zmi )	tipOMP.zmi = gzmi;
		////	if ( gzmx >= tipOMP.zmx )	tipOMP.zmx = gzmx;
		////}
	}
	//after exit parallel region results reduced in global
	aabb3d tipOMP = aabb3d( gxmi, gxmx, gymi, gymx, gzmi, gzmx );
	toc = MPI_Wtime();

	cout << "TipOMP AABB identified " << (toc-tic) << " seconds!" << endl;
*/

//DEBUG sequential fallback
	aabb3d tipSEQ = aabb3d();
	for ( size_t b = 0; b < owner->recon->pp3.size(); ++b ) {
		vector<p3d>* bucket = owner->recon->pp3.at(b);
		if ( bucket != NULL ) {
			size_t n = bucket->size();
			for (size_t i = 0; i < n; ++i) {
				if ( bucket->at(i).x <= tipSEQ.xmi )	tipSEQ.xmi = bucket->at(i).x;
				if ( bucket->at(i).x >= tipSEQ.xmx )	tipSEQ.xmx = bucket->at(i).x;
				if ( bucket->at(i).y <= tipSEQ.ymi )	tipSEQ.ymi = bucket->at(i).y;
				if ( bucket->at(i).y >= tipSEQ.ymx )	tipSEQ.ymx = bucket->at(i).y;
				if ( bucket->at(i).z <= tipSEQ.zmi )	tipSEQ.zmi = bucket->at(i).z;
				if ( bucket->at(i).z >= tipSEQ.zmx )	tipSEQ.zmx = bucket->at(i).z;
			}
		}
	}

	tipSEQ.scale();

	tip = tipSEQ;

	cout << "Worldcoordinate system tip bounding box is " << tip << endl;

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "IonPositionExtrema", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "TipSEQ AABB identified " << (toc-tic) << " seconds!" << endl;

/*
	//numerical difference
	std::cout << "X = " << setprecision(18) << (tipOMP.xmi-tipSEQ.xmi) << " nm" << std::endl;
	std::cout << "X = " << setprecision(18) << (tipOMP.xmx-tipSEQ.xmx) << " nm" << std::endl;
	std::cout << "Y = " << setprecision(18) << (tipOMP.ymi-tipSEQ.ymi) << " nm" << std::endl;
	std::cout << "Y = " << setprecision(18) << (tipOMP.ymx-tipSEQ.ymx) << " nm" << std::endl;
	std::cout << "Z = " << setprecision(18) << (tipOMP.zmi-tipSEQ.zmi) << " nm" << std::endl;
	std::cout << "Z = " << setprecision(18) << (tipOMP.zmx-tipSEQ.zmx) << " nm" << std::endl;
*/

//END DEBUG

	//////##MK::substract guardzone at all sides for safe binning ##MK::consider improving
	////tipOMP.xmi -= TIPAABB_GUARDZONE;	tipOMP.xmx += TIPAABB_GUARDZONE;
	////tipOMP.ymi -= TIPAABB_GUARDZONE;	tipOMP.ymx += TIPAABB_GUARDZONE;
	////tipOMP.zmi -= TIPAABB_GUARDZONE;	tipOMP.zmx += TIPAABB_GUARDZONE;
}


void decompositor::loadpartitioning()
{
	double tic = MPI_Wtime();

	unsigned int PARAPROBE_NUM_THREADS = static_cast<unsigned int>(omp_get_max_threads());	
cout << "PARAPROBE_NUM_THREADS = " << PARAPROBE_NUM_THREADS << endl;
    //into how many snippets approximately to decompose?

cout << "Analyzing point spatial distribution in z for load balancing..." << endl;
	//read z coordinates based on which we spatially decompose
	vector<apt_real> zcoordinates;

	//##MK::for larger datasets (>1.0e9 ions) the effort to sort the list zcoordinates is O(NlogN) may be too way to high, than better random sampling to get approximate quantile locations in z
	//however APT datasets with strong heterogeneity, should always be fully probed, i.e. no random probing

	reconstructor* here = owner->recon;
	size_t nb = here->pp3.size();
	for(size_t b = 0; b < nb; ++b) {
		vector<p3d>* these = here->pp3.at(b);
		if ( these != NULL ) {
			size_t n = these->size();
			for(size_t i = 0; i < n; ++i)
				zcoordinates.push_back( these->at(i).z ); //we partition along the zdirection
		}
	}

cout << zcoordinates.size() << " z-coordinate values extracted, now sorting them..." << endl;

	//full sort ascendingly of zcoordinates, speed up by sampling each n-th only
	sort(zcoordinates.begin(), zcoordinates.end() );

	//when there is a stack of N domains we have to split N-1 times
	vector<apt_real> qqss;
	for(unsigned int tid = 1; tid < PARAPROBE_NUM_THREADS; tid++ ) {
		apt_real thisq = static_cast<apt_real>(tid) / static_cast<apt_real>(PARAPROBE_NUM_THREADS);
		qqss.push_back( thisq );
		//cout << tid << "\t\t" << qqss.back() << endl;
	}
	vector<apt_real> qq = mymath.quantiles_nosort( zcoordinates, qqss );

	for(unsigned int i = 0; i < qq.size(); ++i) {
		cout << "Loadbalancing zmx via z-coordinate quantiles for thread " << i << "\t\t" << qq.at(i) << endl;
	}

	if ( Settings::VolumeTessellation != E_NOTESS ) {
		//compute additional bounds when any tessellation is desired for which eventually halo regions are required
		apt_real lambda = static_cast<apt_real>(zcoordinates.size()) / (tip.xsz*tip.ysz*tip.zsz);
		//z guard is five times the average distance between points
		halothickness.x = 0.f; 	halothickness.y = 0.f;
		apt_real halowidth = Settings::TessellationGuardWidth * pow( (-1.f * log(0.5) * 3.f/(4.f*PI*lambda)), (1.f/3.f) );
		if ( PARAPROBE_NUM_THREADS > 1 ) { //halo in case of one thread only is only AABBINCLUSION_EPSILON wide to prevent leakage of ions at boundary
			halothickness.z = halowidth;
		}
		else { //SINGLETHREADED
			halothickness.z = static_cast<apt_real>(AABBINCLUSION_EPSILON);
		}
		cout << "Threaded " << omp_get_thread_num() << " tessellating halothickness " << halothickness.z << endl;
	}

	//build threadlocal partitioning objects
	for(size_t mt = MASTER; mt < PARAPROBE_NUM_THREADS; mt++) {
		db.push_back(NULL);
		spatialsplits.push_back( p6d64() );
	}


	#pragma omp parallel shared(qq,qqss, zcoordinates) //shared but only read
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads()); //##MK::what is the default behavior of parallel region start with all?
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		p6d64 mii = p6d64();
		mii.zmi = (mt == MASTER) ? (static_cast<apt_real>(zcoordinates.at(0)) - static_cast<apt_real>(AABBINCLUSION_EPSILON)) : qq.at(mt-1);
		mii.zmx = (mt == (nt-1)) ? (static_cast<apt_real>(zcoordinates.back()) + static_cast<apt_real>(AABBINCLUSION_EPSILON)) : qq.at(mt);
		//epsilon environment expansion to catch points exactly on boundary
		//##MK::thereby formally thread local domains at least formally slightly overlap
		//but cells to ion ID mapping is disjoint hence even potential extrem cases double computed VC cells
		//cells are always referred to only once in postprocessable of tessellation result

		bool last = (mt == (nt-1)) ? true : false;

		#pragma omp critical
		{
			spatialsplits.at(mt) = mii;
			cout << "Thread " << omp_get_thread_num() << " working on ["
					<< spatialsplits.at(mt).zmi << ";" << spatialsplits.at(mt).zmx << ") is last? " << last << endl;
		}

		//MK::points belong to mt if their z position in [zmii, zmxx) for which we probe as follows left boundary z < zmii ? not included : (z < zmxx) ? included : not included

		threadmemory* memlocal = NULL;
		try {
			memlocal = new threadmemory; //allocate surplus first-touch
			memlocal->owner = this;

			p3d64 miihalo = halothickness;
			memlocal->init_localmemory( mii, miihalo, last );

			bool status = memlocal->read_localions();

			#pragma omp critical
			{
				db.at(mt) = memlocal;
			}
		}
		catch (bad_alloc &ompexec) {
			string mess = "Thread " + to_string(mt) + " unable to allocate memory region!";
			stopping( mess );
		}
	}

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "Loadpartitioning", APT_BVH, APT_IS_PAR, mm, tic, toc);
	cout << "Decompositor computed a balanced loadpartitioning along Z in " << (toc-tic) << " seconds" << endl;
}


void decompositor::reportpartitioning()
{
	double tic, toc;
	tic = MPI_Wtime();

//report geometry and spatial partitioning of the KDTree
	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".KDTreeNodesAllThreads.csv";
	ofstream plog;
	plog.open( fn.c_str() );

	//plog << "ThreadID;KDTreeNodeID;XMin;YMin;ZMin;XMax;YMax;ZMax;NodeID;SplitDimXYZ012\n";
	plog << "ThreadID;KDTreeNodeID;XMin;XMax;YMin;YMax;ZMin;ZMax;NodeID;SplitDimXYZ012\n"; //more human intuitive layout when splitting boxes along dimensions
	unsigned int nt = db.size();
	for ( unsigned int mt = 0; mt < nt; mt++ ) {
		threadmemory* threadinfo = db.at(mt);
		kd_tree* thistree = threadinfo->threadtree;

		vector<scan_task> boxes;
		thistree->get_allboundingboxes( boxes );

		size_t n = boxes.size();

		for( size_t i = 0; i < n; ++i ) {
			scan_task& t = boxes.at(i);

			plog << mt << ";" << i << ";";
			//plog << boxes[i].bx.min.x << ";" << boxes[i].bx.min.y << ";" << boxes[i].bx.min.z << ";" << boxes[i].bx.max.x << ";" << boxes[i].bx.max.y << ";" << boxes[i].bx.max.z << ";";
			plog << boxes[i].bx.min.x << ";" << boxes[i].bx.max.x << ";" << boxes[i].bx.min.y << ";" << boxes[i].bx.max.y << ";" << boxes[i].bx.min.z << ";" << boxes[i].bx.max.z << ";"; //more human intuitive layout
			plog << t.node << ";" << t.dim << "\n";
		}
	}
	plog.flush();
	plog.close();


//report geometry of threads and meta data
	fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".SpatialDistributionAllThreads.csv";
	ofstream tlog;
	tlog.open( fn.c_str() );

	tlog << "ThreadID;ZMin;ZMax;NIons;KDTreeNodes;MemoryConsumption\n"; //more human intuitive layout when splitting boxes along dimensions
	tlog << ";nm;nm;1;1;Bytes\n";
	tlog << "ThreadID;ZMin;ZMax;NIons;KDTreeNodes;MemoryConsumption\n";
	for ( unsigned int mt = 0; mt < nt; mt++ ) {
		threadmemory* threadinfo = db.at(mt);

		tlog << mt << ";" << threadinfo->get_zmi() << ";" << threadinfo->get_zmx() << ";" << threadinfo->ionpp3_kdtree.size() << ";";
		tlog << threadinfo->threadtree->nodes.size() << ";" << threadinfo->get_memory_consumption() << "\n";
	}

	tlog.flush();
	tlog.close();

	toc = MPI_Wtime();
	cout << "Wrote spatial partitioning to file in " << (toc-tic) << " seconds" << endl;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
surfacer::surfacer()
{
	owner = NULL;
	healthy = true;

	bvh = NULL;

	bitmap1 = NULL;
	bitmap2 = NULL;
	rve = aabb3d();
	grid = sqb();

	string prefix = "TipSurfDebug";
	string formatspec = ".vtk";

	fn_alphashape = prefix + "AlphaShape" + formatspec;
	fn_convexhull = prefix + "ConvexHull" + formatspec;
	fn_marchcube = prefix + "MarchingCube" + formatspec;
}

surfacer::~surfacer()
{
	//dop not delete owner, is only backreference
	//vector<tri3d> clears after itself

	delete [] bitmap1; bitmap1 = NULL;
	delete [] bitmap2; bitmap2 = NULL;

	delete bvh; bvh = NULL;
}


bool surfacer::read_vtk_io_tipsurface( const string vtk_io_fn )
{
	//read triangle set of already existent surface triangulation data
	//##MK::currently the naive format, hence just plain reading off from file
	double tic = MPI_Wtime();

	ifstream vtkfile;
	string vtkline;
	string datapiece;

	size_t nvertices = 0;
	vtkfile.open( vtk_io_fn.c_str(), ifstream::in );
	if ( vtkfile.is_open() == true ) {
		//kick four VTK header lines
		for ( unsigned int i = 0; i < 4; ++i )
			getline( vtkfile, vtkline );

		//read next line to get total number of disjoint triangle vertices
		getline( vtkfile, vtkline );
		if ( vtkline.size() > 0 ) {
			istringstream vtklinestream( vtkline );
			getline( vtklinestream, datapiece, ' ');
			getline( vtklinestream, datapiece, ' '); nvertices = std::stoul( datapiece.c_str() );
		}
		else {
			stopping("Unable to read number of unique triangle vertices from VTK file!");
			vtkfile.close(); return false;
		}

		//read through pushing into tipsurface immediately...
		//thats a benefit of the naive outputting of duplicates
		//##MK::binary file or HDF5
		if ( nvertices % 3 != 0 ) {
			stopping("Number of point data is not an integer multiple of 3, so file is likely faulty!");
			vtkfile.close(); return false;
		}
		try {
			size_t ntriangles = nvertices / 3;
			tipsurface.reserve( ntriangles );
		}
		catch ( bad_alloc &surfexc) {
			stopping("Unable to allocate memory to store all triangles!");
			vtkfile.close(); return false;
		}

		//##MK::optimize with pipe operator
		size_t v = 0;
		tri3d tmp;
		while ( v < nvertices && vtkfile.good() == true ) {
			getline( vtkfile, vtkline );
			istringstream vtklinestream1( vtkline );
			getline(vtklinestream1, datapiece, ' ');	tmp.x1 = stof( datapiece.c_str() );
			getline(vtklinestream1, datapiece, ' ');	tmp.y1 = stof( datapiece.c_str() );
			getline(vtklinestream1, datapiece, ' ');	tmp.z1 = stof( datapiece.c_str() );

			getline( vtkfile, vtkline );
			istringstream vtklinestream2( vtkline );
			getline(vtklinestream2, datapiece, ' ');	tmp.x2 = stof( datapiece.c_str() );
			getline(vtklinestream2, datapiece, ' ');	tmp.y2 = stof( datapiece.c_str() );
			getline(vtklinestream2, datapiece, ' ');	tmp.z2 = stof( datapiece.c_str() );

			getline( vtkfile, vtkline );
			istringstream vtklinestream3( vtkline );
			getline(vtklinestream3, datapiece, ' ');	tmp.x3 = stof( datapiece.c_str() );
			getline(vtklinestream3, datapiece, ' ');	tmp.y3 = stof( datapiece.c_str() );
			getline(vtklinestream3, datapiece, ' ');	tmp.z3 = stof( datapiece.c_str() );

//cout << tmp;

			tipsurface.push_back( tmp );

			v = v + 3; //three lines, points parsed for one triangle, next one
		}

		if ( v != nvertices ) {
			stopping("Unable to parse all triangles!");
			vtkfile.close(); return false;
		}
		else {
			reporting("Parsing triangle input from VTK was successful");
			vtkfile.close();
		}
	}

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "VTKFileInputSurfTriangleHull", APT_IO, APT_IS_SEQ, mm, tic, toc);
	cout << "TipSurface with " << tipsurface.size() << " triangles loaded from file " << vtk_io_fn << " successfully in " << (toc-tic) << " seconds" << endl;
	return true;
}


sqb surfacer::define_binning_grid(const apt_xyz edgelen)
{
	//get dimensions of AABB
	sqb b = sqb();
	aabb3d ti = owner->sp->tip;

	//expand
	apt_real guardwidth = static_cast<apt_real>(edgelen);
	ti.blowup( guardwidth );

	b.nx = static_cast<size_t>( ceil((ti.xmx - ti.xmi) / edgelen) ); b.nx += 2; //MK::one layer guardzone each side
	b.ny = static_cast<size_t>( ceil((ti.ymx - ti.ymi) / edgelen) ); b.ny += 2; //MK::one layer guardzone each side
	b.nz = static_cast<size_t>( ceil((ti.zmx - ti.zmi) / edgelen) ); b.nz += 2; //MK::one layer guardzone each side

	b.nxy = b.nx * b.ny;
	b.nxyz = b.nxy * b.nz;

	b.width = edgelen;
	b.box = ti;

	cout << b << endl;

	return b;
}


bool surfacer::new_pruning_mem()
{
	delete [] bitmap1; bitmap1 = NULL;
	delete [] bitmap2; bitmap2 = NULL;

	try {
		bitmap1 = new bool[grid.nxyz];
		bitmap2 = new bool[grid.nxyz];
	}
	catch (bad_alloc &surfexc) {
		stopping("Allocation of new pruning memory");
		return false;
	}

	return true;
}


void surfacer::del_pruning_mem()
{
	delete [] bitmap1; bitmap1 = NULL;
	delete [] bitmap2; bitmap2 = NULL;
}


void surfacer::seq_binarization()
{
	//initial condition if that bitmap[i] = false \forall i
	for(size_t b = 0; b < grid.nxyz; ++b) {
		bitmap1[b] = false;
	}

	//sequentially identify which bins contain at all data
	size_t nt = owner->sp->db.size();
	threadmemory* thisregion = NULL;
	for(size_t mt = MASTER; mt < nt; mt++) {
		thisregion = owner->sp->db.at(mt);

		size_t n = thisregion->ionpp3.size();
		for ( size_t i = 0; i < n; ++i ) {
			size_t thisbin = thisregion->ion2bin.at(i);

			if ( bitmap1[thisbin] == true ) //visited already
				continue; //most likely case as bins contain usually several ions
			else //bin not yet evaluated
				bitmap1[thisbin] = true;
		}
	}

	//##MK::any cluster with falses
}


bool surfacer::seq_hkbased_pruning()
{
	//now we perform a HoshenKopelman percolation analysis on bitmap1 to find the outer percolating shell of false's culling the tip and at the same time
	//identify any potential percolating cluster or individual arbitrarily arranged flicker noise of individual bins physically inside the tip but
	//with the given resolution too small a bin size that it is likely to find any ion inside (lateral resolution, detector efficiency)
	//knowing that the AABB about the tip was extended by a guardzone surplus the fact that APT tips do not touch the upper corner voxel at the AABB in particular
	//not the one in the guardzone we can now include all embedded flicker and separate from the non-connected culling outer vacuum cluster with id as that of the
	//upper corner guardzone voxel label ID and update bitmap1 excluding the guardzone
	//when we are sure that there are no isolated arbitrarily false voxel cluster inside the tip and we can safely now the Moore based surface patch identification

	//##MK::NUMA memory issues! this analyzer is at the moment running in the memory of the master thread only!

	percAnalyzer* hk = NULL;
	try { hk = new percAnalyzer; }
	catch(bad_alloc &surfexc) {
		stopping( owner->owner->get_rank(), "Unable to allocate HoshenKopelman analyzer instance!");
		return false;
	}

	if ( grid.nx >= UINT32MX || grid.ny >= UINT32MX || grid.nz >= UINT32MX || grid.nxyz >= UINT32MX ) {
		stopping( owner->owner->get_rank(), "Number of bins exceeds UINT32MX implement size_t HK!");
		return false;
	} //now cast safe
	unsigned int nx = static_cast<unsigned int>(grid.nx);
	unsigned int ny = static_cast<unsigned int>(grid.ny);
	unsigned int nz = static_cast<unsigned int>(grid.nz);
	unsigned int nxy = static_cast<unsigned int>(grid.nxy);
	unsigned int nxyz = static_cast<unsigned int>(grid.nxyz);

	//pass binarized ion occupancy bin bitfield to the HK analyzer to perform clustering analysis as one would do for percolation studies

	bool hkstatus = hk->initialize( bitmap1, nx, ny, nz	);

	if ( hkstatus == true ) {
		//t2 = omp_get_wtime(); foo.ProfInitializing = t2 - t1; t1 = omp_get_wtime();
		hkstatus = hk->hoshen_kopelman();
	}

	if ( hkstatus == true ) {
		//t2 = omp_get_wtime();foo.ProfHoshenKopeling = t2 - t1; t1 = omp_get_wtime();
		hkstatus = hk->compactify();
	}

	/*if ( hkstatus == true ) {
		//t2 = omp_get_wtime();foo.ProfCompactifying = t2 - t1;t1 = omp_get_wtime();
		hkstatus = hk->checkLabeling();
	}*/

	if ( hkstatus == true ) {
		hkstatus = hk->determine_clustersize_distr();
	}

	if ( hkstatus == true ) {
		//find label of frontmost top left voxel which will be in the guardzone
		unsigned int tfl_corner_b = 0 + 0*nx + (nz-1)*nxy;
		hkstatus = hk->rebinarize( tfl_corner_b, nxyz, bitmap1 );
		//upon exit bitmap1 is true for bins in vacuum and false for tip
	}

	if ( hkstatus == false ) {
		cout << "ERROR::HoshenKopelman unable to segment cluster!" << endl;
		delete hk; hk = NULL;
		return false;
	}
	else {
		//done with HK bitmap1 is now correctly relabeled and potential flicker eliminated..."
		delete hk; hk = NULL;
		return true;
	}
}


void surfacer::seq_find_surface_adjacent_bins()
{
	//bitmap1 is true if bin of vacuum cluster containing no ion
	//bitmap2 will be set to true only if

	unsigned int nx = static_cast<unsigned int>(grid.nx);
	unsigned int ny = static_cast<unsigned int>(grid.ny);
	unsigned int nz = static_cast<unsigned int>(grid.nz);
	unsigned int nxy = static_cast<unsigned int>(grid.nxy);
	unsigned int nxyz = static_cast<unsigned int>(grid.nxyz);

	//MK::REQUIRED! assume bin still embedded deeply inside tip so surrounded by also deeply embedded bins ones
	for ( unsigned int b = 0; b < nxyz; ++b ) { bitmap2[b] = false; }
	//query in Moore environment which bins have not all Moore neighbors trues (occupied with points)
	//MK::remember that bitmap1 is true if this is a bin belonging to the large cluster contacted to vacuum and false if some bin in the tip even if isolated and empty inside the tip volume this is the trick!
	for ( unsigned int bz = 1; bz < nz-1; ++bz ) { //slightly more complex indexing pattern avoids proximity to boundary branches
		for ( unsigned int by = 1; by < ny-1; ++by ) {
			for ( unsigned int bx = 1; bx < nx-1; ++bx ) {
				unsigned int here = bx + by*nx + bz*nxy;
				unsigned int there = 0;

				if( bitmap1[here] == true ) { //MK::two bitmaps to avoid successive overwriting and data invalidating image buffer1 as input and buffer2 as output

					//a bin in vacuum --- is it surrounded by nonempty bins of the tip?

					//if ( ... == false) asks effectively: me as a bin in vacuum is there a neighbor bin that in the tip volume
					//MK::order optimized for cache reutilization using the implicit indexing x+y*nx+z*nxy

					there = (bx-1)	+	(by-1)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by-1)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true; //MK::order optimized for cache reutilization using the implicit indexing x+y*nx+z*nxy
					there = (bx+1)	+	(by-1)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx-1)	+	(by+0)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by+0)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by+0)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx-1)	+	(by+1)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by+1)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by+1)	*nx	+	(bz-1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;

					there = (bx-1)	+	(by-1)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by-1)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by-1)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx-1)	+	(by+0)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					//exclude myself
					there = (bx+1)	+	(by+0)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx-1)	+	(by+1)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by+1)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by+1)	*nx	+	(bz+0)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;

					there = (bx-1)	+	(by-1)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by-1)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by-1)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx-1)	+	(by+0)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by+0)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by+0)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx-1)	+	(by+1)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+0)	+	(by+1)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
					there = (bx+1)	+	(by+1)	*nx	+	(bz+1)	*nxy;		if(	bitmap1[there] == false ) bitmap2[there] = true;
				}
			} //next bin along +x
		} //next xline along +y
	} //next xyslab along +z


//	//##MK::BEGIN DEBUG result ions in bins for which bitmap2 == true are the candidates
//	for(unsigned int b = 0; b < nxyz; ++b) {
//		if ( bitmap2[b] == true ) {
//cout << "Bitmap2 " << b << " has ions!" << endl;
//		}
//	}
	//##MK::END DEBUG
}

void surfacer::seq_filter_candidates( vector<p3d> & c )
{
	size_t nt = owner->sp->db.size();
	threadmemory* thisregion = NULL;
	for(size_t mt = MASTER; mt < nt; mt++) {
		thisregion = owner->sp->db.at(mt);

		size_t n = thisregion->ionpp3.size();
		for ( size_t i = 0; i < n; ++i ) {
//if (i % 1000 == 0) { cout << "Filtering " << i << endl; }

			size_t thisbin = thisregion->ion2bin.at(i);
			if ( bitmap2[thisbin] == true ) {
				p3dm1 ion = thisregion->ionpp3.at(i);
				c.push_back( p3d(ion.x, ion.y, ion.z) ); //we need only the location of the ion for the tip surface reconstruction
			}
			//##MK::potentially optimization possible
		}
	}

	if ( Settings::IOHKFilteredIons == true ) {
		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".AdvPruningIonCands.vtk";
		positions_vtk( c, fn );
	}


	//get volume and count of bins entirely inside
	unsigned int nxyz = static_cast<unsigned int>(grid.nxyz);

	bool* bitmap3 = NULL;
	try {
		bitmap3 = new bool[nxyz];

	}
	catch (bad_alloc &surfexc) {
		stopping("Allocation of bitmap3 pruning memory");
		//input bitmap buffer not longer required
		delete [] bitmap1; bitmap1 = NULL;
		//finally also bitmap2 obsolete
		delete [] bitmap2; bitmap2 = NULL;
		return;
	}
	//MK::bitmap1 is true if bin belonging to the large cluster contacted to vacuum and false if some bin in the tip even if isolated and empty inside the tip volume this is the trick!
	//MK::bitmap2 by now is true if a bin to consider containing candidate ions to use in tip surface reconstruction
	//hence if bitmap1 == false && (in tip potentially surfacecontact) and bitmap2 == false (excluding surface contact) inside
	unsigned int interiorbins = 0;
	for ( unsigned int b = 0; b < nxyz; ++b ) { //slightly more complex indexing pattern avoids proximity to boundary branches
		if (bitmap1[b] == false && bitmap2[b] == false) {
			bitmap3[b] = true;
			interiorbins++;
		}
		else
			bitmap3[b] = false;
	}

	size_t iontotal = 0;
	nt = owner->sp->db.size();
	thisregion = NULL;
	for(size_t mt = MASTER; mt < nt; mt++) {
		thisregion = owner->sp->db.at(mt);
		size_t n = thisregion->ionpp3.size();
		for ( size_t i = 0; i < n; ++i ) {
			size_t thisbin = thisregion->ion2bin.at(i);
			if ( bitmap3[thisbin] == true ) {
				iontotal++;
			}
			//type specific accounting p3dm1 ion = thisregion->ionpp3.at(i);
		}
	}

	cout << "Ions in interior bins in total " << iontotal << endl;
	cout << "Interior bins in total " << interiorbins << endl;
	cout << "Volume per interior bin (" << grid.width << " nm)^3" << endl;

	//input bitmap buffer not longer required
	//finally also bitmap2 obsolete
	delete [] bitmap3; bitmap3 = NULL;
}



bool surfacer::pruning_candidates( const apt_xyz db, vector<p3d> & cand )
{
	double tic = MPI_Wtime();

	//remove potential previous contetn  ##MK::and eventually even reallocate
	cand.clear();

	//tipvolume AABB gridded in voxel of physical extent db, one voxel guardzone at each side into vacuum
	rve = owner->sp->tip;
	grid = define_binning_grid( db );

	//allocate bitmap to speed up cache local testing of which bins contain ions and which any
	if ( new_pruning_mem() == false ) {
		stopping( "Unable to allocate memory for pruning ion cloud");
		return false;
	}

	//thread-local binning
	#pragma omp parallel
	{
		//unsigned int nt = static_cast<unsigned int>(omp_get_num_threads());
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		owner->sp->db.at(mt)->binning( grid );
	}

	reporting( owner->owner->get_rank(), "Ion point cloud binned ready for HoshenKopelman percolation analysis" );

	seq_binarization();

	bool hksuccess = seq_hkbased_pruning();
	if ( hksuccess == false )
		return false;

	seq_find_surface_adjacent_bins();

	seq_filter_candidates( cand );

	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();

	del_pruning_mem();

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "PruningIonsSurfaceTriangulation", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "\t\tWorker " << owner->owner->get_rank() << " " << setprecision(6) << "cand.size() " << cand.size() << endl;
	cout << "\t\t" << (1.0 - ((double) cand.size()) / ((double) (owner->owner->nevt)))*100.0 << "% ions pruned from the measurement with binning " << grid.width << " nm" << endl;
	cout << "Requiring grid nxyz " << grid.nx << ";" << grid.ny << ";" << grid.nz << " bins " << grid.nxyz << " in total in " << setprecision(18) << (toc-tic)  << " seconds" << endl;

	return true;
}


bool surfacer::alphashape_core( vector<p3d> const & inp )
{

#ifdef UTILIZE_CGAL

	//utilize the CGAL library to compute the alpha shapes
	if ( inp.size() < 1 ) {
		complaining( owner->owner->get_rank(), "There are no candidates to construct an alpha shape from!");
		return false;
	}
	if ( inp.size() >= UINT32MX ) {
		complaining( owner->owner->get_rank(), "computing alpha shapes for more than UINT32MX i.e. 4.2 billion ions not implemented!");
		return false;
	}

	string mess = " computing alpha shape from " + to_string(inp.size()) + " candidate points!";
	reporting( owner->owner->get_rank(), mess);

	double tic = MPI_Wtime();

	Delaunay dt;
	unsigned int ncand = inp.size();
	for ( unsigned int c = 0; c < ncand; ++c ) {
		dt.insert( Point_3( inp.at(c).x, inp.at(c).y, inp.at(c).z) ); //##MK::at(c) to [c]
	}

	double toc = MPI_Wtime();
	memsnapshot mm1 = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "SurfTriangulationDelaunayConstruct", APT_GEO, APT_IS_SEQ, mm1, tic, toc);
	cout << "CGAL::Delaunay data structure with total of " << ncand << " points initialized in " << (toc-tic) << endl;

	tic = MPI_Wtime();

	//compute alpha shape
	Alpha_shape_3 as(dt);
	cout << "Alpha shape will be computed in REGULARIZED mode by default"  << endl;

	//find optimal alpha values
	Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
	Alpha_iterator opt = as.find_optimal_alpha(1);
	cout << "Smallest alpha value to get a solid through data points is " << alpha_solid << endl;
	cout << "Optimal alpha value to get one connected component is " << *opt << endl;

	memsnapshot mm2 = owner->owner->tictoc.get_memoryconsumption();

	//set shape to automatically determined optimum
	if ( Settings::SurfaceAShapeAValue == E_ASHAPE_SMALLEST_SOLID) {
		cout << "Taking the smallest alpha value to get a solid through data points is " << alpha_solid << " for triangulation" << endl;
		as.set_alpha(alpha_solid);
	}
	else if ( Settings::SurfaceAShapeAValue == E_ASHAPE_CGAL_OPTIMAL) {
		cout << "Taking the CGAL optimal value to get a solid through data points is " << *opt << " for triangulation" << endl;
		as.set_alpha(*opt);
	}
	else {
		cout << "Taking the CGAL optimal value to get a solid through data points is " << *opt << " for triangulation" << endl;
		as.set_alpha(alpha_solid);
	}

	cout << "The number of solid components of the CGAL alpha shape is " << as.number_of_solid_components() << endl;

	toc = MPI_Wtime();
	memsnapshot mm3 = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "SurfTriangulationAlphaShapeConstruct", APT_GEO, APT_IS_SEQ, mm3, tic, toc);
	cout << "CGAL::Delaunay-based alpha shape generated and assessed in " << (toc-tic) << " seconds" << endl;

	//extract triangulation
	tic = MPI_Wtime();

	//output triangulation of (optimal) alpha shape
	//in CGAL-4.11 the "outer" surface is what we are looking for, i.e. all facets of kind REGULAR
	//https://stackoverflow.com/questions/15905833/saving-cgal-alpha-shape-surface-mesh
	vector<Alpha_shape_3::Facet> facets;
	as.get_alpha_shape_facets(back_inserter(facets), Alpha_shape_3::REGULAR); // Alpha_shape_3::REGULAR);
	size_t nf = facets.size();

	//##MK::DEBUG
	map<int, int> vunique;
	map<int, int>::iterator which;
	vector<p3d> vuni;
	vector<triref3d> tri;

	for ( size_t f = 0; f < nf; ++f ) {
		//To have a consistent orientation of the facet, always consider an exterior cell
		if ( as.classify( facets[f].first ) != Alpha_shape_3::EXTERIOR ) {
			facets[f] = as.mirror_facet( facets[f] );
		}
		CGAL_assertion( as.classify(facets[f].first) == Alpha_shape_3::EXTERIOR );

		int indices[3] = { (facets[f].second+1)%4, (facets[f].second+2)%4, (facets[f].second+3)%4 };

		//according to the encoding of vertex indices, this is needed to get a consistent orientation
		if ( facets[f].second%2 != 0 ) {
			//collect on the fly disjoint vertices
			for ( unsigned int i = 0; i < 3; ++i ) {
				if ( vunique.find(indices[i]) == vunique.end() ) {
					vuni.push_back( p3d(facets[f].first->vertex(indices[i])->point().x(), facets[f].first->vertex(indices[i])->point().y(), facets[f].first->vertex(indices[i])->point().z()) );
					vunique[indices[i]] = vuni.size()-1;
				}
			}

			tipsurface.push_back( tri3d(
					facets[f].first->vertex(indices[0])->point().x(), facets[f].first->vertex(indices[0])->point().y(), facets[f].first->vertex(indices[0])->point().z(),
					facets[f].first->vertex(indices[1])->point().x(), facets[f].first->vertex(indices[1])->point().y(), facets[f].first->vertex(indices[1])->point().z(),
					facets[f].first->vertex(indices[2])->point().x(), facets[f].first->vertex(indices[2])->point().y(), facets[f].first->vertex(indices[2])->point().z() ) );
			//std::cout << facets[f].first->vertex(indices[0])->point() << ";" << facets[f].first->vertex(indices[1])->point() << ";" << facets[f].first->vertex(indices[2])->point() << std::endl;

			tri.push_back( triref3d(indices[0], indices[1], indices[2]) ); //key maps to position in vuni holding the coordinates
		}
		else {
			swap(indices[0], indices[1]);

			//collect on the fly disjoint vertices
			for ( unsigned int i = 0; i < 3; ++i ) {
				if ( vunique.find(indices[i]) == vunique.end() ) {
					vuni.push_back( p3d(facets[f].first->vertex(indices[i])->point().x(), facets[f].first->vertex(indices[i])->point().y(), facets[f].first->vertex(indices[i])->point().z()) );
					vunique[indices[i]] = vuni.size()-1;
				}
			}

			tipsurface.push_back( tri3d(
				facets[f].first->vertex(indices[0])->point().x(), facets[f].first->vertex(indices[0])->point().y(), facets[f].first->vertex(indices[0])->point().z(),
				facets[f].first->vertex(indices[1])->point().x(), facets[f].first->vertex(indices[1])->point().y(), facets[f].first->vertex(indices[1])->point().z(),
				facets[f].first->vertex(indices[2])->point().x(), facets[f].first->vertex(indices[2])->point().y(), facets[f].first->vertex(indices[2])->point().z() ) );

			tri.push_back( triref3d(indices[0], indices[1], indices[2]) );
		}
	}

//cout << "VUni/Tri.size() " << vuni.size() << "\t\t" << tri.size() << endl;
	toc = MPI_Wtime();
	memsnapshot mm4 = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "SurfTriangulationAlphaShapeTriExtract", APT_GEO, APT_IS_SEQ, mm4, tic, toc);
	cout << "Extracted triangularized tip surface (" << nf << "/" << tri.size() << " triangles in total) in " << (toc-tic) << " seconds" << endl;

	/*
	//##MK::deprecated
	if ( Settings::IOTriangulation == true ) {
		tic = MPI_Wtime();
		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".TipSurface.vtk";
		triangulation_vtk_naive( tipsurface, fn );
		toc = MPI_Wtime();
		memsnapshot mm = memsnapshot();
		owner->owner->tictoc.prof_elpsdtime_and_mem( "VTKFileOutputSurfaceTriangleHull", APT_IO, APT_IS_SEQ, mm, tic, toc);
	}*/

#endif

	return true;
}


bool surfacer::alphashape()
{
	//triangularizes surface of ion point cloud utilizing that feeding all ions in the alpha or
	//poisson shape reconstruction is usually too expensive namely different to convex hull algorithms where most points
	//are rejected very earlier, alpha shape building works differently thereby under pressure for larger point count
	//key idea of advanced pruning for surface extraction is as follows
	//most ions are deep inside the tip, so prune search space first by identifying candidate points closer to the surface
	//one could do so by binning or kernel density averaging over a voxelgrid and check if neighboring voxel contain ions or not
	//albeit this generates often geometrical cases in which volume within the tip volume would be detected
	//close to surface, local voids of ions for instance instead we perform a percolation analysis on successively refined subgrid
	//we build the grid about the tip aabb surplus a one layer thick voxel shell about the aabb
	//by construction then at the upper corner there are no ions
	//we then binarize the bin counts, if counts true, if not false, our upper corner is false
	//then we to percolation analysis (HoshenKopelman) to identify the percolating cluster of false of which our corner is a member
	//the inverse i.e. all other voxel are in the tip, this works as long as there is no continuous crack of falses through the tip
	//lastly we sculpt the adjoining moore voxel to this cluster and take all these points as candidates
	//false cracks will happen numerically only if binning is too fine
	//if binning too coarse we prune too few points, so successive refinement
	//however in total we dont want to spent too much time for pruning
	//we could always downsample randomly the number of points in the sculping phase, eventually arriving
	//at the classical naive random sampling method

	/*
	//##MK::adaptive binning
	for(apt_xyz dbin = Settings::AdvIonPruneBinWidthMax; dbin >= Settings::AdvIonPruneBinWidthMin;
			dbin -= Settings::AdvIonPruneBinWidthIncr) {
		//successively refine binning to probe pruning efficiency
		vector<p3d> probing;
		if ( pruning_candidates( dbin, probing ) == false )
			return false;
		//##MK::discard results, DEBUG
	}
	*/

	//##MK::implement automatic choice of cutoff criterion, as for now set explicit value

	//##MK::use only one binning value
	apt_xyz dbin = 0.5 * (Settings::AdvIonPruneBinWidthMin + Settings::AdvIonPruneBinWidthMax); //binning width in nanometer
	if ( pruning_candidates( dbin, candidates ) == false )
		return false;

	//else ready to CGAL...
	if ( alphashape_core( candidates ) == false )
		return false;
	//else
	return true;
}


bool surfacer::convexhull()
{
	return false;
}


bool surfacer::marchingcubes()
{
	return false;
}


bool surfacer::build_rtree()
{
	const size_t n = tipsurface.size();
	if ( n >= static_cast<size_t>(UINT32MX-1) ) {
		stopping("Implementation of triangle hull supports currently at most UINT32MX individuals!");
		return false;
	}

	unsigned int ntriangles = static_cast<unsigned int>(n);
	aabb3d thistip = owner->sp->tip;
	trpl tipbox( thistip.xmx - thistip.xmi, thistip.ymx - thistip.ymi, thistip.zmx - thistip.zmi );

	try {
		bvh = new class Tree( ntriangles );
	}
	catch (bad_alloc &surfexc) {
		stopping("Unable to allocate memory for BVH in surfacer!");
		return false;
	}

	vector<tri3d>::iterator it;
	apt_real xmi,ymi,zmi,xmx,ymx,zmx;
	unsigned int triid = 0;
	for ( it = tipsurface.begin(); it != tipsurface.end(); ++it ) {
		//##MK::fattening in lohedges/aabbcc code only for purpose of improving rebalancing, however this tree is static,
		//hence no fattening required... xmi -= Settings::MaximumKernelRadius

		xmi = fmin( fmin(it->x1, it->x2), it->x3 );
		xmx = fmax( fmax(it->x1, it->x2), it->x3 );
		ymi = fmin( fmin(it->y1, it->y2), it->y3 );
		ymx = fmax( fmax(it->y1, it->y2), it->y3 );
		zmi = fmin( fmin(it->z1, it->z2), it->z3 );
		zmx = fmax( fmax(it->z1, it->z2), it->z3 );

		//https://github.com/lohedges/aabbcc
		//Here index is a key that is used to create a map between particles and nodes in the AABB tree.
		//The key should be unique to the particle and can take any value between 0 and std::numeric_limits<unsigned int>::max() - 1

		bvh->insertParticle( triid, trpl(xmi, ymi, zmi), trpl(xmx, ymx, zmx) );
		triid++;
	}

	if ( Settings::IOTriangulationBVH == true ) {
		//report all nodes of the resulting tree
		string csvfn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".AABBTree.csv";
		bvh->report_tree( csvfn );
	}

	return true;
}


void surfacer::chop_rtree()
{
	delete bvh; bvh = NULL;
}


void surfacer::report_tipsurface()
{
	double tic = MPI_Wtime();

	vector<unsigned int> wuibuf;
	wuibuf.resize((1+1+3)*tipsurface.size()); //total number of elements 5E6*5*4B so <=100MB
	//+1 first index identifies geometric primitive, here triangle polygon
	//+1 second index tells how many elements to excepttypefor XDMF is XDMF keyword what we see
	size_t i = 0;
	size_t j = 0;
	size_t v = 0;
	for( i = 0; i < tipsurface.size(); i++ ) {
		wuibuf[j+0] = 3; //primitive key
		wuibuf[j+1] = 3; //number of vertices to form the primitive
		wuibuf[j+2] = v+0; //vertex IDs
		wuibuf[j+3] = v+1;
		wuibuf[j+4] = v+2;
		j = j + 5;
		v = v + 3; //##MK::multiple counting
	}
	int status = 0;
	h5iometa ifo = h5iometa( PARAPROBE_SURFRECON_ASHAPE_HULL_TOPO, 5*tipsurface.size(), 1 );
	status = owner->owner->resultsh5Hdl.create_contiguous_matrix_u32le( ifo );
	h5offsets offs = h5offsets( 0, 5*tipsurface.size(), 0, 1, 5*tipsurface.size(), 1);
	status = owner->owner->resultsh5Hdl.write_contiguous_matrix_u32le_hyperslab( ifo, offs, wuibuf );
	wuibuf.clear();

	vector<float> wfbuf;
	wfbuf.resize(3*3*tipsurface.size()); //total number of elements 5E6*9*4B so <=180MB
	i = 0;
	j = 0;
	for( i = 0; i < tipsurface.size(); i++ ) {
		wfbuf[j+0] = tipsurface[i].x1; //reinterleave
		wfbuf[j+1] = tipsurface[i].y1;
		wfbuf[j+2] = tipsurface[i].z1;
		wfbuf[j+3] = tipsurface[i].x2;
		wfbuf[j+4] = tipsurface[i].y2;
		wfbuf[j+5] = tipsurface[i].z2;
		wfbuf[j+6] = tipsurface[i].x3;
		wfbuf[j+7] = tipsurface[i].y3;
		wfbuf[j+8] = tipsurface[i].z3;
		j = j + 9;
	}
	ifo = h5iometa( PARAPROBE_SURFRECON_ASHAPE_HULL_GEOM, 3*tipsurface.size(), 3 );
	status = owner->owner->resultsh5Hdl.create_contiguous_matrix_f32le( ifo );
	offs = h5offsets( 0, 3*tipsurface.size(), 0, 3, 3*tipsurface.size(), 3);
	status = owner->owner->resultsh5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, wfbuf );
	//peak memory consumption
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	wfbuf.clear();

	string xdmffn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".TriangleHull.xdmf";
	owner->owner->xdmfmetaHdl.create_tipsurface_file( xdmffn, tipsurface.size(),
			5*tipsurface.size(), 9*tipsurface.size(), owner->owner->resultsh5Hdl.h5resultsfn );

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "HDF5XDMFReportTriangleHull", APT_IO, APT_IS_SEQ, mm, tic, toc);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
rndlabeler::rndlabeler()
{
	owner = NULL;
	//https://codereview.stackexchange.com/questions/109260/seed-stdmt19937-from-stdrandom-device
	urng.seed( Settings::RndDescrStatsPRNGSeed );
	urng.discard( Settings::RndDescrStatsPRNGDiscard );
	shuffled = false;
	applied = false;
}

rndlabeler::~rndlabeler()
{
	//MK::do not delete owner, is only backreference!
	shuffled = false;
	applied = false;
}

void rndlabeler::reset()
{
	//MK::resetting the class object allows it to reutilize it for both general descriptive statistics as well as a posteriori statistics for clustering
	//MK::do not reset owner
	//do not reset or re-seed PRNG
	initial_iontypid.clear();
	temporary_iontypid.clear();
	shuffled = false;
	applied = false;
}


void rndlabeler::learn_original_types_global()
{
	double tic = MPI_Wtime();

	if ( shuffled == false ) {
		initial_iontypid.clear();
		size_t nions = 0;
		for( size_t mt = MASTER; mt < owner->sp->db.size(); mt++) {
			nions += owner->sp->db.at(mt)->ionpp3_kdtree.size();
		}

		try {
			initial_iontypid.reserve( nions );
		}
		catch (bad_alloc &rndcroak) {
			stopping( "Unable to allocate enough memory to learn original types!" );
			shuffled = false;
			return;
		}

		for( size_t mt = MASTER; mt < owner->sp->db.size(); mt++) {
			vector<p3dm1> & these = owner->sp->db.at(mt)->ionpp3_kdtree;
			size_t ni = these.size();
			for( size_t i = 0; i < ni; ++i) {
				initial_iontypid.push_back( these.at(i).m );
			}
		}
	}
	else {
		complaining("Refusing to randomize iontypes as they currently are randomized already!");
	}

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "RndLabelingLearnOldLabelsGlobal", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "Learning initial iontypes took " << (toc-tic) << " seconds" << endl;
}

void rndlabeler::learn_original_types_cluster()
{

}


void rndlabeler::shuffle_types_mt19937()
{
	double tic = MPI_Wtime();

	//https://meetingcpp.com/blog/items/stdrandom_shuffle-is-deprecated.html
	if ( shuffled == false ) { //shuffle only once
		//build temporary with original iontypes
		temporary_iontypid.clear();
		temporary_iontypid.reserve( initial_iontypid.size() );
		temporary_iontypid = initial_iontypid;

		//randomize order on array via uniform Mersenne twister PRNG shuffling, maintaining total number of ions of specific type
		shuffle(temporary_iontypid.begin(), temporary_iontypid.end(), urng);

		shuffled = true; //shuffled but not applied yet ie overwritten on data containers
	}

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "RndLabelingShuffleTypesMT19937Global", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "Randomizing iontypes took " << (toc-tic) << " seconds" << endl;
}


void rndlabeler::apply_shuffled_types_global()
{
	double tic = MPI_Wtime();

	if ( shuffled == true && applied == false ) {
		size_t nions = 0;
		for( size_t mt = MASTER; mt < owner->sp->db.size(); mt++) {
			nions += owner->sp->db.at(mt)->ionpp3_kdtree.size();
		}
		//do we know type of all ions?
		if ( temporary_iontypid.size() == nions ) { //yes
			size_t implicitid = 0; //MK::order maintaining for each thread, threads orders ascendingly based on threadid
			for( size_t mt = MASTER; mt < owner->sp->db.size(); mt++) {
				vector<p3dm1> & these = owner->sp->db.at(mt)->ionpp3_kdtree;
				size_t ni = these.size();
				for( size_t i = 0; i < ni; ++i) {
					these.at(i).m = temporary_iontypid.at( implicitid + i );
				}
				implicitid += these.size();
			} //proceed with next thread
			applied = true;
		}
		else { //no
			complaining( "We have not consistently learned type of all ions, therefore performing no randomization!");
			applied = false;
			return;
		}
	}
	//else {}
		//either shuffled == true and applied == true, so was shuffled already
		//or shuffled == false with any value of applied, was not even shuffled

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "RndLabelingApplyNewLabelsGlobal", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "Applying randomized iontypes took " << (toc-tic) << " seconds" << endl;
}

void rndlabeler::apply_shuffled_types_cluster()
{

}


void rndlabeler::reset_original_types_global()
{
	double tic = MPI_Wtime();

	//actively resets iontype for all threads again to original values
	if ( applied == true ) {
		//####actively reset iontype for all threads to original physically measured values
		size_t nions = 0;
		for( size_t mt = MASTER; mt < owner->sp->db.size(); mt++) {
			nions += owner->sp->db.at(mt)->ionpp3_kdtree.size();
		}
		//##MK::not safe what about temp and init dissimiarly long? do we know type of all ions?
		if ( initial_iontypid.size() == nions ) { //yes
			size_t implicitid = 0; //MK::order maintaining for each thread, threads orders ascendingly based on threadid
			for( size_t mt = MASTER; mt < owner->sp->db.size(); mt++) {
				vector<p3dm1> & these = owner->sp->db.at(mt)->ionpp3_kdtree;
				size_t ni = these.size();
				for( size_t i = 0; i < ni; ++i) {
					these.at(i).m = initial_iontypid.at( implicitid + i );
				}
				implicitid += these.size();
			} //proceed with next thread
			applied = false;
		}
		else { //no
			stopping( "Reestablishing original labels was not possible, dataset is now inconsistent!");
			return;
		}
	}

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "RndLabelingResetOldLabelsGlobal", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "Resetting to original iontypes took " << (toc-tic) << " seconds" << endl;
}


void rndlabeler::reset_original_types_cluster()
{

}


inline bool rndlabeler::is_shuffled() const
{
	return shuffled;
}

inline bool rndlabeler::is_applied() const
{
	return applied;
}



horderdist::horderdist()
{
	owner = NULL;
	healthy = true;
}


horderdist::~horderdist()
{
}




void horderdist::initialize_descrstat_tasks()
{
	double tic = MPI_Wtime();

	parse_tasks( Settings::DescrStatTasksCode, spatstat_tasks, owner->owner->mypse );

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "DescrStatsInitTasks", APT_UTL, APT_IS_SEQ, mm, tic, toc);
}




inline size_t str2sizet( const string str )
{
	return stoul( str );
}


inline bool isanumber( const string str )
{
    std::string::const_iterator it = str.begin();
    while ( it != str.end() && std::isdigit(*it) )
    	++it;

    return !str.empty() && it == str.end();
}


inline bool SortSizetDescOrder( const size_t & a1, const size_t &a2)
{
	return a1 > a2;
}


void horderdist::compute_generic_spatstat2()
{
	//thread team machines off successively the regions, within each regions all members of the thread team share the load
	//each thread machine off his ions binning in local summary statistics which is afterwards dumped to MASTER thread i critical region
	//threads may require to probe KDTrees of their bottom and top neighbors --- potentially multiple
	//(for very flat z slabs---many threads---) can do so in parallel because KDTrees are queried only in local function call but not written to

	if ( 	Settings::SpatStatDoRDF == false &&
			Settings::SpatStatDo1NN == false &&
			Settings::SpatStatDoRIPK == false &&
			Settings::SpatStatDoMKNN == false &&
			Settings::SpatStatDoNPCorr == false &&
			Settings::SpatStatDoCountNbors == false ) {
		complaining( "No spatial statistics task to do" );
		return;
	}

	double tic = MPI_Wtime();
	cout << "Threadparallel general higher order spatial distribution function..." << endl;

	//avoid that in case of 2-point statistics containers get too large
	if ( Settings::SpatStatDoNPCorr == true ) {
		npc3d test = npc3d(Settings::SpatStatRadiusMax, Settings::SpatStatRadiusIncr, 0.f, false);
		if ( test.get_support().nxyz > CUBE(501) ) {
			complaining( "Computation of 2-point statistics would require too large buffers, skipping!" );
			return;
		}
	}

	vector<histogram> g_res_rdf; //##MK::all threads know Settings so the implicit order of the histograms in globalres is well defined
	vector<histogram> g_res_1nn;
	vector<histogram> g_res_rpk;
	vector<npc3d> g_res_npcorr;
	vector<discrhistogram> g_res_cntnb;

	vector<size_t> DescrStatMKNNCandidates;
	//vector<vector<histogram>> g_res_mknn;
	vector<histogram> g_res_mknn;
	if ( Settings::SpatStatDoRDF == true )
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			g_res_rdf.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
	if ( Settings::SpatStatDo1NN == true )
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			g_res_1nn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
	if ( Settings::SpatStatDoRIPK == true )
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			g_res_rpk.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
	if ( Settings::SpatStatDoNPCorr == true )
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			g_res_npcorr.push_back( npc3d(Settings::SpatStatRadiusMax, Settings::SpatStatRadiusIncr, 0.f, true) ); //bin width in nanometer motivated by finite spatial resolution
	if ( Settings::SpatStatDoCountNbors == true )
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			g_res_cntnb.push_back( discrhistogram() );

	if ( Settings::SpatStatDoMKNN == true ) {
		//parse off candidates from thiscode
		if ( Settings::DescrStatMKNNCode.empty() == false && Settings::DescrStatMKNNCode.find('-') == string::npos ) { //if not empty and containing only nonzero numbers
			istringstream line( Settings::DescrStatMKNNCode );
			string datapiece;
			//find how many semicola if any
			int nsemicola = std::count( Settings::DescrStatMKNNCode.begin(), Settings::DescrStatMKNNCode.end(), ';' );
cout << "--->nsemicola=" << nsemicola << endl;
			for( int i = 0; i <= nsemicola; ++i ) { //interpret numeral kth order values, <= because one value trailing last semicolon
				getline( line, datapiece, ';');
				if ( isanumber( datapiece ) == true ) {
					size_t val = str2sizet( datapiece );
					if ( val < 1 ) { complaining( "Instructing a kth order value smaller than 1" ); return; }
					else { val -= 1; }
					DescrStatMKNNCandidates.push_back( val );
cout << "--->" << val << "<---" << endl;
				}
			}
			//MK::sort in descending order because as such when querying and partially sorting via n_th element first with the highest kth order value
			//the array is in most cases sorted already on [0,kth order element]  and subsequent successive sorts less and less costly if necessary at all
			std::sort( DescrStatMKNNCandidates.begin(), DescrStatMKNNCandidates.end(), SortSizetDescOrder );
cout << "--->After sorting" << endl;
for(auto it = DescrStatMKNNCandidates.begin(); it != DescrStatMKNNCandidates.end(); ++it) cout << "--->" << *it << "<---" << endl;
		}

		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
			//g_res_mknn.push_back( vector<histogram>() );
			for(size_t kth = 0; kth < DescrStatMKNNCandidates.size(); ++kth) {
				//g_res_mknn.at(tsk).push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
				g_res_mknn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
			}
		} //in case of kth order candidates one histogram per candidate and task
	}

	#pragma omp parallel
	{
		double mytic, mytoc;
		int jnt = omp_get_num_threads();
		unsigned int nt = static_cast<unsigned int>(jnt);
		int jmt = omp_get_thread_num();
		unsigned int mt = static_cast<unsigned int>(jmt);

		//all distances squared instead of sqrt computation for efficiency unless nbor objects
		//basic operation is as follows take each ion, if it is sufficiently distant from tip take into consideration

		//set thread-local task-local histograms collecting to improve parallel efficiency
		vector<histogram> m_res_rdf;
		vector<histogram> m_res_1nn;
		vector<histogram> m_res_rpk;
		vector<npc3d> m_res_npcorr; //MK::BE CAREFUL underlying container for the 3d can be very large, so make sure to set always ulimit -s unlimited !
		vector<discrhistogram> m_res_cntnb;
		vector<size_t> m_mknn_cand;
		//vector<vector<histogram>> m_res_mknn; //a vector for the tasks, for each task a vector of histograms for the individual kth order values
		vector<histogram> m_res_mknn;
		if ( Settings::SpatStatDoRDF == true )
			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
				m_res_rdf.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDo1NN == true )
			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
				m_res_1nn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDoRIPK == true )
			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
				m_res_rpk.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDoNPCorr == true )
			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
				m_res_npcorr.push_back( npc3d(Settings::SpatStatRadiusMax, Settings::SpatStatRadiusIncr, 0.f, true) );
		if ( Settings::SpatStatDoCountNbors == true )
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
					m_res_cntnb.push_back( discrhistogram() );
		if ( Settings::SpatStatDoMKNN == true ) {
			//make threadlocal copy of, ##MK::potentially threadlocal copies of task list might also be beneficial to increase locality
			for( size_t i = 0; i < DescrStatMKNNCandidates.size(); ++i)
				m_mknn_cand.push_back( DescrStatMKNNCandidates.at(i) );

			//#pragma omp critical
			//{
			//	cout << "Thread/m_mknn_cand.size() " << jmt << "\t\t" << m_mknn_cand.size() << endl;
			//}

			//init threadlocal results arrays per task per kth each
			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
				//m_res_mknn.push_back( vector<histogram>() );
				for ( size_t kth = 0; kth < m_mknn_cand.size(); ++kth ) {
					//m_res_mknn.at(tsk).push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
					m_res_mknn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
				}
			}
		}

		//#pragma omp critical
		//{
		//	cout << "Thread/m_mknn_cand.size()\t\t" << jmt << "\t\t" << m_mknn_cand.size() << endl;
		//}

		//initialization done
		for( int thr = MASTER; thr < nt; thr++ ) { //thread team is instructed to machine off region by region

			mytic = omp_get_wtime();

			//each member in the team reads details of specific region and initializes threadlocal comparator variables
			threadmemory* currentregiondata = owner->sp->db.at(thr);
			vector<p3dm1> const & theseions = currentregiondata->ionpp3_kdtree;
			vector<apt_xyz> & thesedistances = currentregiondata->ion2surf_kdtree;

			apt_xyz R = Settings::SpatStatRadiusMax;
			apt_xyz RSQR = SQR(R);
			kd_tree* curr_kauri = currentregiondata->threadtree;

			apt_real currdata_zmi = currentregiondata->get_zmi();
			apt_real currdata_zmx = currentregiondata->get_zmx();

			size_t MyIonsCurrRegionConsider = 0; //deep inside the tip
			size_t MyIonsCurrRegionDiscard = 0; //close to the tip

			bool NonNPointCorrelationAnalysesDesired = false;
			if ( 	Settings::SpatStatDoRDF == true ||
					Settings::SpatStatDo1NN == true ||
					Settings::SpatStatDoRIPK == true ||
					Settings::SpatStatDoMKNN == true ||
					Settings::SpatStatDoCountNbors == true )
				NonNPointCorrelationAnalysesDesired = true;

			#pragma omp for schedule(dynamic,1) nowait
			for( size_t i = 0; i < theseions.size(); ++i ) { //all members of thread team grap an ion from the current region process it and continue
				p3dm1 me = theseions.at(i);
				//pre-screening is the ion part of any of the spatstat tasks I want to perform?
				bool considerme = false;
				for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
					if ( considerme == false ) {
						for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
							if ( me.m == kt->second ) {
								considerme = true; break;
							} //at least appears in one task, so no need to check other tasks
						}
					} //continue if not in another task
				}
				/*alternative formulation
				for(auto tsk = spatstat_tasks.begin(); tsk != spatstat_tasks.end(); ++tsk) {
					if ( considerme == false ) {
						for( auto kt = tsk->trgcandidates.begin(); kt != tsk->trgcandidates.end(); ++kt ) {
							if ( me.m == kt->second ) {
								considerme = true; break;
							} //at least appears in one task, so no need to check other tasks
						}
					} //continue if not in another task
				}*/
				if ( thesedistances.at(i) >= RSQR && considerme == true ) {
					MyIonsCurrRegionConsider++;

					//probe local environment up to R, get all neighbors, only for 1NN this is inefficient, for RDF and high k kNN it is required
					//and considering that once computed the neighbors are reutilized for all spatstat tasks and multiple distribution functions this is a superior strategy

					//if no directional values are of interest it is sufficient to report distance and nature of the neighbors only
					vector<nbor> neighbors1dm1; neighbors1dm1.clear();
					//if directional values are of interest, ie. to compute n-point spatial correlation also the position of the neighbor is required
					vector<p3dm1> neighbors3dm1; neighbors3dm1.clear();

					if ( Settings::SpatStatDoNPCorr == false ) { //only 1dm1 information required
						if ( (me.z - R) > currdata_zmi && (me.z + R) < currdata_zmx ) { //we have to probe only in curr_kauri this is the most likely case
							curr_kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors1dm1 );
						}
						else {
							//we have to probe the curr_kauri and trees of neighboring regions, potentially multiple to the top and bottom
							curr_kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors1dm1 );
							//probe top neighbors and climb up, when i am not the topmost thread
							if ( thr < (jnt-1) ) {
								for( int nb = (thr+1); nb < jnt; nb++) { //will eventually climb up to the topmost threadregion, also in practice this will never happen
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors1dm1 );
									//##MK::add break criterion
								}
							}
							//probe bottom neighbors and climb down, when i am not the bottommost thread already
							if ( thr > MASTER ) {
								for( int nb = (thr-1); nb > -1; nb-- ) {
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors1dm1 );
									//##MK::add break
								}
							}
						}
					}
					else { //3dm1 information required
						//pull this information first by scanning the environment
						if ( (me.z - R) > currdata_zmi && (me.z + R) < currdata_zmx ) { //we have to probe only in curr_kauri this is the most likely case
							curr_kauri->range_rball_noclear_nosort_p3dm1( i, theseions, RSQR, neighbors3dm1 );
						}
						else {
							//we have to probe the curr_kauri and trees of neighboring regions, potentially multiple to the top and bottom
							curr_kauri->range_rball_noclear_nosort_p3dm1( i, theseions, RSQR, neighbors3dm1 );
							//probe top neighbors and climb up, when i am not the topmost thread
							if ( thr < (jnt-1) ) {
								for( int nb = (thr+1); nb < jnt; nb++) { //will eventually climb up to the topmost threadregion, also in practice this will never happen
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external_p3dm1( me, nborthreadions, RSQR, neighbors3dm1 );
								}
							}
							//probe bottom neighbors and climb down, when i am not the bottommost thread already
							if ( thr > MASTER ) {
								for( int nb = (thr-1); nb > -1; nb-- ) {
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external_p3dm1( me, nborthreadions, RSQR, neighbors3dm1 );
								}
							}
						}

						//utilize the environmental information to execute n-point correlation function related characterization
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) {
									//are there at all neighbors
									if ( neighbors3dm1.size() > Settings::SpatStatKNNOrder ) { //MK::even getting the zeroth element requires neighbors to contain at least one element!
										//filter target envcandidates before sorting
										vector<npnbor> tmp;
										for( size_t nbb = 0; nbb < neighbors3dm1.size(); ++nbb) {
											for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
												if ( neighbors3dm1.at(nbb).m != jt->second ) { continue; }//type of the only one and desired candidate is different
												else {
													p3dm1 thisnb = neighbors3dm1.at(nbb);
													d3d diffv = d3d( thisnb.x - me.x, thisnb.y - me.y, thisnb.z - me.z );
													apt_xyz dd = sqrt( SQR(diffv.u)+SQR(diffv.v)+SQR(diffv.w) );
													tmp.push_back( npnbor(diffv, dd) ); //, thisnb.m) );
												}
											} //done checking all keywords for tsk
										}
										//having now our envcandidates we make a partial sort to get the n_th element
										//still enough envcandidates to name a KNNOrder?
										if ( tmp.size() > Settings::SpatStatKNNOrder ) { //okay there is no way around O(n) at least partial sorting the envcandidates
											nth_element( tmp.begin(), tmp.begin() + Settings::SpatStatKNNOrder, tmp.end(), SortNPNeighborsForAscDistance );
											npnbor luckyone = tmp.at(Settings::SpatStatKNNOrder);
											m_res_npcorr.at(tsk).add( luckyone ); //.diffvector, 1.f );
										}
									}
									break;
								}
							}
						}

						//potentially apart from npoint correlation functions additional tasks are desired
						//in this case though, reutilize information in p3dm1 but flatten for cache efficiency of preceeding analyses
						//flatten next 3dm1 to 1dm1 if further analysis demanded
						if ( NonNPointCorrelationAnalysesDesired == true ) {
							for( size_t ii = 0; ii < neighbors3dm1.size(); ++ii ) {
								p3dm1 thisone = neighbors3dm1.at(ii);
								apt_xyz distance = sqrt( SQR(thisone.x - me.x) + SQR(thisone.y - me.y) + SQR(thisone.z - me.z) );
								neighbors1dm1.push_back( nbor( distance, thisone.m ) );
							}
						}
						else {
							continue; //##MK::I think here is possibility for continuing with next ion already in case of else
						}
					}


					//use the environment for all NonNPointCorrelation-related descriptive statistics tasks and tasks of this ion
					if ( Settings::SpatStatDoRDF == true ) {
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) { //consider me central ion only if of type required by that specific task tsk but only account for once
									for( size_t nbb = 0; nbb < neighbors1dm1.size(); ++nbb) { //scan cache locality efficiently threadlocal neighbor array for envcandidates
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors1dm1.at(nbb).m != jt->second ) { continue; }
											else { m_res_rdf.at(tsk).add( neighbors1dm1.at(nbb).d ); break; }
										} //done checking all keywords for tsk
									} //done checking all neighboring ions of me for that specific task
									break; //we break, because in cases of multiple central ions we must not account for neighbors multiple times!
								} //done performing task tsk
							} //potentially multiple central ion types request the neighbors to be accounted for
						} //proceed with the next task, reutilize the extracted spatial environment of me again to get histogram of other tasks
					}
					if ( Settings::SpatStatDo1NN == true ) {
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) { //spatstat_tasks.at(tsk).target.second ) {
									apt_real closest = F32MX;
									size_t nborid = SIZETMX;
									for( size_t nbb = 0; nbb < neighbors1dm1.size(); ++nbb) {
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors1dm1.at(nbb).m != jt->second ) { continue; }
											else {
												if ( neighbors1dm1.at(nbb).d > closest ) { continue; }//most likely most ions in Settings::SpatStatRadiusMax
												else { closest = neighbors1dm1.at(nbb).d; nborid = nbb; }
											}
										} //done checking all keywords for tsk
									} //done checking all neighboring ions of me for that specific task

									if ( nborid != SIZETMX )
										m_res_1nn.at(tsk).add( closest ); //if ( closest < 0.021 ) {//	cout << closest << "\t\t" << me.x << ";" << me.y << ";" << me.z << endl; //}
									else
										m_res_1nn.at(tsk).add( R + EPSILON );

									break;
								}
							}
						}
					}
					if ( Settings::SpatStatDoRIPK == true ) {
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) { //consider me central ion only if of type required by that specific task tsk
									for( size_t nbb = 0; nbb < neighbors1dm1.size(); ++nbb) { //scan cache locality efficiently threadlocal neighbor array for envcandidates
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors1dm1.at(nbb).m != jt->second ) { continue; }
											else { m_res_rpk.at(tsk).add( neighbors1dm1.at(nbb).d ); break; }
										} //done checking all keywords for tsk
									} //done checking all neighboring ions of me for that specific task
									break;
								} //done performing task tsk
							}
						} //reutilize the extracted spatial environment of me again to get histogram of other tasks
					}
					if ( Settings::SpatStatDoCountNbors == true ) {
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) {
									//how many neighbors of specific nature to an ion of also a certain nature?
									size_t cnt_of_specific_nature = 0;
									for( size_t nbb = 0; nbb < neighbors1dm1.size(); ++nbb) {
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors1dm1.at(nbb).m != jt->second ) { continue; }//type of the only one and desired candidate is different
											else { ++cnt_of_specific_nature; }
										} //done checking all keywords for tsk
									}
									m_res_cntnb.at(tsk).add( cnt_of_specific_nature, 1.f );
									break;
								}
							}
						}
					}
					if ( Settings::SpatStatDoMKNN == true && m_mknn_cand.size() > 0 ) { //if at least k elements exist at all do n_th element partial sorting for distances
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) {
									//filter target envcandidates before sorting
									vector<nbor> tmp;
									for( size_t nbb = 0; nbb < neighbors1dm1.size(); ++nbb) {
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors1dm1.at(nbb).m != jt->second ) { continue; }//type of the only one and desired candidate is different
											else { tmp.push_back( neighbors1dm1.at(nbb) ); }
										} //done checking all keywords for tsk
									}
									//attempt to process through all kth order candidates
									for(size_t kth = 0; kth < m_mknn_cand.size(); ++kth ) {
										size_t SpatStatCurrKNNOrder = m_mknn_cand.at(kth);
										size_t thishist = tsk*m_mknn_cand.size() + kth;
										if ( tmp.size() > SpatStatCurrKNNOrder ) {
											//okay there is no way around O(n) at least partial sorting the envcandidates
											//MK::given that kth values descrease also time spent in nth_element reduces for successive elements

											nth_element( tmp.begin(), tmp.begin() + SpatStatCurrKNNOrder, tmp.end(), SortNeighborsForAscDistance );

											nbor luckyone = tmp.at(SpatStatCurrKNNOrder);
											m_res_mknn.at(thishist).add( luckyone.d );
											//m_res_mknn.at(tsk).at(kth).add( luckyone.d );
										}
										else {
											m_res_mknn.at(thishist).add( R + EPSILON );
											//m_res_mknn.at(tsk).at(kth).add( R + EPSILON );
										}
										//do not break next kth order candidate given that m_mknn_cand is sorted in descending order n_th element sort will be done
									}
									//break now because we do not account for the neighbors within a single task tsk multiple times if me.m == kt->second is true several times
									break;
								}
							} //analyze next kth order candidate on current task
						} //next mknn task
					}
				} //done with all tasks for this single ion, proceed with the next
				else {
					MyIonsCurrRegionDiscard++;
				}
			} //me working through the queue of the current region

			mytoc = omp_get_wtime();
			#pragma omp critical
			{
				cout << "Thread " << mt << " participated in descriptive spatial statistics of region " << thr << " and considered/discarded " << MyIonsCurrRegionConsider << ";" << MyIonsCurrRegionDiscard << " took " << (mytoc-mytic) << " seconds" << endl;
			}

			#pragma omp barrier
		} //thread team processes the next region

		//thread team members reduce local results into global
		#pragma omp critical
		{
			double mytictic = omp_get_wtime();
			if ( Settings::SpatStatDoRDF == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_rdf.size(); ++tsk) {
					g_res_rdf.at(tsk).cnts_lowest += m_res_rdf.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_rdf.at(tsk).bincount(); ++b)
						g_res_rdf.at(tsk).cnts.at(b) += m_res_rdf.at(tsk).cnts.at(b);
					g_res_rdf.at(tsk).cnts_highest += m_res_rdf.at(tsk).cnts_highest;
				}
			}
			if ( Settings::SpatStatDo1NN == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_1nn.size(); ++tsk) {
					g_res_1nn.at(tsk).cnts_lowest += m_res_1nn.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_1nn.at(tsk).bincount(); ++b)
						g_res_1nn.at(tsk).cnts.at(b) += m_res_1nn.at(tsk).cnts.at(b);
					g_res_1nn.at(tsk).cnts_highest += m_res_1nn.at(tsk).cnts_highest;
				}
			}
			if ( Settings::SpatStatDoRIPK == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_rpk.size(); ++tsk) {
					g_res_rpk.at(tsk).cnts_lowest += m_res_rpk.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_rpk.at(tsk).bincount(); ++b)
						g_res_rpk.at(tsk).cnts.at(b) += m_res_rpk.at(tsk).cnts.at(b);
					g_res_rpk.at(tsk).cnts_highest += m_res_rpk.at(tsk).cnts_highest;
				}
			}
			if ( Settings::SpatStatDoNPCorr == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_npcorr.size(); ++tsk) {
//apt_real dummy = 0.f;
					for( size_t vxl = 0; vxl < m_res_npcorr.at(tsk).cnts.size(); ++vxl) {
//dummy = dummy + m_res_npcorr.at(tsk).cnts.at(vxl);
//cout << dummy << endl;
						g_res_npcorr.at(tsk).cnts[vxl] = g_res_npcorr.at(tsk).cnts[vxl] + m_res_npcorr.at(tsk).cnts[vxl];
					}
				}
			}
			if ( Settings::SpatStatDoCountNbors == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_cntnb.size(); ++tsk) {
					for( size_t b = 0; b < m_res_cntnb.at(tsk).cnts.size(); ++b)
						g_res_cntnb.at(tsk).add( b, m_res_cntnb.at(tsk).cnts.at(b) );
				}
			}
			if ( Settings::SpatStatDoMKNN == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
					for( size_t kth = 0; kth < m_mknn_cand.size(); ++kth) {
						size_t thishist = tsk*m_mknn_cand.size()+kth;
						g_res_mknn.at(thishist).cnts_lowest += m_res_mknn.at(thishist).cnts_lowest;
						for( size_t b = 0; b < m_res_mknn.at(thishist).bincount(); ++b)
							g_res_mknn.at(thishist).cnts.at(b) += m_res_mknn.at(thishist).cnts.at(b);
 						g_res_mknn.at(thishist).cnts_highest += m_res_mknn.at(thishist).cnts_highest;
						/*
						/g_res_mknn.at(tsk).at(kth).cnts_lowest += m_res_mknn.at(tsk).at(kth).cnts_lowest;
						for( size_t b = 0; b < m_res_mknn.at(tsk).at(kth).bincount(); ++b)
							g_res_mknn.at(tsk).at(kth).cnts.at(b) += m_res_mknn.at(tsk).at(kth).cnts.at(b);
						g_res_mknn.at(tsk).at(kth).cnts_highest += m_res_mknn.at(tsk).at(kth).cnts_highest;
						*/
					}
				}
			}
			double mytoctoc = omp_get_wtime();
			cout << "Thread " << mt << " finished general spatstat results reduction took " << (mytoctoc-mytictic) << " seconds" << endl;
		} //end of critical region
	}
	//explicit barrier end of parallel region

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "DescrStatsThreadparallelCompute", APT_PPP, APT_IS_PAR, mm, tic, toc);
	cout << "Computing general spatial statistics completed took " << (toc-tic) << " seconds" << endl;

	tic = MPI_Wtime();

	//report results
	//select which
	if ( Settings::SpatStatDoRDF == true ) {
		for(size_t tsk = 0; tsk < g_res_rdf.size(); ++tsk) {
			string what = "RDF";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( E_RDF, what, whichtarg, whichcand, g_res_rdf.at(tsk) ); //##MK::more elegant cast from enum to long
		}
	}
	if ( Settings::SpatStatDo1NN == true ) {
		for(size_t tsk = 0; tsk < g_res_1nn.size(); ++tsk) {
			string what = "1NN";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( E_NEAREST_NEIGHBOR, what, whichtarg, whichcand, g_res_1nn.at(tsk) ); //##MK
		}
	}
	if ( Settings::SpatStatDoRIPK == true ) {
		for(size_t tsk = 0; tsk < g_res_rpk.size(); ++tsk) {
			string what = "RIPK";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( E_RIPLEYK, what, whichtarg, whichcand, g_res_rpk.at(tsk) ); //##MK
		}
	}
	if ( Settings::SpatStatDoMKNN == true ) {
		for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
			for(size_t kth = 0; kth < DescrStatMKNNCandidates.size(); ++kth ) {
				string what = "MKNN" + to_string(DescrStatMKNNCandidates.at(kth)+1); //+1 to change from C style to intuitive neighbor accounting
				string whichtarg = "";
				for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
					whichtarg += jt->first;
				string whichcand = "";
				for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
					whichcand += jt->first;

				size_t thishist = tsk*DescrStatMKNNCandidates.size()+kth;

				report_apriori_descrstat2( E_MKNN, what, whichtarg, whichcand, g_res_mknn.at(thishist) ); //g_res_mknn.at(tsk).at(kth) );
			}
		}
	}
	if ( Settings::SpatStatDoCountNbors == true ) {
		for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
			string whichtarg = "";
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat3( whichtarg, whichcand, g_res_cntnb.at(tsk) );
		}
	}

	if ( Settings::SpatStatDoNPCorr == true ) {
		//generate 3D grid of bin center positions only once
		bool BinningIsForAllTasksTheSame = true;
		if ( g_res_npcorr.size() > 0 ) {
			size_t NumberOfBinsForAllTasks = g_res_npcorr.back().get_support().nxyz;

			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
				if ( g_res_npcorr.at(tsk).get_support().nxyz == NumberOfBinsForAllTasks)
					continue;
				else {
					BinningIsForAllTasksTheSame = false;
					break;
				}
			}
		}

		if ( BinningIsForAllTasksTheSame == true ) {
			//get bincenter only once
			report_npc3d_bincenters_hdf5( g_res_npcorr.back() );

			//but histogram values for every task
			for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
				string whichtarg = "";
				for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
					whichtarg += jt->first;
				string whichcand = "";
				for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
					whichcand += jt->first;

				//report_apriori_descrstat4(  whichtarg, whichcand, g_res_npcorr.at(tsk) );
				report_npc3d_histvalues_hdf5( whichtarg, whichcand, g_res_npcorr.at(tsk) );
			}
		}
		else {
			//##MK::remains to be implemented
		}
	}

	toc = MPI_Wtime();
	mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ReportingDescrStatsResults", APT_IO, APT_IS_SEQ, mm, tic, toc);
	cout << "Reporting general spatial statistics completed took " << (toc-tic) << " seconds" << endl;
}


void horderdist::compute_generic_spatstat3()
{
	//thread team machines off successively the regions, within each regions all members of the thread team share the load
	//each thread machine off his ions binning in local summary statistics which is afterwards dumped to MASTER thread i critical region
	//threads may require to probe KDTrees of their bottom and top neighbors --- potentially multiple
	//(for very flat z slabs---many threads---) can do so in parallel because KDTrees are queried only in local function call but not written to

	if ( 	Settings::SpatStatDoRDF == false &&
			Settings::SpatStatDo1NN == false &&
			Settings::SpatStatDoRIPK == false &&
			Settings::SpatStatDoMKNN == false &&
			Settings::SpatStatDoNPCorr == false &&
			Settings::SpatStatDoCountNbors == false ) {
		complaining( "No spatial statistics task to do" );
		return;
	}

	double tic = MPI_Wtime();
	cout << "Threadparallel general higher order spatial distribution function..." << endl;

	//avoid that in case of 2-point statistics containers get too large
	if ( Settings::SpatStatDoNPCorr == true ) {
		npc3d test = npc3d(Settings::SpatStatRadiusMax, Settings::SpatStatRadiusIncr, 0.f, false);
		if ( test.get_support().nxyz > CUBE(MAXIMUM_NPC3D_BIN_RESOLUTION) ) {
			complaining( "Computation of 2-point statistics would require too large buffers, skipping!" );
			return;
		}
	}

	vector<histogram> g_res_rdf; //##MK::all threads know Settings so the implicit order of the histograms in globalres is well defined
	vector<histogram> g_res_1nn;
	vector<histogram> g_res_rpk;
	vector<npc3d> g_res_npcorr;
	vector<discrhistogram> g_res_cntnb;

	vector<size_t> DescrStatMKNNCandidates;
	vector<histogram> g_res_mknn;
	size_t ntsk = spatstat_tasks.size();
	if ( Settings::SpatStatDoRDF == true )
		for( size_t tsk = 0; tsk < ntsk; ++tsk )
			g_res_rdf.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
	if ( Settings::SpatStatDo1NN == true )
		for( size_t tsk = 0; tsk < ntsk; ++tsk )
			g_res_1nn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
	if ( Settings::SpatStatDoRIPK == true )
		for( size_t tsk = 0; tsk < ntsk; ++tsk )
			g_res_rpk.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
	if ( Settings::SpatStatDoNPCorr == true )
		for( size_t tsk = 0; tsk < ntsk; ++tsk )
			g_res_npcorr.push_back( npc3d(Settings::SpatStatRadiusMax, Settings::SpatStatRadiusIncr, 0.f, true) ); //bin width in nanometer motivated by finite spatial resolution
	if ( Settings::SpatStatDoCountNbors == true )
		for( size_t tsk = 0; tsk < ntsk; ++tsk )
			g_res_cntnb.push_back( discrhistogram() );

	if ( Settings::SpatStatDoMKNN == true ) {
		//parse off candidates from thiscode
		if ( Settings::DescrStatMKNNCode.empty() == false && Settings::DescrStatMKNNCode.find('-') == string::npos ) { //if not empty and containing only nonzero numbers
			istringstream line( Settings::DescrStatMKNNCode );
			string datapiece;
			//find how many semicola if any
			int nsemicola = std::count( Settings::DescrStatMKNNCode.begin(), Settings::DescrStatMKNNCode.end(), ';' );
cout << "--->nsemicola=" << nsemicola << endl;
			for( int i = 0; i <= nsemicola; ++i ) { //interpret numeral kth order values, <= because one value trailing last semicolon
				getline( line, datapiece, ';');
				if ( isanumber( datapiece ) == true ) {
					size_t val = str2sizet( datapiece );
					if ( val < 1 ) { complaining( "Instructing a kth order value smaller than 1" ); return; }
					else { val -= 1; }
					DescrStatMKNNCandidates.push_back( val );
cout << "--->" << val << "<---" << endl;
				}
			}
			//MK::sort in descending order because as such when querying and partially sorting via n_th element first with the highest kth order value
			//the array is in most cases sorted already on [0,kth order element]  and subsequent successive sorts less and less costly if necessary at all
			std::sort( DescrStatMKNNCandidates.begin(), DescrStatMKNNCandidates.end(), SortSizetDescOrder );
cout << "--->After sorting" << endl;
for(auto it = DescrStatMKNNCandidates.begin(); it != DescrStatMKNNCandidates.end(); ++it) cout << "--->" << *it << "<---" << endl;
		}

		for( size_t tsk = 0; tsk < ntsk; ++tsk ) {
			for(size_t kth = 0; kth < DescrStatMKNNCandidates.size(); ++kth)
				g_res_mknn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		} //in case of kth order candidates one histogram per candidate and task
	}

	#pragma omp parallel
	{
		double mytic, mytoc;
		int jnt = omp_get_num_threads();
		int jmt = omp_get_thread_num();
		unsigned int MaximumNumberOfIontypes = owner->owner->mypse.get_maxtypeid();

		//all distances squared instead of sqrt computation for efficiency unless nbor objects
		//basic operation is as follows take each ion, if it is sufficiently distant from tip take into consideration

		//set thread-local task-local histograms collecting to improve parallel efficiency
		vector<histogram> m_res_rdf;
		vector<histogram> m_res_1nn;
		vector<histogram> m_res_rpk;
		vector<npc3d> m_res_npcorr; //MK::BE CAREFUL underlying container for the 3d can be very large, so make sure to set always ulimit -s unlimited !
		vector<discrhistogram> m_res_cntnb;
		vector<size_t> m_mknn_cand;
		vector<histogram> m_res_mknn;
		if ( Settings::SpatStatDoRDF == true )
			for( size_t tsk = 0; tsk < ntsk; ++tsk )
				m_res_rdf.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDo1NN == true )
			for( size_t tsk = 0; tsk < ntsk; ++tsk )
				m_res_1nn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDoRIPK == true )
			for( size_t tsk = 0; tsk < ntsk; ++tsk )
				m_res_rpk.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDoNPCorr == true )
			for( size_t tsk = 0; tsk < ntsk; ++tsk )
				m_res_npcorr.push_back( npc3d(Settings::SpatStatRadiusMax, Settings::SpatStatRadiusIncr, 0.f, true) );
		if ( Settings::SpatStatDoCountNbors == true )
				for( size_t tsk = 0; tsk < ntsk; ++tsk )
					m_res_cntnb.push_back( discrhistogram() );
		if ( Settings::SpatStatDoMKNN == true ) {
			//make threadlocal copy of, ##MK::potentially threadlocal copies of task list might also be beneficial to increase locality
			for( size_t i = 0; i < DescrStatMKNNCandidates.size(); ++i)
				m_mknn_cand.push_back( DescrStatMKNNCandidates.at(i) );

			//init threadlocal results arrays per task per kth each
			for( size_t tsk = 0; tsk < ntsk; ++tsk ) {
				for ( size_t kth = 0; kth < m_mknn_cand.size(); ++kth )
					m_res_mknn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
			}
		}
		//initialization done

		for( int thr = MASTER; thr < jnt; thr++ ) { //thread team is instructed to machine off region by region

			mytic = omp_get_wtime();

			//each member in the team reads details of specific region and initializes threadlocal comparator variables
			threadmemory* currentregiondata = owner->sp->db.at(thr);
			vector<p3dm1> const & theseions = currentregiondata->ionpp3_kdtree;
			vector<apt_xyz> & thesedistances = currentregiondata->ion2surf_kdtree;

			apt_xyz R = Settings::SpatStatRadiusMax;
			apt_xyz RSQR = SQR(R);
			size_t KNNOrder = Settings::SpatStatKNNOrder;
			kd_tree* curr_kauri = currentregiondata->threadtree;

			apt_real currdata_zmi = currentregiondata->get_zmi();
			apt_real currdata_zmx = currentregiondata->get_zmx();

			size_t MyIonsCurrRegionConsider = 0; //deep inside the tip
			size_t MyIonsCurrRegionDiscard = 0; //close to the tip

			bool NonNPointCorrelationAnalysesDesired = false;
			if ( 	Settings::SpatStatDoRDF == true ||
					Settings::SpatStatDo1NN == true ||
					Settings::SpatStatDoRIPK == true ||
					Settings::SpatStatDoMKNN == true ||
					Settings::SpatStatDoCountNbors == true )
				NonNPointCorrelationAnalysesDesired = true;

			#pragma omp for schedule(dynamic,1) nowait
			for( size_t i = 0; i < theseions.size(); ++i ) { //all members of thread team grap an ion from the current region process it and continue
				p3dm1 me = theseions.at(i);
				//pre-screening is the ion part of any of the spatstat tasks I want to perform?
				bool considerme = false;
				for( auto tskit = spatstat_tasks.begin(); tskit != spatstat_tasks.end(); ++tskit) {
					if ( considerme == false ) {
						for( auto kt = tskit->trgcandidates.begin(); kt != tskit->trgcandidates.end(); ++kt ) {
							if ( me.m == kt->second ) {
								considerme = true; break;
							} //at least appears in one task, so no need to check other tasks
						}
					} //continue if not in another task
				}

				if ( thesedistances.at(i) >= RSQR && considerme == true ) {
					MyIonsCurrRegionConsider++;

					//probe local environment up to R, get all neighbors, only for 1NN this is inefficient, for RDF and high k kNN it is required
					//and considering that once computed the neighbors are reutilized for all spatstat tasks and multiple distribution functions this is a superior strategy

					//if no directional values are of interest it is sufficient to report distance and nature of the neighbors only
					vector<nbor> neighbors1dm1; neighbors1dm1.clear();
					//if directional values are of interest, ie. to compute n-point spatial correlation also the position of the neighbor is required
					vector<p3dm1> neighbors3dm1; neighbors3dm1.clear();

					if ( Settings::SpatStatDoNPCorr == false ) { //only 1dm1 information required
						if ( (me.z - R) > currdata_zmi && (me.z + R) < currdata_zmx ) { //we have to probe only in curr_kauri this is the most likely case
							curr_kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors1dm1 );
						}
						else {
							//we have to probe the curr_kauri and trees of neighboring regions, potentially multiple to the top and bottom
							curr_kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors1dm1 );
							//probe top neighbors and climb up, when i am not the topmost thread
							if ( thr < (jnt-1) ) {
								for( int nb = (thr+1); nb < jnt; nb++) { //will eventually climb up to the topmost threadregion, also in practice this will never happen
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors1dm1 );
									//##MK::add break criterion
								}
							}
							//probe bottom neighbors and climb down, when i am not the bottommost thread already
							if ( thr > MASTER ) {
								for( int nb = (thr-1); nb > -1; nb-- ) {
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors1dm1 );
									//##MK::add break
								}
							}
						}
					}
					else { //3dm1 information required
						//pull this information first by scanning the environment
						if ( (me.z - R) > currdata_zmi && (me.z + R) < currdata_zmx ) { //we have to probe only in curr_kauri this is the most likely case
							curr_kauri->range_rball_noclear_nosort_p3dm1( i, theseions, RSQR, neighbors3dm1 );
						}
						else {
							//we have to probe the curr_kauri and trees of neighboring regions, potentially multiple to the top and bottom
							curr_kauri->range_rball_noclear_nosort_p3dm1( i, theseions, RSQR, neighbors3dm1 );
							//probe top neighbors and climb up, when i am not the topmost thread
							if ( thr < (jnt-1) ) {
								for( int nb = (thr+1); nb < jnt; nb++) { //will eventually climb up to the topmost threadregion, also in practice this will never happen
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external_p3dm1( me, nborthreadions, RSQR, neighbors3dm1 );
								}
							}
							//probe bottom neighbors and climb down, when i am not the bottommost thread already
							if ( thr > MASTER ) {
								for( int nb = (thr-1); nb > -1; nb-- ) {
									threadmemory* nbordata = owner->sp->db.at(nb);
									vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
									kd_tree* nbor_kauri = nbordata->threadtree;
									nbor_kauri->range_rball_noclear_nosort_external_p3dm1( me, nborthreadions, RSQR, neighbors3dm1 );
								}
							}
						}

						//utilize the environmental information to execute n-point correlation function related characterization
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) {
									//are there at all neighbors
									if ( neighbors3dm1.size() > KNNOrder ) { //MK::even getting the zeroth element requires neighbors to contain at least one element!
										//filter target envcandidates before sorting
										vector<npnbor> tmp;
										for( size_t nbb = 0; nbb < neighbors3dm1.size(); ++nbb) {
											for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
												if ( neighbors3dm1[nbb].m != jt->second ) { continue; }//type of the only one and desired candidate is different
												else {
													p3dm1 thisnb = neighbors3dm1[nbb];
													d3d diffv = d3d( thisnb.x - me.x, thisnb.y - me.y, thisnb.z - me.z );
													apt_xyz dd = sqrt( SQR(diffv.u)+SQR(diffv.v)+SQR(diffv.w) );
													tmp.push_back( npnbor(diffv, dd) ); //, thisnb.m) );
												}
											} //done checking all keywords for tsk
										}
										//having now our envcandidates we make a partial sort to get the n_th element
										//still enough envcandidates to name a KNNOrder?
										if ( tmp.size() > KNNOrder ) { //okay there is no way around O(n) at least partial sorting the envcandidates
											nth_element( tmp.begin(), tmp.begin() + KNNOrder, tmp.end(), SortNPNeighborsForAscDistance );
											npnbor luckyone = tmp[KNNOrder];
											m_res_npcorr.at(tsk).add( luckyone );
										}
									}
									break;
								}
							}
						}

						//potentially apart from npoint correlation functions additional tasks are desired
						//in this case reuse information in p3dm1 but flatten for cache efficiency of preceeding analyses
						if ( NonNPointCorrelationAnalysesDesired == true ) {
							for( size_t ii = 0; ii < neighbors3dm1.size(); ++ii ) {
								p3dm1 thisone = neighbors3dm1.at(ii);
								apt_xyz distance = sqrt( SQR(thisone.x - me.x) + SQR(thisone.y - me.y) + SQR(thisone.z - me.z) );
								neighbors1dm1.push_back( nbor( distance, thisone.m ) );
							}
						}
						else {
							continue; //##MK::I think here is possibility for continuing with next ion already in case of else
						}
					}

					//if (  NonNPointCorrelationAnalysesDesired == true ) {
						//in case many different analyses are desired it is useful at a 3*O(n) cost to restructure the array
						vector<io_bounds> reorganizer( MaximumNumberOfIontypes, io_bounds() );
						for( auto iit = neighbors1dm1.begin(); iit != neighbors1dm1.end(); ++iit ) {
	//cout << "MaxNumberOfTypes/m\t\t" << MaximumNumberOfIontypes << "\t\t" << iit->m << endl;
							reorganizer[iit->m].n++;
						}
						for( size_t jjt = 0; jjt < reorganizer.size(); ++jjt ) {
							reorganizer[jjt].s = 0;
							for( size_t trailjjt = 0; trailjjt < jjt; ++trailjjt ) {
								reorganizer[jjt].s += reorganizer[trailjjt].n;
							}
							reorganizer[jjt].e = reorganizer[jjt].s + reorganizer[jjt].n;
						}
						for ( size_t jjt = 0; jjt < reorganizer.size(); ++jjt) {
							reorganizer[jjt].n = 0; //reset to use as offsetter
						}
						vector<apt_xyz> type_ordered_dists( neighbors1dm1.size(), F32MX );
						for( auto iit = neighbors1dm1.begin(); iit != neighbors1dm1.end(); ++iit ) {
	//cout << "Length\t\t" << type_ordered_dists.size() << "\t\t" << (reorganizer[iit->m].s + reorganizer[iit->m].n) << endl;
							type_ordered_dists.at(reorganizer[iit->m].s + reorganizer[iit->m].n) = iit->d;
							reorganizer[iit->m].n++;
						}
						//now type_ordered_dists.size() >> reorganizer.size() size in memory half of neighbors1dm1 so slightly better efficiency no probing of incorrect types for any task!

						//use the environment for all NonNPointCorrelation-related descriptive statistics tasks and tasks of this ion
						if ( Settings::SpatStatDoRDF == true ) {
							for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
								for( auto kt = spatstat_tasks[tsk].trgcandidates.begin(); kt != spatstat_tasks[tsk].trgcandidates.end(); ++kt ) {
									if ( me.m == kt->second ) { //consider me central ion only if of type required by that specific task tsk but only account for once
										for ( auto jt = spatstat_tasks[tsk].envcandidates.begin(); jt != spatstat_tasks[tsk].envcandidates.end(); jt++ ) {
											size_t ival_s = reorganizer[jt->second].s;
											size_t ival_e = reorganizer[jt->second].e;
											for( size_t nbb = ival_s; nbb < ival_e; ++nbb) { //scan cache locality efficientlthreadlocal neighbor array for envcandidates
												m_res_rdf[tsk].add( type_ordered_dists[nbb] );
											}
										} //done checking all neighboring ions of me for that specific task
										break;
									} //done performing task tsk
								} //potentially multiple central ion types request the neighbors to be accounted for
							} //proceed with the next task, reutilize the extracted spatial environment of me again to get histogram of other tasks
						}
	//cout << "RDF" << endl;
						if ( Settings::SpatStatDo1NN == true ) {
							for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
								for( auto kt = spatstat_tasks[tsk].trgcandidates.begin(); kt != spatstat_tasks[tsk].trgcandidates.end(); ++kt ) {
									if ( me.m == kt->second ) { //spatstat_tasks[tsk].target.second ) {
										apt_real closest = F32MX;
										size_t nborid = SIZETMX;
										for ( auto jt = spatstat_tasks[tsk].envcandidates.begin(); jt != spatstat_tasks[tsk].envcandidates.end(); jt++ ) {
											size_t ival_s = reorganizer[jt->second].s;
											size_t ival_e = reorganizer[jt->second].e;
											for( size_t nbb = ival_s; nbb < ival_e; ++nbb) {
												if ( type_ordered_dists[nbb] > closest ) { continue; }//most likely most ions in Settings::SpatStatRadiusMax
												else { closest = type_ordered_dists[nbb]; nborid = nbb; }
											} //done checking all keywords for tsk
										} //done checking all neighboring ions of me for that specific task
										if ( nborid != SIZETMX )
											m_res_1nn[tsk].add( closest ); //if ( closest < 0.021 ) {//	cout << closest << "\t\t" << me.x << ";" << me.y << ";" << me.z << endl; //}
										else
											m_res_1nn[tsk].add( R + EPSILON );

										break;
									}
								}
							}
						}
	//cout << "1NN" << endl;
						if ( Settings::SpatStatDoRIPK == true ) {
							for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
								for( auto kt = spatstat_tasks[tsk].trgcandidates.begin(); kt != spatstat_tasks[tsk].trgcandidates.end(); ++kt ) {
									if ( me.m == kt->second ) { //consider me central ion only if of type required by that specific task tsk
										for ( auto jt = spatstat_tasks[tsk].envcandidates.begin(); jt != spatstat_tasks[tsk].envcandidates.end(); jt++ ) {
											size_t ival_s = reorganizer[jt->second].s;
											size_t ival_e = reorganizer[jt->second].e;
											for( size_t nbb = ival_s; nbb < ival_e; ++nbb) { //scan cache locality efficiently threadlocal neighbor array for envcandidates
												m_res_rpk[tsk].add( type_ordered_dists[nbb] );
											} //done checking all keywords for tsk
										} //done checking all neighboring ions of me for that specific task
										break;
									} //done performing task tsk
								}
							} //reutilize the extracted spatial environment of me again to get histogram of other tasks
						}
						if ( Settings::SpatStatDoCountNbors == true ) {
							for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
								for( auto kt = spatstat_tasks[tsk].trgcandidates.begin(); kt != spatstat_tasks[tsk].trgcandidates.end(); ++kt ) {
									if ( me.m == kt->second ) {
										for ( auto jt = spatstat_tasks[tsk].envcandidates.begin(); jt != spatstat_tasks[tsk].envcandidates.end(); jt++ ) {
											size_t cnt_of_specific_nature = reorganizer[jt->second].e - reorganizer[jt->second].s;
											m_res_cntnb[tsk].add( cnt_of_specific_nature, 1.f );
										}
										break;
									}
								}
							}
						}
						if ( Settings::SpatStatDoMKNN == true && m_mknn_cand.size() > 0 ) { //if at least k elements exist at all do n_th element partial sorting for distances
							for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
								for( auto kt = spatstat_tasks[tsk].trgcandidates.begin(); kt != spatstat_tasks[tsk].trgcandidates.end(); ++kt ) {
									if ( me.m == kt->second ) {
										//filter target envcandidates before sorting
										vector<apt_xyz> tmp;
										for ( auto jt = spatstat_tasks[tsk].envcandidates.begin(); jt != spatstat_tasks[tsk].envcandidates.end(); jt++ ) {
											size_t ival_s = reorganizer[jt->second].s;
											size_t ival_e = reorganizer[jt->second].e;
											for( size_t nbb = ival_s; nbb < ival_e; ++nbb ) {
												tmp.push_back( type_ordered_dists[nbb] );
											} //done checking all keywords for tsk
										}
										//attempt to process through all kth order candidates
										for(size_t kth = 0; kth < m_mknn_cand.size(); ++kth ) {
											size_t SpatStatCurrKNNOrder = m_mknn_cand[kth];
											size_t thishist = tsk*m_mknn_cand.size() + kth;
											if ( tmp.size() > SpatStatCurrKNNOrder ) { //average O(n) time complex partitioning instead of full O(nlogn) sorting
												nth_element( tmp.begin(), tmp.begin() + SpatStatCurrKNNOrder, tmp.end() ); //ascending order sorting by default
												m_res_mknn[thishist].add( tmp[SpatStatCurrKNNOrder] );
											}
											else {
												m_res_mknn[thishist].add( R + EPSILON );
											}
											//do not break next kth order candidate given that m_mknn_cand is sorted in descending order n_th element sort will be done
										}
										//break now because we do not account for the neighbors within a single task tsk multiple times if me.m == kt->second is true several times
										break;
									}
								} //analyze next kth order candidate on current task
							} //next mknn task
						}
					//}
					//else { //only NPCorr desired go to next ion
					//	continue;
					//}
				} //done with all tasks for this single ion, proceed with the next
				else {
					MyIonsCurrRegionDiscard++;
				}
			} //me working through the queue of the current region

			mytoc = omp_get_wtime();
			#pragma omp critical
			{
				cout << "Thread " << jmt << " participated in descriptive spatial statistics of region " << thr << " and considered/discarded " << MyIonsCurrRegionConsider << ";" << MyIonsCurrRegionDiscard << " took " << (mytoc-mytic) << " seconds" << endl;
			}

			#pragma omp barrier
		} //thread team processes the next region

		//thread team members reduce local results into global
		#pragma omp critical
		{
			double mytictic = omp_get_wtime();
			if ( Settings::SpatStatDoRDF == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_rdf.size(); ++tsk) {
					g_res_rdf[tsk].cnts_lowest += m_res_rdf[tsk].cnts_lowest;
					for( size_t b = 0; b < m_res_rdf[tsk].bincount(); ++b)
						g_res_rdf[tsk].cnts.at(b) += m_res_rdf[tsk].cnts.at(b);
					g_res_rdf[tsk].cnts_highest += m_res_rdf[tsk].cnts_highest;
				}
			}
			if ( Settings::SpatStatDo1NN == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_1nn.size(); ++tsk) {
					g_res_1nn[tsk].cnts_lowest += m_res_1nn[tsk].cnts_lowest;
					for( size_t b = 0; b < m_res_1nn[tsk].bincount(); ++b)
						g_res_1nn[tsk].cnts.at(b) += m_res_1nn[tsk].cnts.at(b);
					g_res_1nn[tsk].cnts_highest += m_res_1nn[tsk].cnts_highest;
				}
			}
			if ( Settings::SpatStatDoRIPK == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_rpk.size(); ++tsk) {
					g_res_rpk[tsk].cnts_lowest += m_res_rpk[tsk].cnts_lowest;
					for( size_t b = 0; b < m_res_rpk[tsk].bincount(); ++b)
						g_res_rpk[tsk].cnts.at(b) += m_res_rpk[tsk].cnts.at(b);
					g_res_rpk[tsk].cnts_highest += m_res_rpk[tsk].cnts_highest;
				}
			}
			if ( Settings::SpatStatDoNPCorr == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
					size_t nvxl = m_res_npcorr[tsk].cnts.size();
					for( size_t vxl = 0; vxl < nvxl; ++vxl) {
						g_res_npcorr[tsk].cnts[vxl] = g_res_npcorr[tsk].cnts[vxl] + m_res_npcorr[tsk].cnts[vxl];
					}
				}
			}
			if ( Settings::SpatStatDoCountNbors == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) { //##MK::m_res_cntnb.size(); ++tsk) {
					for( size_t b = 0; b < m_res_cntnb[tsk].cnts.size(); ++b)
						g_res_cntnb[tsk].add( b, m_res_cntnb[tsk].cnts.at(b) );
				}
			}
			if ( Settings::SpatStatDoMKNN == true ) {
				for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
					for( size_t kth = 0; kth < m_mknn_cand.size(); ++kth) {
						size_t thishist = tsk*m_mknn_cand.size()+kth;
						g_res_mknn[thishist].cnts_lowest += m_res_mknn[thishist].cnts_lowest;
						for( size_t b = 0; b < m_res_mknn[thishist].bincount(); ++b)
							g_res_mknn[thishist].cnts.at(b) += m_res_mknn[thishist].cnts.at(b);
 						g_res_mknn[thishist].cnts_highest += m_res_mknn[thishist].cnts_highest;
					}
				}
			}
			double mytoctoc = omp_get_wtime();
			cout << "Thread " << jmt << " finished general spatstat results reduction took " << (mytoctoc-mytictic) << " seconds" << endl;
		} //end of critical region
	}
	//explicit barrier end of parallel region

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "DescrStatsThreadparallelCompute", APT_PPP, APT_IS_PAR, mm, tic, toc);
	cout << "Computing general spatial statistics completed took " << (toc-tic) << " seconds" << endl;

	tic = MPI_Wtime();

	//report results
	//select which
	if ( Settings::SpatStatDoRDF == true ) {
		for(size_t tsk = 0; tsk < g_res_rdf.size(); ++tsk) {
			string what = "RDF";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( E_RDF, what, whichtarg, whichcand, g_res_rdf.at(tsk) ); //##MK::more elegant cast from enum to long
		}
	}
	if ( Settings::SpatStatDo1NN == true ) {
		for(size_t tsk = 0; tsk < g_res_1nn.size(); ++tsk) {
			string what = "1NN";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( E_NEAREST_NEIGHBOR, what, whichtarg, whichcand, g_res_1nn.at(tsk) ); //##MK
		}
	}
	if ( Settings::SpatStatDoRIPK == true ) {
		for(size_t tsk = 0; tsk < g_res_rpk.size(); ++tsk) {
			string what = "RIPK";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( E_RIPLEYK, what, whichtarg, whichcand, g_res_rpk.at(tsk) ); //##MK
		}
	}
	if ( Settings::SpatStatDoMKNN == true ) {
		for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
			for(size_t kth = 0; kth < DescrStatMKNNCandidates.size(); ++kth ) {
				string what = "MKNN" + to_string(DescrStatMKNNCandidates.at(kth)+1); //+1 to change from C style to intuitive neighbor accounting
				string whichtarg = "";
				for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
					whichtarg += jt->first;
				string whichcand = "";
				for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
					whichcand += jt->first;

				size_t thishist = tsk*DescrStatMKNNCandidates.size()+kth;

				report_apriori_descrstat2( E_MKNN, what, whichtarg, whichcand, g_res_mknn.at(thishist) ); //g_res_mknn.at(tsk).at(kth) );
			}
		}
	}
	if ( Settings::SpatStatDoCountNbors == true ) {
		for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
			string whichtarg = "";
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat3( whichtarg, whichcand, g_res_cntnb.at(tsk) );
		}
	}

	if ( Settings::SpatStatDoNPCorr == true ) {
		//generate 3D grid of bin center positions only once
		bool BinningIsForAllTasksTheSame = true;
		if ( g_res_npcorr.size() > 0 ) {
			size_t NumberOfBinsForAllTasks = g_res_npcorr.back().get_support().nxyz;

			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk ) {
				if ( g_res_npcorr.at(tsk).get_support().nxyz == NumberOfBinsForAllTasks)
					continue;
				else {
					BinningIsForAllTasksTheSame = false;
					break;
				}
			}
		}

		if ( BinningIsForAllTasksTheSame == true ) {
			//get bincenter only once
			report_npc3d_bincenters_hdf5( g_res_npcorr.back() );

			//but histogram values for every task
			for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
				string whichtarg = "";
				for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
					whichtarg += jt->first;
				string whichcand = "";
				for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
					whichcand += jt->first;

				//report_apriori_descrstat4(  whichtarg, whichcand, g_res_npcorr.at(tsk) );
				report_npc3d_histvalues_hdf5( whichtarg, whichcand, g_res_npcorr.at(tsk) );
			}
		}
		else {
			//##MK::remains to be implemented
		}
	}

	toc = MPI_Wtime();
	mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ReportingDescrStatsResults", APT_IO, APT_IS_SEQ, mm, tic, toc);
	cout << "Reporting general spatial statistics completed took " << (toc-tic) << " seconds" << endl;
}


void horderdist::report_apriori_descrstat2( const long tsktype, const string whichmetric,
		const string whichtarget, const string againstwhich, histogram & hist )
{
//report geometry and spatial partitioning of the KDTree
	//are we randomized or not?
	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + "." + whichmetric + "." + whichtarget + "." + againstwhich;
	if ( owner->rndmizer->is_shuffled() == true && owner->rndmizer->is_applied() == true )
		fn += ".Rndmized.csv";
	else
		fn += ".Original.csv";

	ofstream sslog;
	sslog.open( fn.c_str() );
	if ( sslog.is_open() == false ) {
		string mess = "Unable to open " + fn;
		stopping(mess);
		return;
	}

	sslog.precision(18);
	double bend = hist.start() + hist.width();

	//global meta information

	//##MK::add global status information here

	if ( tsktype == E_RDF ) { //mode-specific file layout
		sslog << "NumberOfIonsInside" << ";" << owner->binner->metadata.nions_inside << "\n";
		sslog << "VolumeInside" << ";" << owner->binner->metadata.volume_inside << " (nm^3)\n";

		if ( owner->binner->metadata.volume_inside > EPSILON ) { //now division is safnm^3...

			double globaldensity = static_cast<double>(owner->binner->metadata.nions_inside);
			globaldensity /= static_cast<double>(owner->binner->metadata.volume_inside);
			sslog << "InsideAtomicDensityEst" << ";" << globaldensity << "\n";

			double norm = (1.f / globaldensity) * (1.f / ((4.f/3.f)*PI));
			sslog << "NormalizationFactor" << ";" << norm << "\n";

			//following the definitions on page 281ff of B. Gault, M. P. Moody, J. M. Cairney and S. P. Ringer
			//Atom Probe Microscopy, dx.doi.org/10.1007/978-1-4614-3436-8
			//$RDF(r) = \frac{1}{\bar{\rho}} \frac{n_{RDF}(r)}{\frac{4}{3}\pi(r+0.5\Delta r)^3-(r-0.5\Delta r)^3}$
			//mind that n_{RDF}(r) in between we find that accumulated counts estimates RipleyK but not the RDF
			sslog << "r;AccumulatedCount;nRDF(r);ShellVolume;RDF(r)\n";
			sslog << "nm;1;1;nm^3;1\n";
			sslog << "r;AccumulatedCount;nRDF(r);ShellVolume;RDF(r)\n";

			//everything below hist.start() is not of interest to us just report lower tail dump
			sslog << hist.start() << ";" << hist.cnts_lowest << ";;;\n";

			//##MK::according to A. Baddeley

			double half_dr = 0.5*hist.width();

			for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
				double r =  hist.start() + (0.5 + static_cast<double>(b)) * hist.width();
				double diff = (b > 0) ? (hist.report(b) - hist.report(b-1)) : 0.f;
				double sphvol = CUBE(r + half_dr) - CUBE(r - half_dr);
				double rdfval = (sphvol > EPSILON) ? norm*diff/sphvol : 0.f;

				sslog << r << ";" << hist.report(b) << ";" << diff << ";" << sphvol << ";" << rdfval << "\n";
			}
			//everything above hist.end() is also not of interest to us just report upper tail dump
			sslog << "" << ";" << hist.cnts_highest << ";;;\n";
		}
	}
	if ( tsktype == E_NEAREST_NEIGHBOR || tsktype == E_MKNN ) {
		if ( tsktype == E_NEAREST_NEIGHBOR )	sslog << "1NN\n";
		if ( tsktype == E_MKNN )				sslog << whichmetric << "\n";

		sslog << "BinEnd(r);Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n";
		sslog << "nm;1;1;1\n";
		sslog << "BinEnd(r);Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n"; //we report binends and accumulated counts

		//below hist.start() lower tail dump
		sslog << hist.start() << ";" << hist.cnts_lowest << ";;\n";

		//get cumulative sum on [ ) interval
		double cum_sum = 0.f;
		for( unsigned int b = 0; b < hist.bincount(); ++b) { cum_sum += hist.report(b); }

		double sum = 0.f;
		for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
			sum += hist.report(b);
			double csum = (cum_sum > EPSILON) ? (sum / cum_sum) : 0.f;
			sslog << bend << ";" << hist.report(b) << ";" << sum << ";" << csum << "\n";
		}

		//above hist.end() upper tail dump
		sslog << "" << ";" << hist.cnts_highest << ";;\n";
	}
	if ( tsktype == E_RIPLEYK ) {
		sslog << "BinEnd(r);AccumulatedCounts\n";
		sslog << "nm;1\n";
		sslog << "BinEnd(r);AccumulatedCounts\n"; //we report binends and accumulated counts

		//below hist.start() lower tail dump
		sslog << hist.start() << ";" << hist.cnts_lowest << "\n";
		//on[ ) interval
		double sum = 0;
		for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
			sum += hist.report(b);
			sslog << bend << ";" << sum << "\n";
		}
		//above hist.end() upper tail dump
		sslog << "" << ";" << hist.cnts_highest << "\n";
	}

	sslog.flush();
	sslog.close();
}


void horderdist::report_apriori_descrstat3( const string whichtarget, const string againstwhich, discrhistogram & hist )
{
//report geometry and spatial partitioning of the KDTree
	//are we randomized or not?
	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".NBCNTinRMAX." + whichtarget + "." + againstwhich;
	if ( owner->rndmizer->is_shuffled() == true && owner->rndmizer->is_applied() == true )
		fn += ".Rndmized.csv";
	else
		fn += ".Original.csv";

	ofstream sslog;
	sslog.open( fn.c_str() );
	if ( sslog.is_open() == false ) {
		string mess = "Unable to open " + fn;
		stopping(mess);
		return;
	}

	sslog.precision(18);

	sslog << "NeighborCountInRMax\n";
	sslog << "Value;Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n";
	sslog << "1;1;1;1\n";
	sslog << "Value;Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n";

	//get cumulative sum on [ ) interval
	double cum_sum = 0.f;
	for( size_t b = 0; b < hist.cnts.size(); ++b) {
		cum_sum += hist.cnts.at(b);
	}

	double sum = 0.f;
	for( size_t b = 0; b < hist.cnts.size(); ++b ) {
		sum += hist.cnts.at(b);
		double csum = (cum_sum > EPSILON) ? (sum / cum_sum) : 0.f;
		sslog << b << ";" << hist.cnts.at(b) << ";" << sum << ";" << csum << "\n";
	}

	sslog.flush();
	sslog.close();
}


/*
void horderdist::report_apriori_descrstat4( const string whichtarget, const string againstwhich, npc3d & pcf )
{
	sqb pcfinfo = pcf.get_support();

	//##MK::int bounds checks
	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".NPCORR." + whichtarget + "." + againstwhich;
	fn += ".NEDGCUBE." + to_string(pcfinfo.nx);
	if ( owner->rndmizer->is_shuffled() == true && owner->rndmizer->is_applied() == true )
		fn += ".Rndmized.raw";
	else
		fn += ".Original.raw";

	MPI_File msFileHdl;
	MPI_Status msFileStatus;
	//in mpi.h MPI_Offset is defined as an __int64 which is long long, thus we can jump much more than 2^32 directly when unsigned int would be utilized
	MPI_File_open( MPI_COMM_SELF, fn.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &msFileHdl );
	long long totalOffset = 0;
	MPI_File_seek( msFileHdl, totalOffset, MPI_SEEK_SET );

	double* rawdata = NULL;
	try { rawdata = new double[pcfinfo.nxy]; }
	catch (bad_alloc &ioexc) {
		stopping("Unable to allocate memory for writing in descrstats4!");
		return;
	}

	int NX = pcfinfo.nx;
	int NY = pcfinfo.ny;
	int NZ = pcfinfo.nz;
	int NXY = pcfinfo.nxy;

	for( int z = 0; z < NZ; ++z ) {
		for ( int i = 0; i < NXY; ++i ) { rawdata[i] = 0.f; } //debug label
		double dummy = 0.f;
		for ( int y = 0; y < NY; ++y ) {
			for ( int x = 0; x < NX; ++x ) {
				dummy = dummy + pcf.cnts.at(x + y*NX + z*NXY);
				rawdata[x + y*NX] = static_cast<double>(pcf.cnts.at(x + y*NX + z*NXY));
			}
		}
		cout << "Section z " << z << " has " << dummy << " counts" << endl;
		//xy layer at once
		MPI_File_write(msFileHdl, rawdata, NXY, MPI_DOUBLE, &msFileStatus); //implicit advancement of fp
	} //next region z with regions on stacked top of one another in y

	delete [] rawdata; rawdata = NULL;
	MPI_File_close(&msFileHdl); //no Barrier as MPI_COMM_SELF
}
*/


void horderdist::report_npc3d_bincenters_hdf5( npc3d & pcf )
{
	sqb pcfinfo = pcf.get_support();
	vector<size_t> u64buf;
	u64buf.push_back( pcfinfo.nx );
	u64buf.push_back( pcfinfo.ny );
	u64buf.push_back( pcfinfo.nz );

	int status = 0;
	h5iometa ifo = h5iometa( PARAPROBE_DESCRSTATS_NCORR_BINNING, u64buf.size(), 1 );
	status = owner->owner->resultsh5Hdl.create_contiguous_matrix_u64le( ifo );
	h5offsets offs = h5offsets( 0, u64buf.size(), 0, 1, u64buf.size(), 1 );
	cout << "Reporting 2-point spatial statistics histogram bin extent" << endl;
	status = owner->owner->resultsh5Hdl.write_contiguous_matrix_u64le_hyperslab( ifo, offs, u64buf );

	vector<float> wf32buf;
	try {
		wf32buf.reserve( 3*pcfinfo.nxyz ); //xyz are three position values per bin
	}
	catch (std::bad_alloc &croak) {
		return;
	}
	for( size_t z = 0; z < pcfinfo.nz; ++z) {
		for( size_t y = 0; y < pcfinfo.ny; ++y) {
			for( size_t x = 0; x < pcfinfo.nx; ++x) {
				p3d thisone = pcf.get_bincenter(x, y, z);
				wf32buf.push_back( thisone.x );
				wf32buf.push_back( thisone.y );
				wf32buf.push_back( thisone.z );
			}
		}
	}

	//const string thish5fn = owner->owner->resultsh5Hdl.h5resultsfn;
	ifo = h5iometa( PARAPROBE_DESCRSTATS_NCORR_CELLCENTER, (wf32buf.size() / 3), 3 );
	status = owner->owner->resultsh5Hdl.create_contiguous_matrix_f32le( ifo );
	offs = h5offsets( 0, (wf32buf.size() / 3), 0, 3, (wf32buf.size() / 3), 3 );
	cout << "Reporting 2-point spatial statistics histogram bin center positions" << endl;
	status = owner->owner->resultsh5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, wf32buf );
	cout << status << endl;
	wf32buf = vector<float>();
}


void horderdist::report_npc3d_histvalues_hdf5( const string whichtarget, const string againstwhich, npc3d & pcf )
{

	string thiscombination = whichtarget + "_" + againstwhich;
	if ( owner->rndmizer->is_shuffled() == true && owner->rndmizer->is_applied() == true )
		thiscombination += "_Rndm";
	else
		thiscombination += "_Orig";

	//##MK::check if npc3d results for such specific combination exist already

	int status = 0;
	string fwslash = "/";
	string pcfdsetnm = PARAPROBE_DESCRSTATS_NCORR + fwslash + thiscombination;
	cout << "Writing 2-point spatial statistics for __" << thiscombination << "__ to H5 file" << endl;
	h5iometa ifo = h5iometa( pcfdsetnm, pcf.cnts.size(), 1 );
	h5offsets offs = h5offsets( 0, pcf.cnts.size(), 0, 1, pcf.cnts.size(), 1);
	status = owner->owner->resultsh5Hdl.create_contiguous_matrix_u32le( ifo );
	cout << "Creation " << status << endl;
	status = owner->owner->resultsh5Hdl.write_contiguous_matrix_u32le_hyperslab( ifo, offs, pcf.cnts );
	cout << "Write " << status << endl;
}



unsigned int flag_as_guarded( const unsigned int typid, const unsigned int maxtypid )
{
	if ( typid <= maxtypid)
		return (typid + ION_AT_CLUSTER);
	else
		return typid;
}

unsigned int flag_as_normal( const unsigned int typid, const unsigned int maxtypid )
{
	//##MK::existent types
	//20000XX guarded - ION_AT_CLUSTER if this difference remains positive or zero ion is guarded
	//10000XX cluster else it may be clustered
	//00000XX normal
	long thisone = static_cast<long>(typid);
	long val = thisone - static_cast<long>(ION_AT_CLUSTER);
	if ( val >= 0 )
		return val; //-->XX

	val = thisone - static_cast<long>(ION_IN_CLUSTER);
	if ( val >= 0 ) //-->XX
		return val;

	return thisone; //-->XX
}


bool is_clustered( const unsigned int typid, const unsigned int maxtypid )
{
	if ( typid <= (maxtypid + ION_IN_CLUSTER) ) { //on [0,maxtypid+ION_IN_CLUSTER] either cluster or normal
		if ( typid <= maxtypid )
			return false; //normal
		else
			return true; //cluster
	}
	else
		return false; //in guard zone
}

bool is_guardzone( const unsigned int typid, const unsigned int maxtypid )
{
	if ( typid <= (maxtypid + ION_IN_CLUSTER)) //on [0,maxtypid+ION_IN_CLUSTER] either cluster or normal
		return false;
	else
		return true; //in guard zone
}


clustertask::clustertask()
{
	boss = NULL;
	globaltree = NULL;
	tskid = -1;
}

clustertask::~clustertask()
{
	//do not delete boss only backreference
	delete globaltree;
}


void clustertask::initialize()
{
	double tic, toc;
	tic = MPI_Wtime();

	//get information about current task ##MK::deprecated unsigned int central_iontypid = mission.target.second;
	map<string, unsigned int> & targ_iontypid = mission.trgcandidates;
	map<string, unsigned int> & cand_iontypid = mission.envcandidates;
	//p3dm1 tmp = p3dm1();

cout << "MaximumSeparationMethod Task with centraltype ";
for(auto it = targ_iontypid.begin(); it != targ_iontypid.end(); it++) { cout << it->second << " "; }
cout << "\n with types in environment ";
for(auto it = cand_iontypid.begin(); it != cand_iontypid.end(); it++) { cout << it->second << " "; }
cout << endl;

	ions_filtered.clear();
	ion2surf_filtered.clear();

	//stack-local whitelist for better pruning performance
	vector<unsigned int> whitelist;
	for(auto it = targ_iontypid.begin(); it != targ_iontypid.end(); it++ )
		whitelist.push_back( it->second );
	for(auto it = cand_iontypid.begin(); it != cand_iontypid.end(); it++)
		whitelist.push_back( it->second );


	//start filtering out ions, ##MK::make thread parallel
	size_t nt = boss->owner->sp->db.size();
	for( size_t mt = MASTER; mt < nt; mt++) {
		vector<p3dm1> & theseions = boss->owner->sp->db.at(mt)->ionpp3_kdtree;
		vector<apt_xyz> & thesedists = boss->owner->sp->db.at(mt)->ion2surf_kdtree;
		assert( theseions.size() == thesedists.size() );
		size_t ni = theseions.size();
		for ( size_t i = 0; i < ni; ++i) {
			for ( size_t j = 0; j < whitelist.size(); ++j ) {
				if ( theseions.at(i).m != whitelist[j] ) { //most likely different type
					continue;
				}
				else { //not different type but one to include
					ions_filtered.push_back( theseions.at(i) );
					ion2surf_filtered.push_back( thesedists.at(i) );

					//##MK::error handling
					//because types are exclusive ##MK::poly ion species we dont need ot check others in whitelist
					break;
				}
			} //next typid to check if to include or not
		} //next ion of threadregion
	} //next threadregion

	toc = MPI_Wtime();
	cout << "MaximumSeparationMethod filtered out " << ions_filtered.size() << "/" << ion2surf_filtered.size() << " took " << (toc-tic) << " seconds" << endl;
}


dbscanres clustertask::hpdbscan( const apt_real d, const size_t Nmin, const unsigned int runid )
{
	double tic, toc;
	tic = omp_get_wtime();

	//pass to HPDBScan and do clustering analysis
	HPDBParms parameters = HPDBParms( d, 1 ); //##MK::was Nmin, however classical 1NN ie nearest neighbor should be within dmax was Nmin );
	omp_set_num_threads(parameters.threads);

	//##MK::further improvement potential by spreading point data to threadlocal memory, will become important for ccNUMA issues

	HPDBSCAN scanner( ions_filtered );

	toc = omp_get_wtime();
	cout << "\t\tHPDBSCAN parallel initialization took " << (toc-tic) << " seconds" << endl;
	tic = omp_get_wtime();

	cout << "\t\tHPDBSCAN parallel scanning with epsilon " << parameters.epsilon << " and minPoints " << parameters.minPoints << endl;
	scanner.scan(parameters.epsilon, parameters.minPoints);

	toc = omp_get_wtime();
	cout << "\t\tHPDBSCAN parallel scanning took " << (toc-tic) << " seconds" << endl;
	tic = omp_get_wtime();
	//grab an already vxlized representation of the tip and use it to remove bias in the precip size distro
	//by not reporting precipitates that span beyond the tip inside

	//summarize clustering analysis
	dbscanres N = scanner.summarize2( 	mission,
										boss->owner->binner->IsInside,
										boss->owner->binner->vxlgrid,
										ions_filtered,
										boss->owner->owner->mypse.get_maxtypeid(),
										tskid, runid,
										boss->owner->owner->clusth5Hdl );
	N.Nmin = Nmin;
	N.Dmax = d;

	toc = omp_get_wtime();
	cout << "\t\tHPDBSCAN sequential post-processing took " << (toc-tic) << " seconds" << endl;

	//##MK::report based on the clustering analysis which ions with respect to the input ioncloud to exclude
	//in a subsequent descriptive statistical analysis where all clusters and a guarding region about them
	//is cut out -- implemented for Priyanshu Bajaj

	return N;
}

bool clustertask::build_kdtree()
{
	//global KDTreelaying in MASTER thread memory
	double tic, toc;
	tic = omp_get_wtime();

	size_t ioncnt = ions_filtered.size();
	vector<p3d> ioninput;
	vector<size_t> permutations;
	try {
		ioninput.reserve(ioncnt);
		permutations.reserve(ioncnt);
		ions_kdtree.reserve(ioncnt);
		ion2surf_kdtree.reserve(ioncnt);
	}
	catch (bad_alloc &clcroak) {
		stopping( "Allocation error during building global KDTree!" );
		return false;
	}

	for(size_t i = 0; i < ioncnt; ++i) {
		ioninput.push_back( p3d( ions_filtered.at(i).x, ions_filtered.at(i).y, ions_filtered.at(i).z) );
	}

	globaltree = new kd_tree; //this calls the default constructor

	//##MK::error handling

	//this really constructs the nodes of the tree
	globaltree->build( ioninput, permutations );

cout << "Task specific global KDTree has " << globaltree->nodes.size() << " nodes in total consuming at least " << globaltree->get_treememory_consumption() << " Bytes" << endl;

	globaltree->pack_p3dm1_d( permutations, ions_filtered, ion2surf_filtered, ions_kdtree, ion2surf_kdtree );

	//permutation array no longer necessary
	//ioninput as well, both are temporaries will be freed automatically upon exiting subroutine

	if ( globaltree->verify( ions_kdtree ) == false ) {
		stopping( "Task-specific global KDTree for has overlapping indices!");
		return false;
	}

	toc = omp_get_wtime();
	cout << "Task-specific global KDTree successful took " << (toc-tic) << " seconds";
	return true;
}


void clustertask::flag2exclude_guard( const apt_real dmx )
{
	//##MK::if an ion is within spherical ball environment to any clustered ion with respect to a
	//previously conducted inclusion analysis it is also flagged for exclusion

	//avoid revisiting and recursive exclusion
	apt_real Guard = Settings::SpatStatRadiusMax; //##MK::Priyanshu or ClusterDmax current here?
	apt_real SQRGuard = SQR(Guard);
	size_t guarded = 0;
	unsigned int mxtypid = boss->owner->owner->mypse.get_maxtypeid();

	//##MK::performance make thread parallel
	//################################

	size_t ni = ions_kdtree.size();
	for(size_t i = 0; i < ni; ++i) { //go over all ions that are cluster
		if ( is_clustered( ions_kdtree.at(i).m, mxtypid ) == true ) {

			vector<size_t> candidates;
			globaltree->range_rball_noclear_nosort_indices( i, ions_kdtree, SQRGuard, candidates );

			//######MK::collect set of thread local rules by passing the range query function a mxtypid
			for(size_t cand = 0; cand < candidates.size(); ++cand) {
				if ( is_clustered( ions_kdtree.at(cand).m, mxtypid) == true ) {
					continue; //cluster ions remain labeled as cluster ions
				}
				else { //normal ions get relabelled as guardzone about cluster
					ions_kdtree.at(cand).m = flag_as_guarded( ions_kdtree.at(i).m, mxtypid );
					++guarded;
				}
			} //process all candidate ions in environment
		}
	}

	//#####MK::threadlocal critical applying of rules

cout << endl << "Flag2Exclude guardzone " << guarded << endl;

	//now all ions of ions_kdtree which cluster or are close to a cluster ie guarded have flagged ion types,
	//while all others still have their original iontype
}


void clustertask::generic_spatstat( const unsigned int runid )
{
	//now get for this point cloud descriptive spatial statistics
	//we remove bias by inspecting
	//a) only ions close enough from tip surface / pole near features, technically by inspecting distance of ion to outer or internal tip surface
	//b) only "normal-typed" ions i.e. those that are not cluster and at least a spherical guard zone of radius about each cluster
	//get all neighbors within Settings::SpatStatRadius
	//each thread machine off his ions binning in local summary statistics which is afterwards dumped to MASTER thread i critical region
	//threads may require to probe KDTrees of their bottom and top neighbors --- potentially multiple
	//(for very flat z slabs---many threads---) can do so in parallel because KDTrees are queried only in local function call but not written to
	double tic, toc;
	tic = MPI_Wtime();

	cout << "Thread-parallel general descriptive statistics on eroded point set..." << endl;
	histogram globalres = histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax);

	#pragma omp parallel shared(globalres)
	{
		double mytic, mytoc;
		mytic = omp_get_wtime();
		int jnt = omp_get_num_threads(); //unsigned int nt = static_cast<unsigned int>(jnt);
		int jmt = omp_get_thread_num(); //unsigned int mt = static_cast<unsigned int>(jmt);

		//all distances squared for efficiency unless nbor objects
		//basic operation is as follows take each ion, if it is sufficiently distant from tip take into consideration
		vector<p3dm1> & theseions = ions_kdtree;
		vector<apt_xyz> & thesedistances = ion2surf_kdtree;
		size_t ni = theseions.size();
		apt_xyz R = Settings::SpatStatRadiusMax;
		apt_xyz RSQR = SQR(R);
		unsigned int mxtypid = boss->owner->owner->mypse.get_maxtypeid();

		//set thread-local task-local histogram collecting to improve parallel efficiency
		histogram myres = histogram( Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax);

		size_t myworkload = 0;
		size_t mychunk = 0;
		size_t mylost = 0;

		//thread process cooperatively descriptive statistics
		#pragma omp for schedule(dynamic,1) nowait
		for(size_t i = 0; i < ni; ++i) {
			p3dm1 me = theseions.at(i);

			//##MK::the descriptive spatial statistics tasks is the same as the mission of this clustertask
			//e.g. when clustering of Sc-Sc, we get descriptive statistics as of now for only Sc-Sc
			if ( thesedistances.at(i) >= RSQR && is_clustered(me.m, mxtypid) == false && is_guardzone(me.m, mxtypid) == false ) {
				vector<nbor> neighbors;
				neighbors.clear();

				globaltree->range_rball_noclear_nosort( i, theseions, RSQR, neighbors );

				//now that we know all neighbors do something useful with it, sort is not required
				mychunk++;
				myworkload += neighbors.size();
				if ( Settings::SpatialDistributionTask == E_RIPLEYK || Settings::SpatialDistributionTask == E_RDF ) {
					for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
						if ( is_clustered( neighbors.at(nbb).m, mxtypid ) == false )
							myres.add( neighbors.at(nbb).d );
					}
				}
				else if ( Settings::SpatialDistributionTask == E_NEAREST_NEIGHBOR ) {
					apt_real closest = F32MX;
					size_t nborid = numeric_limits<size_t>::max();
					for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
						if ( is_clustered( neighbors.at(nbb).m, mxtypid ) == true )
							continue; //do not consider ##MK
						if ( neighbors.at(nbb).d > closest )
							continue; //most likely most ions in Settings::SpatStatRadiusMax
						closest = neighbors.at(nbb).d;
						nborid = nbb;
					}

					if ( nborid != numeric_limits<size_t>::max() ) //myres.at(tsk).add_nodump( neighbors.at(nbb).d );
						myres.add( closest );
					else
						mylost++;
				}
				else { //all other specific spatial statistics tasks
					continue;
				}
//cout << i << "\t\t" << neighbors.size() << endl;
			} //done checking a single ion
//cout << jmt << "/" << i << endl;
		} //next ion

		//MK::implicit barrier at pragma omp for unless nowait clause

		mytoc = omp_get_wtime();

		#pragma omp critical
		{
			//we a accumulating on the master thread
			globalres.cnts_lowest += myres.cnts_lowest;
			for( size_t b = 0; b < myres.bincount(); ++b) {
				globalres.cnts.at(b) += myres.cnts.at(b);
			}
			globalres.cnts_highest += myres.cnts_highest;

			cout << "Thread " << omp_get_thread_num() << " finished general spatstat " << mychunk << "/" << mylost << " took " << (mytoc-mytic) << " seconds" << endl;
		}

	} //end of parallel region

	//report results
	string what = "";
	if ( Settings::SpatialDistributionTask == E_RDF )	what = "RDF";
	else if ( Settings::SpatialDistributionTask == E_NEAREST_NEIGHBOR ) what = "1NN";
	else if ( Settings::SpatialDistributionTask == E_RIPLEYK ) what = "RIPK";
	else what = "";

	string whichtarg = "";
	for( auto jt = mission.trgcandidates.begin(); jt != mission.trgcandidates.end(); jt++)
		whichtarg += jt->first;

	string whichcand = "";
	for( auto jt = mission.envcandidates.begin(); jt != mission.envcandidates.end(); jt++)
		whichcand += jt->first;

	report_aposteriori_descrstat( what, whichtarg, whichcand, runid, globalres );

	toc = MPI_Wtime();
	cout << "Computing aposteriori clustering general spatial statistics completed took " << (toc-tic) << " seconds" << endl;
}


void clustertask::report_aposteriori_descrstat( const string whichmetric,
			const string whichtarget, const string againstwhich,
			const unsigned int rid, histogram & hist )
{
	double tic, toc;
	tic = MPI_Wtime();

	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + "." + whichmetric;
	fn += "." + whichtarget + "." + againstwhich + ".TskID." + to_string(tskid) + ".RunID." + to_string(rid) + ".PostClust.csv";

	ofstream sslog;
	sslog.open( fn.c_str() );
	sslog.precision(18);


	//##MK::improve here for all cases
	sslog << "BinEnd(r);Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n";
	sslog << "nm;1;1;1\n";
	sslog << "BinEnd(r);Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n"; //we report binends and accumulated counts
	apt_real bend = hist.start() + hist.width();

	//below hist.start() lower tail dump
	sslog << hist.start() << ";" << hist.cnts_lowest << ";;\n";

	//get cumulative sum on [ ) interval
	apt_real cum_sum = 0.f;
	for( unsigned int b = 0; b < hist.bincount(); ++b) { cum_sum += hist.report(b); }

	apt_real sum = 0.f;
	for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
		sum += hist.report(b);
		apt_real csum = (cum_sum > EPSILON) ? (sum / cum_sum) : 0.f;
		sslog << bend << ";" << hist.report(b) << ";" << sum << ";" << csum << "\n";
	}

	//above hist.end() upper tail dump
	sslog << "" << ";" << hist.cnts_highest << ";;\n";

	sslog.flush();
	sslog.close();

	toc = MPI_Wtime();
}


void clustertask::chop_kdtree()
{
	if ( globaltree != NULL ) {
		delete globaltree; globaltree = NULL;
	}
}


void clustertask::unflag_all()
{
	//resetting labels to physical value for next iteration of dmax
	//vector<p3dm1> ions_filtered;		//local copy of guys to work on, copy requires memory but allows to process more task characteristics specifically
	size_t ni = ions_filtered.size();
	unsigned int maxtypid = boss->owner->owner->mypse.get_maxtypeid();
	for(size_t i = 0; i < ni; ++i) {
		ions_filtered.at(i).m = flag_as_normal( ions_filtered.at(i).m, maxtypid );
	}
}

void clustertask::reset()
{
	//resetting internal containers for next iteration of dmax
	unflag_all();

	//MK::do not clear or delete **_filtered
	ions_kdtree.clear();
	ion2surf_kdtree.clear();

	chop_kdtree();
}



vxlizer::vxlizer()
{
	owner = NULL;

	vxlgrid = sqb();
	metadata = occupancy();

	IsVacuum = NULL;
	IsSurface = NULL;
	IsInside = NULL;
}


vxlizer::~vxlizer()
{
	//MK::do not delete owner only backreference!
	delete [] IsVacuum; IsVacuum = NULL;
	delete [] IsSurface; IsSurface = NULL;
	delete [] IsInside; IsInside = NULL;
}


void vxlizer::characterize_binning()
{
	//MK::characterize binning, do so partially thread-parallel
	size_t NumberOfIonsSurface = 0;
	size_t NumberOfIonsInside = 0;

	#pragma omp parallel //by default all global variables are shared across threads
	{
		size_t myNumberOfIonsSurface = 0; //every thread processes local results
		size_t myNumberOfIonsInside = 0;
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		vector<size_t> & myone = owner->sp->db.at(mt)->ion2bin;
		for(auto it = myone.begin(); it != myone.end(); ++it) {
			if ( IsInside[*it] == true ) {
				myNumberOfIonsInside++;
				continue;
			}
			if ( IsSurface[*it] == true ) {
				myNumberOfIonsSurface++;
			}
		}

		#pragma omp critical
		{
			NumberOfIonsSurface += myNumberOfIonsSurface;
			NumberOfIonsInside += myNumberOfIonsInside;
		}
	}
	//implicit barrier

	//##MK::sequential accounting of binary containers is usually too few worker to start threading
	unsigned int nxyz = static_cast<unsigned int>(vxlgrid.nxyz);

	size_t VacuumBins = 0;
	size_t SurfaceBins = 0;
	size_t InsideBins = 0;
	for ( unsigned int b = 0; b < nxyz; ++b ) {
		if ( IsVacuum[b] == true )
			VacuumBins++;
		if ( IsSurface[b] == true )
			SurfaceBins++;
		if ( IsInside[b] == true )
			InsideBins++;
	}

	apt_real InsideVolume = static_cast<apt_real>(InsideBins);
	InsideVolume *= CUBE(vxlgrid.width); //in nm^3

	metadata = occupancy( vxlgrid.nxyz, VacuumBins, SurfaceBins, InsideBins,
			NumberOfIonsSurface, NumberOfIonsInside, InsideVolume );

	cout << metadata << endl;
}


void vxlizer::identify_inside_bins()
{
	//get volume and count of bins entirely inside
	unsigned int nxyz = static_cast<unsigned int>(vxlgrid.nxyz);

	//MK::IsVacuum is true if bin belonging to the large cluster contacted to vacuum and false if some bin in the tip even if isolated and empty inside the tip volume this is the trick!
	//MK::IsSurface by now is true if a bin to consider containing candidate ions to use in tip surface reconstruction
	//hence if IsVacuum == false && (in tip potentially surfacecontact) and IsSurface == false (excluding surface contact) inside
	for ( unsigned int b = 0; b < nxyz; ++b ) {
		if (IsVacuum[b] == false && IsSurface[b] == false)
			IsInside[b] = true;
		else
			IsInside[b] = false;
	}
}


void vxlizer::identify_surface_adjacent_bins()
{
	//IsVacuum is true if bin is part of the vacuum cluster enclosing the tip, these are bins that do not containing an ion
	unsigned int nx = static_cast<unsigned int>(vxlgrid.nx);
	unsigned int ny = static_cast<unsigned int>(vxlgrid.ny);
	unsigned int nz = static_cast<unsigned int>(vxlgrid.nz);
	unsigned int nxy = static_cast<unsigned int>(vxlgrid.nxy);
	unsigned int nxyz = static_cast<unsigned int>(vxlgrid.nxyz);

	//MK::assume first each bin is not part of the surface, i.e. initialize...
	for ( unsigned int b = 0; b < nxyz; ++b )
		IsSurface[b] = false;

	//now reset IsSurface to true only if any of the bins' Moore neighbors IsVacuum[i] == false...
	//query in Moore environment which bins have not all Moore neighbors trues ie. are not deep in tip (occupied with points)
	//MK::remember that IsVacuum is true if this is a bin belonging to the large cluster contacted to vacuum and false only if it is some bin in the tip even if isolated and empty inside the tip volume this is the trick!
	//start querying at 1 because guardzone is definately not in tip surplus we need to probe non-periodic Moore neighbors
	for ( unsigned int bz = 1; bz < nz-1; ++bz ) {
		for ( unsigned int by = 1; by < ny-1; ++by ) {
			for ( unsigned int bx = 1; bx < nx-1; ++bx ) {
				unsigned int here = bx + by*nx + bz*nxy;
				unsigned int there = 0;

				if( IsVacuum[here] == true ) { //MK::two bitmaps to avoid successive overwriting and data invalidating image buffer1 as input and buffer2 as output

					//a bin in vacuum --- is it surrounded by nonempty bins i.e. bins that obviously intrude the tip volume?

					//MK::mind that access order such to improve cache temporal and spatial locality
					there = (bx-1)	+	(by-1)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by-1)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true; //MK::order optimized for cache reutilization using the implicit indexing x+y*nx+z*nxy
					there = (bx+1)	+	(by-1)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx-1)	+	(by+0)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by+0)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by+0)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx-1)	+	(by+1)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by+1)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by+1)	*nx	+	(bz-1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;

					there = (bx-1)	+	(by-1)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by-1)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by-1)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx-1)	+	(by+0)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					//exclude myself
					there = (bx+1)	+	(by+0)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx-1)	+	(by+1)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by+1)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by+1)	*nx	+	(bz+0)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;

					there = (bx-1)	+	(by-1)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by-1)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by-1)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx-1)	+	(by+0)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by+0)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by+0)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx-1)	+	(by+1)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+0)	+	(by+1)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
					there = (bx+1)	+	(by+1)	*nx	+	(bz+1)	*nxy;		if(	IsVacuum[there] == false ) IsSurface[there] = true;
				}
			} //next bin along +x
		} //next xline along +y
	} //next xyslab along +z
}


bool vxlizer::identify_vacuum_bins()
{
	//now we perform a HoshenKopelman percolation analysis on IsVacuum to find the outer percolating shell of false's culling the tip and at the same time
	//identify any potential percolating cluster or individual arbitrarily arranged flicker noise of individual bins physically inside the tip but
	//with the given resolution too small a bin size that it is likely to find any ion inside (lateral resolution, detector efficiency)
	//knowing that the AABB about the tip was extended by a guardzone surplus the fact that APT tips do not touch the upper corner voxel at the AABB in particular
	//not the one in the guardzone we can now include all embedded flicker and separate from the non-connected culling outer vacuum cluster with id as that of the
	//upper corner guardzone voxel label ID and update bitmap1 excluding the guardzone
	//when we are sure that there are no isolated arbitrarily false voxel cluster inside the tip and we can safely now the Moore based surface patch identification

	//##MK::NUMA memory issues! this analyzer is at the moment running in the memory of the master thread only!
	//##empirical evidence so far even for 0.5nm binning and 1 billion ions though showed that total core time taken
	//is too insignificant to be a prime target for code optimization...

	percAnalyzer* hk = NULL;
	try { hk = new percAnalyzer; }
	catch(bad_alloc &surfexc) {
		stopping( owner->owner->get_rank(), "Unable to allocate HoshenKopelman analyzer instance!");
		return false;
	}

	if ( vxlgrid.nx >= UINT32MX || vxlgrid.ny >= UINT32MX || vxlgrid.nz >= UINT32MX || vxlgrid.nxyz >= UINT32MX ) {
		stopping( owner->owner->get_rank(), "Number of bins exceeds UINT32MX implement size_t HK!");
		return false;
	}

	//now cast safe
	unsigned int nx = static_cast<unsigned int>(vxlgrid.nx);
	unsigned int ny = static_cast<unsigned int>(vxlgrid.ny);
	unsigned int nz = static_cast<unsigned int>(vxlgrid.nz);
	unsigned int nxy = static_cast<unsigned int>(vxlgrid.nxy);
	unsigned int nxyz = static_cast<unsigned int>(vxlgrid.nxyz);

	//pass binarized ion occupancy bin bitfield to the HK analyzer to perform clustering analysis as one would do for percolation studies

	if ( hk->initialize( IsVacuum, nx, ny, nz) == false ) {
		stopping( "HoshenKopelman initialization failed during binning attempt!" );
		delete hk; return false;
	}

	if ( hk->hoshen_kopelman() == false ) {
		stopping( "HoshenKopelman run failed during binning attempt!" );
		delete hk; return false;
	}

	if ( hk->compactify() == false ) {
		stopping( "HoskenKopelman label compactification failed during binning attempt!");
		delete hk; return false;
	}

	/*if ( hk->checkLabeling() == false ) {
		stopping( "HoshenKopelman found labeling inconsistency");
		delete hk; return false;
	}*/

	if ( hk->determine_clustersize_distr() == false ) {
		stopping( "HoshenKopelman cluster size distribution failed during binning attempt!");
		delete hk; return false;
	}

	//if not yet return everything worked out so far so
	//MK::find HK label of frontmost top left voxel which is guaranteed to be in the guardzone, therefore has no ion
	//to use it to reassign which bins are vacuum and which not
	unsigned int tfl_corner_b = 0 + 0*nx + (nz-1)*nxy;

	if ( hk->rebinarize( tfl_corner_b, nxyz, IsVacuum ) == false ) {
		stopping( "HoshenKopelman rebinarization attempt failed!");
		delete hk; return false;
		//upon exit bitmap1 is true for bins in vacuum and false for tip
	}

	reporting( "HoshenKopelman clustering analysis was successful");
	delete hk; return true;
}


void vxlizer::binarization()
{
	//initial condition is that IsVacuum[i] = false \forall i
	for(size_t b = 0; b < vxlgrid.nxyz; ++b)
		IsVacuum[b] = false;

	//sequentially identify which bins contain at all data
	size_t nt = owner->sp->db.size();
	threadmemory* thisregion = NULL;
	for(size_t mt = MASTER; mt < nt; mt++) {
		thisregion = owner->sp->db.at(mt);
		//size_t n = thisregion->ionpp3.size();
		for ( auto it = thisregion->ion2bin.begin(); it != thisregion->ion2bin.end(); ++it ) {
			if ( IsVacuum[*it] == true ) //visited already
				continue;
			else
				IsVacuum[*it] = true;
		}
	} //for all thread
}


bool vxlizer::allocate_binaries()
{
	//do not reallocate!
	if ( IsVacuum != NULL )
		return false;
	if ( IsSurface != NULL )
		return false;
	if ( IsInside != NULL )
		return false;

	try {
		IsVacuum = new bool[vxlgrid.nxyz];
		IsSurface = new bool[vxlgrid.nxyz];
		IsInside = new bool[vxlgrid.nxyz];
	}
	catch (bad_alloc &surfexc) {
		stopping("Allocation of binary containers failed during binning attempt");
		return false;
	}
	return true;
}


sqb vxlizer::define_vxlgrid( const apt_xyz bwidth )
{
	sqb res = sqb();

	//get dimensions of tip AABB
	aabb3d ti = owner->sp->tip;

	//expand by one voxel guardzone of edgelen bwidth
	apt_real guardwidth = static_cast<apt_real>(bwidth);
	ti.blowup( guardwidth );

	res.nx = static_cast<size_t>( ceil((ti.xmx - ti.xmi) / bwidth) ); res.nx += 2; //MK::one layer guardzone on EACH side
	res.ny = static_cast<size_t>( ceil((ti.ymx - ti.ymi) / bwidth) ); res.ny += 2;
	res.nz = static_cast<size_t>( ceil((ti.zmx - ti.zmi) / bwidth) ); res.nz += 2;

	res.nxy = res.nx * res.ny;
	res.nxyz = res.nxy * res.nz;

	res.width = bwidth;
	res.box = ti;

	return res;
}


void vxlizer::rectangular_binning()
{
	double tic = MPI_Wtime();

	//tipvolume AABB binned into voxel of physical edge length db surplus a one voxel guardzone at each side into vacuum
	aabb3d domain = owner->sp->tip;

	apt_xyz binwidth = 0.5*(Settings::AdvIonPruneBinWidthMin + Settings::AdvIonPruneBinWidthMax);
	vxlgrid = define_vxlgrid( binwidth );

	string mess = "TipBinning with a nxyz grid " + to_string(vxlgrid.nx) + ";" + to_string(vxlgrid.ny) + ";" + to_string(vxlgrid.nz) + " binwidth " + to_string(binwidth) + " nm";
	reporting( mess );

	//allocate bitmap to speed up cache local testing of which bins contain ions and which any
	if ( allocate_binaries() == false ) {
		stopping( "Unable to allocate memory for binning ion cloud"); return;
	}

	//thread-local binning
	#pragma omp parallel
	{
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());
		owner->sp->db.at(mt)->binning( vxlgrid );
	}
	reporting( owner->owner->get_rank(), "Rectangular binning ion point cloud was binned successfully" );

	binarization();

	if ( identify_vacuum_bins() == false ) {
		stopping( "Unable to perform HoshenKopelman clustering analysis while binning ion cloud"); return;
	}

	identify_surface_adjacent_bins();

	identify_inside_bins();

	characterize_binning();

	//MK::we do not deallocate the binning memory!

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "TipBinning", APT_UTL, APT_IS_PAR, mm, tic, toc);
	cout << "TipBinning completed took " << (toc-tic) << " seconds" << endl;
}


clusterer::clusterer()
{
	owner = NULL;
	healthy = true;
}

clusterer::~clusterer()
{
	//MK::do not delete owner only backreference!
}


void clusterer::initialize_clustering_tasks()
{
	double tic = MPI_Wtime();

	parse_tasks( Settings::ClusteringTasksCode, clustering_tasks, owner->owner->mypse );

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ClusteringInitTasks", APT_UTL, APT_IS_SEQ, mm, tic, toc);
}


void clusterer::maximum_separation_method()
{
	//performs for each task, parameter sweeping of maximum separation method for dmax within range [ClustMSDmaxMin,ClustMSDmaxIncr,ClustMSDmaxMax]
	//without dilatation and erosion and setting Nmin to ClustMSNmin

	//distinguish time spent in clustering itself and a posteriori spat stat
	double tic = MPI_Wtime();

	//##MK::clustering tasks are executed sequentially having always all threads in parallel working on each
	cout << "Threadparallel maximum separation method using Goetz et al. HPDBScan..." << endl;

	//final result of analysis is diagram number of clusters = f(dmax) and cluster size distribution
	//##MK::are clusters with boundary contact excluded from the analysis? i.e. unbiased cluster size distribution

	for(size_t tsk = 0; tsk < clustering_tasks.size(); ++tsk) {

		//build a new class instance which performs a task-specific clustering analysis
		clustertask* cworker = NULL;
		try {
			cworker = new clustertask;
		}
		catch (std::bad_alloc &clcroak) {
			stopping("Unable to allocate memory to perform cluster analysis");
			return;
		}

		cworker->boss = this;
		cworker->mission = clustering_tasks.at(tsk);
		cworker->tskid = tsk;

		//prepare memory to store results
		vector<dbscanres> tskresults;

		//filter subset of ions, ##MK::benchmark performance, potentially improvements required
		cworker->initialize();

		if ( cworker->ions_filtered.size() > 0 ) {
			//perform task-specific HPDBScan clustering analysis with various dmax at fixed Nmin
			unsigned int dmaxid = 0;

			for( apt_real dmax = Settings::ClustMSDmaxMin; dmax <= Settings::ClustMSDmaxMax; dmax += Settings::ClustMSDmaxIncr, dmaxid++ ) {

				double mtic = omp_get_wtime();

				tskresults.push_back( cworker->hpdbscan( dmax, Settings::ClustMSNmin, dmaxid ) );
				//MK::all mark labels in the task-specific ions_filtered vector have now changed for non-noise ions
				//perform post-clustering analysis descriptive spatial statistics

				double mtoc = omp_get_wtime();
				cout << "\t\tTask/Dmax/NClusterFound took " << tsk << "\t\t" << dmax << " took " << (mtoc-mtic) << " seconds" << endl;

				if ( Settings::ClustPostSpatStat == true ) {
					double dtic = omp_get_wtime();

					if ( cworker->build_kdtree() == true ) {
						cworker->flag2exclude_guard( dmax );
						//task- and dmax-specific descriptive spatial statistics
						cworker->generic_spatstat( dmaxid );
						cworker->reset(); //MK::unflagging and resetting to reutilize parts of the datastructure in next dmax iteration
					}
					else {
						complaining("Unable to build task-specific KDTree for reduced ion cloud descriptive statistics!");
						continue;
					}

					double dtoc = omp_get_wtime();
					cout << "\t\tTask-specific spatial statistics done took " << (dtoc-dtic) << " seconds" << endl;
				} //done a posteriori spat stat
			} //probe next dmax of still the same task

			//Dmax parameter sweeping for task tsk done, report
			string whichtarg = ""; //##MK::deprecated, clustering_tasks.at(tsk).target.first;
			for( auto jt = clustering_tasks.at(tsk).trgcandidates.begin(); jt != clustering_tasks.at(tsk).trgcandidates.end(); jt++)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = clustering_tasks.at(tsk).envcandidates.begin(); jt != clustering_tasks.at(tsk).envcandidates.end(); jt++)
				whichcand += jt->first;

			if ( tskresults.size() > 0) { //task-specific results exist
				maximum_separation_report( whichtarg, whichcand, tskresults );
			}

			if ( cworker != NULL ) { //release tsk-specific resources
				delete cworker; cworker = NULL;
			}
		}
		else {
			string mess = "Processing maximum separation task " + to_string(tsk) + " was not processable!";
			complaining( mess); continue;
		}
	} //proceed with next task

	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ClusteringMaximumSeparationDmaxStudy", APT_CLU, APT_IS_PAR, mm, tic, toc);
	cout << "Conducting thread-parallelized HPDBScan-power Maximum Separation Method completed took " << (toc-tic) << " seconds" << endl;
}


void clusterer::maximum_separation_report( const string whichtarget, const string againstwhich,
		vector<dbscanres> const & results )
{
	double tic, toc;
	tic = MPI_Wtime();

//report geometry and spatial partitioning of the KDTree
	//are we randomized or not?
	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".MaxSepDmaxScan."
			+ whichtarget + "." + againstwhich + ".csv";

	ofstream sslog;
	sslog.open( fn.c_str() );
	sslog.precision(18);

	//##MK::as maximum separation method is not DBScan the output of DBScan itself is useless the cluster population is filtered afterwards
	//##MK::check however whether the HPDBScan implementation checks at all whether the neighborhood query function returns neighbors excluding or including the query point

	sslog << "DMaxID;DMax;NMin;Ncluster\n"; //;Nlinkpoints;Ncorepoints;Nnoisepoints;Nclusteredpoints\n";
	sslog << "1;nm;1;1\n"; //;1;1;1;1\n";
	sslog << "DMaxID;DMax;NMin;Ncluster\n"; //;Nlinkpoints;Ncorepoints;Nnoisepoints;Nclusteredpoints\n";
	size_t id = 0;
	for( auto it = results.begin(); it != results.end(); ++it, id++ ) {
		sslog << id << ";" << it->Dmax << ";" << it->Nmin << ";" << it->nClusterFound << "\n";
		//sslog << ";" << (it->nClustered - it->nCore) << ";" << it->nCore << ";" << it->nNoise << ";" << it->nClustered << "\n";
	}
	sslog.flush();
	sslog.close();

	toc = MPI_Wtime();
	cout << "Reporting MaximumSeparation dmax range scan results took " << (toc-tic) << " seconds" << endl;
}


vics_materialpoint_result::vics_materialpoint_result()
{
	MatPointPos = p3d();
	FFTSummary = vicsfftsummary();
}


vics_materialpoint_result::vics_materialpoint_result( const p3d here )
{
	MatPointPos = here;
	FFTSummary = vicsfftsummary();
}


vics_materialpoint_result::~vics_materialpoint_result()
{
}



crystindexer::crystindexer()
{
	owner = NULL;
	info = vicsmeta();
}


crystindexer::~crystindexer()
{
	//do not delete owner is only a backreference
}


void crystindexer::configure()
{
	//configure threadlocal worker and initiate FFT plans
	//for every point we need to scan through elevation and azimuth space
	int NumberOfBins = pow( static_cast<int>(2), static_cast<int>(Settings::CrystalloHistoM) );
	apt_real R = Settings::CrystalloRadiusMax;
	apt_real dR = 2.f*R / (static_cast<apt_real>(NumberOfBins) - 2.f);
	apt_real binner = static_cast<apt_real>(NumberOfBins) / (2.f*R+2.f*dR);
	info = vicsmeta( R, dR, binner, NumberOfBins );

	#pragma omp critical
	{
		cout << "Thread/NumberOfBins/R/dR/Binner = " << omp_get_thread_num() << ";" << info.NumberOfBins << ";" << info.R << ";" << info.dR << ";" << info.binner << endl;
	}


	//precompute windowing coefficients
	//J. F. Kaiser, R. W. Schafer, 1980, On the use of the I_0-sinh window for spectrum analysis,
	//doi::10.1109/TASSP.1980.1163349
	window_coeff = vector<float>( NumberOfBins, 1.f );
	//window_coeff = vector<double>( NumberOfBins, 1.f );
	if ( Settings::WindowingMethod == E_KAISER_WINDOW ) {
		double alpha = Settings::WindowingAlpha;
		double zero = 0.f;
		double I0a = boost::math::cyl_bessel_i( zero, alpha );

		for ( int i = 0; i < NumberOfBins; i++ ) {
			//w(n) modified Bessel function of first kind with shape parameter WindowingAlpha
			//w(n) = I_0(alpha*sqrt(1-(n-(N/2)/(N/2)))^2)) / I_0(alpha) for 0<= n <= N-1  otherwise w(n) = 0

			double nN = (static_cast<double>(i) - (static_cast<double>(NumberOfBins) / 2.f)) / (static_cast<double>(NumberOfBins) / 2.f);
			double x = alpha * sqrt( 1.f - SQR(nN) );
			double I0x = boost::math::cyl_bessel_i( zero, x );

			window_coeff.at(i) = static_cast<float>(I0x / I0a); //high precision bessel then cut precision
cout << setprecision(32) << "Computing modified Bessel function of first kind I0(v,x) coefficients " << i << "\t\t" << window_coeff.at(i) << "\n";
		}
cout << "Kaiser windowing defined" << endl;
	}
	else { //E_RECTANGULAR_WINDOW
		window_coeff.at(0) = 0.f;
		window_coeff.back() = 0.f;
cout << "Rectangular windowing defined" << endl;
	}
}


/*
vicsresult crystindexer::ExecutePeaksFinding1D( vector<pair<double,double>> & CntsVsFreq )
{
	//##MK::implement 1D peak finding algorithm
	//MK::quick and dirty here O(N), divide and conquer O(lg(N)) solution should be possible for future...

	std::sort( CntsVsFreq.begin(), CntsVsFreq.end() );

	vicsresult tmp = vicsresult();


//	singleresult.fft_allc_success = true;
//	singleresult.fft_plac_success = true;
//	singleresult.fft_init_success = true;
//	singleresult.fft_comp_success = true;
//	singleresult.fft_free_success = true;
	tmp.max1 = CntsVsFreq.at(CntsVsFreq.size()-1);
	tmp.max2 = CntsVsFreq.at(CntsVsFreq.size()-2);
	tmp.max3 = CntsVsFreq.at(CntsVsFreq.size()-3);

	return tmp;
}


vicsfftstatus crystindexer::ExecuteSingleDiscreteFFT(
		vector<unsigned int> const & in, vector<pair<double,double>> & out )
{
	vicsfftstatus whatsapp = vicsfftstatus();

	////vector<float> dhistogram = vector<float>(info.NumberOfBins, 0.f);
	size_t NN = info.NumberOfBins;
	vector<float> dhistogram( NN, 0.f);

	//windowing of the histogram and normalizing
	//float norm = static_cast<float>(cand.size());
	for( size_t i = 1; i < NN-1; i++ ) { //+1 and -1 because left and right pad position
		////dhistogram[i] = window_coeff[i] * static_cast<float>(uihistogram[i]); // / norm;
		//float cnts = in[i];
		dhistogram[i] = window_coeff[i] * in[i]; //cnts; // / norm;
//##cout << i << "\t\t" << dhistogram[i] << endl;
	}

	rfftn* anfft = NULL;
	try {
		anfft = new rfftn;
	}
	catch (std::bad_alloc &ompcroak)
	{
		return whatsapp;
	}

	//cout << "Attempting1 an FFT anfft __" << anfft << "___" << endl;
	anfft->init( NN );
	//cout << "Attempting2 an FFT anfft __" << anfft << "___" << endl;
	anfft->fill( dhistogram );
	//cout << "Attempting3 an FFT anfft __" << anfft << "___" << endl;
	anfft->forwardFFT();
	//cout << "Done with an FFT anfft __" << anfft << "___" << endl;

	//create plan for executing a 1D discrete fast Fourier transform using the IntelMKL library
	//##MK::in the future can reutilize this plan FFT, histograms have the same size for every case
	vector<char>* fftTemp = new vector<char>;
//	MKL_LONG NI = 512; //##MK::
//	MKL_LONG requiredSize = NI * sizeof(MKL_Complex16);
//	fftTemp->resize(requiredSize);
//	MKL_Complex16* fftTempP = (MKL_Complex16*) &(*fftTemp)[0];



//	DFTI_CONFIG_VALUE precision = DFTI_DOUBLE; //DFTI_SINGLE;
//	DFTI_DESCRIPTOR_HANDLE descr;
//	MKL_LONG dimensions = dhistogram.size();

//	//perform single precision discrete Fourier transform of this histogram using single threaded Intel MKL
//	////vector<std::complex<float>> out( dhistogram.size() );
//	vector<std::complex<double>> tmp( dimensions, std::complex<double>(0.f,0.f) );

//	MKL_LONG status;
//	status = DftiCreateDescriptor(&descr, precision, DFTI_COMPLEX, 1, dimensions);
//	if ( status == DFTI_NO_ERROR )
//		whatsapp.fft_allc_success = true;

//	status = DftiSetValue(descr, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
//	if ( status == DFTI_NO_ERROR )
//		whatsapp.fft_plac_success = true;

//	status = DftiCommitDescriptor(descr);
//	if ( status == DFTI_NO_ERROR )
//		whatsapp.fft_init_success = true;

//	status = DftiComputeForward(descr, dhistogram.data(), tmp.data() ); //&dhistogram[0], fftTempP );the FFT computation itself
//	if ( status == DFTI_NO_ERROR )
//		whatsapp.fft_comp_success = true;

//	status = DftiFreeDescriptor(&descr);
//	if ( status == DFTI_NO_ERROR )
//		whatsapp.fft_free_success = true;

//	//compute absolute value of Fourier transform results in first part of conjugate-even halfspace
//	////size_t NN = static_cast<size_t>(ceil( 0.5*static_cast<float>(out.size())));
//	//size_t NN = 256; //######static_cast<size_t>(ceil( 0.5*static_cast<double>(out.size()) ));
//	////vector<pair<apt_real, apt_real>> ea_result( NN, make_pair(0.f, 0.f) );
//	//vector<pair<double, double>> ea_result( NN, make_pair(0.f, 0.f) );
//	////float mult = 1.f / static_cast<float>(NN) * (info.R + info.dR); //nm ##MK
//	size_t NN = out.size();


	float mult = static_cast<float>(info.R + info.dR) / static_cast<float>(NN);
//##cout << "FFT result component values" << endl;
	size_t ni = out.size();
	for ( size_t i = 0; i < ni; ++i ) {
		////float reci_pos = static_cast<float>(i) * mult;
		////float fft_abs = sqrt( SQR(out[i].real()) + SQR(out[i].imag()) );
		double reci_pos = static_cast<double>(i) * mult;
		double fft_abs = sqrt( SQR(anfft->fftTempP[i].real) + SQR(anfft->fftTempP[i].imag) );
		//double fft_abs = sqrt( SQR(fftTempP[i].real) + SQR(fftTempP[i].imag) );

		//if ( isnan(fft_abs) )
		//	cout << i << "\t\t" << tmp[i].real() << "\t\t" << tmp[i].imag() << "\t\t" << reci_pos << ";" << fft_abs << endl;

		out[i] = make_pair( reci_pos, fft_abs ); //-1.f*fft_abs );
	}

	//free resources
	delete anfft;

	//cout << "FreedAnFFT" << endl;

	//we are interest in finding those pairs with the n largest signal magnitude
	//if working naively were could achieve this by dumping magnitude values in our ea_result buffer
	//sort the buffer ascendingly and pick the largest n i.e. last n values
	//however this not efficient
	//better we multiply the magnitude with -1, given that its real and >= 0 this factually
	//reverse the sorting, for instance results like 1.2, 4.5, 3.1, 7.8, 5.1 will be ordered
	//like so -7.8, -5.1, -4.5, -3.1, -1.2 now we just need to pick the first n objects using nth_element
	//with usually O(N) instead of O(NlgN) complexity
	//and reverse sign (-1.f) only for the n such identified candidates

	//ea_result.at(j) = make_pair( static_cast<apt_real>(j), -1.f * sqrt(SQR(it->real())+SQR(it->imag())) );

	//std::complex<float> PconjP = *it * std::conj(*it);
	//ea_result.at(j) = make_pair( static_cast<apt_real>(j), -1.f*PconjP.real() );

	////std::complex<float> P = *it;
	////std::complex<float> PconjP = P * std::conj(*it);
	////Pyy.push_back( P * std::conj(*it) );
	///cout << P.real() << "\t\t" << P.imag() << "\t\t" << Pyy.real() << "\t\t" << Pyy.imag() << endl;


	return whatsapp;
}


vicsresult crystindexer::SingleElevationAzimuth( const apt_real ee, const apt_real aa,
		vector<d3d> const & ions, vicsfftsummary & diary )
{
	//reset temporaries
	vicsresult singleresult = vicsresult();

	vector<unsigned int> uihistogram( info.NumberOfBins, 0 );
	size_t NN = static_cast<size_t>(ceil( 0.5*static_cast<float>(info.NumberOfBins) ));

	vector<pair<double,double>> ea_result( NN, make_pair(0.f,0.f) );
	//vector<double> ea_result(2*NN, 0.f);
	//compute signed distance of every candidate point to current plane given current potential crystal plane normal vector
	//d = dot(n,cand), //##MK::employ bSIMD here make cand aligned memory
	for ( auto kt = ions.begin(); kt != ions.end(); kt++ ) {
		apt_real d = 	+ sin(ee) * cos(aa) * kt->u
						- sin(ee) * sin(aa) * kt->v
						+           cos(ee) * kt->w;

		apt_real b = floor((d + info.R + info.dR)*info.binner);
		unsigned int bin = static_cast<unsigned int>(b);
		uihistogram[bin]++;
//cout << d << "\t\t" << b << "\t\t" << bin << "\t\t" << uihistogram[bin] << endl;
	}

	//check if histogram front and end are really zero
	if ( uihistogram.at(0) == 0 && uihistogram.back() == 0 ) {
//cout << e << "\t\t" << a << "\t\t-->doing\n";

		vicsfftstatus stats = ExecuteSingleDiscreteFFT( uihistogram, ea_result );

		singleresult = ExecutePeaksFinding1D( ea_result );
		singleresult.elevation = ee;
		singleresult.azimuth = aa;

		//if ( stats.fft_allc_success == true && stats.fft_plac_success == true &&
		//		stats.fft_init_success == true && stats.fft_comp_success == true &&
		//		stats.fft_free_success == true ) {

	//		diary.nFFTsSuccess++;

	//		//post-process result
	//		std::sort( ea_result.begin(), ea_result.end() );


	//		singleresult.fft_allc_success = true;
	//		singleresult.fft_plac_success = true;
	//		singleresult.fft_init_success = true;
	//		singleresult.fft_comp_success = true;
	//		singleresult.fft_free_success = true;
	//		singleresult.max1 = ea_result.at(ea_result.size()-1);
	//		singleresult.max2 = ea_result.at(ea_result.size()-2);
	//		singleresult.max3 = ea_result.at(ea_result.size()-3);
	//	}
	//	else {
	//		diary.nFFTsFailed++;
	//	}
	//
	}
	//else {
	//	diary.nFFTsFailed++;
	//}

	return singleresult;


//		//pull entire histogram
//		//##MK::comment out to avoid storing excessive intermediate results
//		for( size_t i = 0; i < NN; ++i ) {
//			dumpfull[i].cnts = dhistogram[i];
//			dumpfull[i].FFTr = out[i].real();
//			dumpfull[i].FFTi = out[i].imag();
//			dumpfull[i].pw = sqrt( SQR(out[i].real()) + SQR(out[i].imag()) );
//			//std::complex<float> PconjP = out[i] * std::conj(out[i]);
//			//dumpfull[i].pw = PconjP.real();
////##cout << dumpfull[i].cnts << ";" << dumpfull[i].FFTr << ";" << dumpfull[i].FFTi << ";" << dumpfull[i].pw << endl;
//		}

		//cout << e << "\t\t" << a << "\t\t-->success\n";cout << e << "\t\t" << a << "\t\t-->SUCCESS hISTOGRAM padding region contains non-zero values!\n";
}


void crystindexer::vicsmethod_core_incrementalfft( p3d const & target, vector<d3d> const & cand )
{
	//scan through elevation and azimuth, for each pair project candidate
	//MK::the reason to use only power of 2 number of bins is that FFT is most efficient for such cases

	//add a dR guard zone about the histogram to the left and right ends
	//this effectively pads the histogram to zero, alternatively: ##MK::use windowing approaches like Kaiser etc...
	vector<vicsresult> & dumpsummary = res.back().ElevAzimTblSummary;
	//size_t NN = static_cast<size_t>(ceil( 0.5*static_cast<float>(info.NumberOfBins) ));

	//##MK::comment out to avoid storing excessive intermediate results
	//res.back().ElevAzimTblHistogram.push_back( vector<vicshistbin>( NN, vicshistbin()) );
	//vector<vicshistbin> & dumpfull = res.back().ElevAzimTblHistogram.back();

	//temporaries
	//vector<unsigned int> uihistogram( info.NumberOfBins, 0);
	//vector<pair<double,double>> ea_result( NN, make_pair(0.f, 0.f) );
	//vicsresult health = vicsresult();
	vicsfftsummary diary = vicsfftsummary();
	diary.nIons = cand.size();

	//elevation/azimuth loop nest
	for( apt_real e = Settings::ElevationAngleMin; e <= Settings::ElevationAngleMax; e += Settings::ElevationAngleIncr ) {
		//apt_real e = DEGREE2RADIANT(-5.f);
		//apt_real a = DEGREE2RADIANT(0.f);

		//apt_real sin_e = sin(e);
		//apt_real cos_e = cos(e);
		for ( apt_real a = Settings::AzimuthAngleMin; a <= Settings::AzimuthAngleMax; a += Settings::AzimuthAngleIncr ) {
			//apt_real sin_a = sin(a);
			//apt_real cos_a = cos(a);

			vicsresult thisone = SingleElevationAzimuth( e, a, cand, diary );

			dumpsummary.push_back( thisone );
		} //next azimuth value to analyze at fixed elevation
	} //next set of azimuth values for a given elevation

	res.back().FFTSummary = diary;
}
*/


void crystindexer::vicsmethod_core_incrementalfft2( p3d const & target, vector<d3d> const & cand )
{
	//MK::the reason to use only power of 2 number of bins is that FFT relatively most efficient
	vector<vicsresult> & dumpsummary = res.back().ElevAzimTblSummary;
	bool ReportHeavyData = Settings::IOCrystallography;

	//MK::comment out to avoid storing excessive intermediate results
	//res.back().ElevAzimTblHistogram.push_back( vector<vicshistbin>( NN, vicshistbin()) );
	//vector<vicshistbin> & dumpfull = res.back().ElevAzimTblHistogram.back();

	vicsfftsummary diary = vicsfftsummary();
	diary.threadID = omp_get_thread_num();
	diary.nIons = cand.size();

	//elevation/azimuth loop nest
	size_t NumberOfBins = info.NumberOfBins;
	size_t NumberOfBinsHalf = static_cast<size_t>(ceil(0.5*static_cast<double>(NumberOfBins)));
	float RdR = info.R + info.dR;
	float BinMultiplier = info.binner;

	//mapping of FFT bins to reciprocal frequency bins
	vector<float> frequencies( NumberOfBinsHalf, 0.f);
	float NormalizerFreq = RdR / static_cast<float>(NumberOfBinsHalf);
	for( size_t i = 0; i < NumberOfBinsHalf; ++i) {
		frequencies[i] = static_cast<float>(i) * NormalizerFreq;
	}

	//elevation/azimuth nested main loop
	for( apt_real e = Settings::ElevationAngleMin; e <= Settings::ElevationAngleMax; e += Settings::ElevationAngleIncr ) {
		//float e = DEGREE2RADIANT(-5.f); 	float a = DEGREE2RADIANT(0.f);
		float sin_e = sin(e);
		float cos_e = cos(e);
		for ( apt_real a = Settings::AzimuthAngleMin; a <= Settings::AzimuthAngleMax; a += Settings::AzimuthAngleIncr ) {
			float sin_a = sin(a);
			float cos_a = cos(a);
			float se_ca = +1.f*sin_e*cos_a;
			float se_sa = -1.f*sin_e*sin_a;
			float ce = +1.f*cos_e;

			//histogram of projected signed distance of every point to current plane
			//given the current plane normal vector parameterized through elevation and azimuth
			vicsresult thisone = vicsresult();
			thisone.elevation = e;
			thisone.azimuth = a;
			vector<unsigned int> uihistogram( NumberOfBins, 0 );
			for ( auto kt = cand.begin(); kt != cand.end(); ++kt ) { //d = dot(n,cand), ##MK::bSIMD
				float d = se_ca*kt->u + se_sa*kt->v + ce*kt->w;
				float b = floor((d + RdR)*BinMultiplier);
				unsigned int bin = static_cast<unsigned int>(b);
				uihistogram[bin]++;
				//cout << d << "\t\t" << b << "\t\t" << bin << "\t\t" << uihistogram[bin] << endl;
			}
			//vector<uint32_t, bs::allocator<uint32_t> uihistogram(NumberOfBins, 0); //##MK::check if filled with 0

			//did padding bins remain empty, if so we are ready for doing a discrete 1d FFT of dhistogram
			if ( uihistogram[0] == 0 && uihistogram[NumberOfBins-1] == 0 ) {

				//windowing of the histogram
				//potential normalizing, definately translate uint2f32 to avoid accumulation of single precision errors
				vector<float> dhistogram( NumberOfBins, 0.f);
				//float NormalizingMultiplier = 1.f / static_cast<float>(cand.size());
				for( size_t i = 1; i < NumberOfBins-1; ++i ) { //+1 and -1 because left and right pad position are zero already
					dhistogram[i] = window_coeff[i] * static_cast<float>(uihistogram[i]);
				}

				//create plan for executing a 1D discrete fast Fourier transform using the IntelMKL library
				rfftn* anfft = NULL;
				anfft = new rfftn;
				anfft->init( NumberOfBins, dhistogram );
				anfft->forwardFFT(); //execute this plan 1d discrete FFT using IntelMKL library
				//real to complex discrete 1d FFT has conjugate-even symmetry IMKL so work only on [0,NumberOfBinsHalf)
				anfft->getMagnitude( NumberOfBinsHalf, dhistogram );
				//MK::getMagnitude overwrites [0,NumberOfBinsHalf) with sqrt(zR^2+zI^2)
				delete anfft;

				//peak finding on [0,NumberOfBinsHalf), NumberOfBinsHalf is guaranteed = pow(2, (Settings::CrystalloHistoM-1) )
				//##MK::first shot, O(n) finding of all clear localmaxima, first/last value checked separately to avoid if inside for loop,
				vector<localmaximum> peaks;
				//first bin a peak?
				if ( dhistogram[0] > dhistogram[1] )
					peaks.push_back( localmaximum(frequencies[0], dhistogram[0]) ); //frequently the case 0-freq peak
				//peaks somewhere in the middle, ##MK::divide and conquer for better than O(n) complexity
				for ( size_t b = 1; b < NumberOfBinsHalf-1; ++b ) {
					if ( dhistogram[b-1] >= dhistogram[b] )
						continue;
					if ( dhistogram[b] <= dhistogram[b+1] )
						continue;
					//not continued
					peaks.push_back( localmaximum(frequencies[b], dhistogram[b]) );
				}
				//last bin a peak?
				if ( dhistogram[NumberOfBinsHalf-2] < dhistogram[NumberOfBinsHalf-1] )
					peaks.push_back( localmaximum(frequencies[NumberOfBinsHalf-1], dhistogram[NumberOfBinsHalf-1]) );

				//lastly find the three strongest peaks
				//MK::a naive full sort on peaks has O(NlgN) time complexity, consider multiply dhistogram with -1.f and do nth_element partial sort
				std::sort( peaks.begin(), peaks.end(), SortLocalMaximaForStrength );

				//evaluate 3 strongest peaks of current transform and report as result
				size_t idx = 0;
				if ( peaks.size() >= 3) {
					idx = peaks.size()-1;
					thisone.max1 = make_pair(peaks[idx].position, peaks[idx].strength); //1 strongest
					idx = peaks.size()-2;
					thisone.max2 = make_pair(peaks[idx].position, peaks[idx].strength);
					idx = peaks.size()-3;
					thisone.max3 = make_pair(peaks[idx].position, peaks[idx].strength);
				}

				//bookkeeping
				diary.nFFTsSuccess++;
			}
			else {
				diary.nFFTsFailed++;
			}

			if ( ReportHeavyData == true ) {
				dumpsummary.push_back( thisone );
			}

		} //next azimuth value to analyze at fixed elevation
	} //next set of azimuth values for a given elevation

	res.back().FFTSummary = diary;
}


/*void crystindexer::vicsmethod_core_batchfft_simd( p3d const & target, vector<d3d> const & cand )
{
	//MK::the reason to use only power of 2 number of bins is that FFT relatively most efficient
	vector<vicsresult> & dumpsummary = res.back().ElevAzimTblSummary;

	//MK::comment out to avoid storing excessive intermediate results
	//res.back().ElevAzimTblHistogram.push_back( vector<vicshistbin>( NN, vicshistbin()) );
	//vector<vicshistbin> & dumpfull = res.back().ElevAzimTblHistogram.back();

	vicsfftsummary diary = vicsfftsummary();
	diary.nIons = cand.size();


	//reduce plan construction, memory initialization, and other overhead by performing
	//a batch FFT, i.e. have a matrix of histograms and do single-threaded FFT of this at once
	//each row is a histogram, matrix has as many columns info.NumberOfBins
	//at higher memory initialization costs

	//evaluate how many elevation/azimuth pairs to process per material volume point
	size_t NumberOfTransforms = 0; //rows
	for( apt_real e = Settings::ElevationAngleMin; e <= Settings::ElevationAngleMax; e += Settings::ElevationAngleIncr )
		for ( apt_real a = Settings::AzimuthAngleMin; a <= Settings::AzimuthAngleMax; a += Settings::AzimuthAngleIncr )
			NumberOfTransforms++;
	size_t NumberOfBins = info.NumberOfBins; //columns

	//allocate a buffer for storing elevation/azimuth pair results
	vector<vicsresult> results( NumberOfTransforms, vicsresult());

	//evaluate memory costs
	size_t MemoryCosts = NumberOfTransforms * NumberOfBins * (sizeof(float) + 2*sizeof(float)); //initially needs uint32+float frees uint32 add 2xfloat
	if ( MemoryCosts > (1024*1024*1024)) { //1GB per thread
		cout << "Falling back to incrementally FFT" << endl;
		return; //##MK
	}

	//vectorized filling of buffer
	//using pack_f32 = bs::pack<float>; //we need aligned memory for vectorization to be maximally efficient
	//, bs::allocator<unsigned int>
	vector<unsigned int> cnts( NumberOfTransforms*NumberOfBins, 0);
	float doffset = info.R + info.dR;
	float binmult = info.binner;
	size_t TransformID = 0;
	for( apt_real e = Settings::ElevationAngleMin; e <= Settings::ElevationAngleMax; e += Settings::ElevationAngleIncr ) {
		//pack_f32 sin_e{sin(e)};
		//pack_f32 cos_e{cos(e)};
		float sin_e = sin(e);
		float cos_e = cos(e);
		for ( apt_real a = Settings::AzimuthAngleMin; a <= Settings::AzimuthAngleMax; a += Settings::AzimuthAngleIncr ) {
			//pack_f32 sin_a{sin(a)};
			//pack_f32 cos_a{cos(a)};
			float sin_a = sin(a);
			float cos_a = cos(a);
			float doffset = info.R + info.dR;
			size_t TransformOffset = TransformID * NumberOfBins;
			float sin_e__cos_a = +1.f * sin_e * cos_a;
			float sin_e__sin_a = -1.f * sin_e * sin_a;
			for( auto kt = cand.begin(); kt != cand.end(); ++kt ) { //cand is contiguous in cache
				//dot product distance to projected on rotated inclined hkl plane
				float d = sin_e__cos_a*kt->u  + sin_e__sin_a*kt->v  + cos_e*kt->w;
				float b = floor((d + doffset)*binmult);
				size_t bin = TransformOffset + static_cast<size_t>(b);
				cnts[bin]++; //writing to restricted interval of contiguous cache
			}
			results[TransformID].elevation = e;
			results[TransformID].azimuth = a;
			TransformID++;
		}
	}

	//windowing
	vector<float> intensity( NumberOfTransforms*NumberOfBins, 0.f);
	for( size_t trsfrm = 0; trsfrm < NumberOfTransforms; ++trsfrm ) { //contiguous cache again
		size_t PositionOffset = trsfrm * NumberOfBins;
		if ( cnts.at(PositionOffset+0) == 0 && cnts.at(PositionOffset+NumberOfBins-1) == 0 ) {
			//add a dR guard zone about the histogram to the left and right ends
			//this effectively pads the histogram to zero, in addition use windowing
			size_t pos0 = PositionOffset+1;
			size_t pos1 = PositionOffset+NumberOfBins-1;
			for( size_t i = 1; i < NumberOfBins-1; ++i ) { //+1 and -1 because left and right pad position
				intensity[PositionOffset+i] = window_coeff[i] * static_cast<float>(cnts[PositionOffset+bin]);
			}
		}
		else {
			#pragma omp critical
			{
				cerr << "Thread " << omp_get_thread_num() << " point " << target.x << ";" << target.y << ";" << target.z << " the padding region of transform " << trsfrm << " is has non-zero counts!" << endl; return;
			}
		}
	}

	//at this point cnts is no longer necessary so delete it
	cnts.clear();

	//batch discrete 1d FFT using IMKL
	batch1dfft* bfft = NULL;
	try {
		bfft = new batch1dfft;
	}
	catch (std::bad_alloc &ompcroak)
	{
		cerr << "Thread " << omp_get_thread_num() << " unable to allocate BatchFFT processor class object" << endl;
		return;
	}

	bfft->init( NumberOfBins );
	bfft->fill( intensity );
	bfft->comp_fFFT();

	//now we are at peak memory consumption so take a snapshot here
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();

	//pull results inplace back into intensity overwrite intensity inplace
	bfft->get_abs( intensity );

//	///get FrequencySpace positions
//	size_t NumberOfBinsHalf = static_cast<size_t>(ceil(0.5*static_cast<float>(NumberOfBins)));
//	for( size_t trsfrm = 0; trsfrm < NumberOfTransforms; ++trsfrm ) {
//		size_t
//		for ( size_t i = 0; i < ni; ++i ) {
//			intensity[i] = sqrt( SQR(bfft->fftTempP[i].real) + SQR(anfft->fftTempP[i].imag) );
//	}


	//and get corresponding
	size_t NumberOfBinsHalf = bfft->get_nihalf();
	float ReciprocFreqSpaceScaling = static_cast<float>(info.R + info.dR) / static_cast<float>(NumberOfBinsHalf);
	vector<float> reci_pos( NumberOfBinsHalf, 0.f );
	for( size_t b = 0; b < NumberOfBinsHalf; ++b ) {
		reci_pos[b] = static_cast<float>(b) * ReciprocFreqSpaceScaling;
	}

	//we have pull the FFT so can delete bfft
	delete bfft;

	//finally perform peak finding for each histogram
	vector<localmaximum> peaks;
	for( size_t trsfrm = 0; trsfrm < NumberOfTransforms; trsfrm++ ) {
		size_t PositionOffset = trsfrm * NumberOfBins;
		//real to complex discrete 1d FFT has conjugate-even symmetry IMKL so work only on [0,NumberOfBinsHalf)

		//for every elevation/angle pair NumberOfBinsHalf is guaranteed = pow(2, (Settings::CrystalloHistoM-1) )
		peaks.clear();

		//##MK::to begin with stable O(n) time complex finding of all localmaxima
		//check first value and last separately to avoid if in for loop,
		//work on float intensity (which is now FFT value abs) [PositionOffset,PositionOffset+NumberOfBinsHalf)

		//is first value a clear peak?
		if ( intensity[PositionOffset] > intensity[PositionOffset+1] )
			peaks.push_back( localmaximum( reci_pos[0], intensity[PositionOffset] ) );
		//somewhere in a middle clear peaks?
		//##MK::divide and conquer for better than O(n) complexity
		for ( size_t b = PositionOffset+1; b < PositionOffset+NumberOfBinsHalf-1; b++ ) {
			if ( intensity[b-1] < intensity[b] && intensity[b+1] < intensity[b] ) {
				peaks.push_back( localmaximum(reci_pos[b-PositionOffset], intensity[b]) );
			}
		}
		//is last value a clear peak?
		if ( intensity[PositionOffset+NumberOfBinsHalf-2] < intensity[PositionOffset+NumberOfBinsHalf-1] )
			peaks.push_back( localmaximum(reci_pos[NumberOfBinsHalf-1], intensity[PositionOffset+NumberOfBinsHalf-1]) );

		//##MK::find n-th strongest peaks
		//MK::a naive full sort on peaks has O(NlgN) time complexity
		std::sort( peaks.begin(), peaks.end(), SortLocalMaximaForStrength );

		//evaluate 3 strongest peaks of current transform and report as result
		size_t thisone = 0;
		if ( peaks.size() >= 3) {
			thisone = peaks.size()-1;
			results[trsfrm].max1 = make_pair(peaks.at(thisone).position, peaks.at(thisone).strength); //1 strongest
			thisone = peaks.size()-2;
			results[trsfrm].max2 = make_pair(peaks.at(thisone).position, peaks.at(thisone).strength); //2 strongest
			thisone = peaks.size()-3;
			results[trsfrm].max3 = make_pair(peaks.at(thisone).position, peaks.at(thisone).strength); //3 strongest
		}
		else {
			results[trsfrm].max1 = make_pair(F32MX,F32MX);
			results[trsfrm].max2 = make_pair(F32MX,F32MX);
			results[trsfrm].max3 = make_pair(F32MX,F32MX);
		}

		dumpsummary.push_back( results[trsfrm] );
	} //evaluate result of the next elevation / azimuth

	//##MK::for low values of n it might be potentially more efficient to take -1.f*strength and to nth_element

	res.back().FFTSummary = diary;
}
*/


aptcrystHdl::aptcrystHdl()
{
	owner = NULL;
	probehere = cuboidgrid3d();
}


aptcrystHdl::~aptcrystHdl()
{
	//do not delete owner, only backreference
}


void aptcrystHdl::compute_crystallography_definegrid()
{
	//define a grid of cuboidal cells at which vertices we conduct V. J. Araullo-Peters et al., 2015 technique to identify crystallographic data
	probehere.binwidths = p3d( 	Settings::SamplingGridBinWidthX,
								Settings::SamplingGridBinWidthY,
								Settings::SamplingGridBinWidthZ );
	probehere.ve = owner->sp->tip;
	probehere.ve.scale();
	p3d cxyz = probehere.ve.center();

	cout << "Tip center " << cxyz << endl;

	//MK::grid is defined symmetric about world coordinate z axis
	//##MK::dirty debugging
	probehere.gridcells = p6i(
			static_cast<int>(floor((probehere.ve.xmi-cxyz.x)/probehere.binwidths.x)),
			static_cast<int>(ceil((probehere.ve.xmx-cxyz.x)/probehere.binwidths.x)),
			static_cast<int>(floor((probehere.ve.ymi-cxyz.y)/probehere.binwidths.y)),
			static_cast<int>(ceil((probehere.ve.ymx-cxyz.y)/probehere.binwidths.y)),
			static_cast<int>(floor((probehere.ve.zmi-cxyz.z)/probehere.binwidths.z)),
			static_cast<int>(ceil((probehere.ve.zmx-cxyz.z)/probehere.binwidths.z))  );

	//##MK::nicify implementation of the grid definitions
	for( int z = probehere.gridcells.zmi; z <= probehere.gridcells.zmx; ++z) {
		for( int y = probehere.gridcells.ymi; y <= probehere.gridcells.ymx; ++y) {
			for( int x = probehere.gridcells.xmi; x <= probehere.gridcells.xmx; ++x) {
				p3d p = probehere.where( x, y, z);
				if ( p.x < (F32MX-EPSILON) && p.y < (F32MX-EPSILON) && p.z < (F32MX-EPSILON) ) {
					//cout << "-->Adding point " << p << endl;
					samplingpoints.push_back( p );
				}
			}
		}
	}

	cout << "APTCrystallography Vic's method grid defined" << endl;
	cout << probehere << endl;
	cout << "Number of samplingpoints " << samplingpoints.size() << endl;

	//define which elevation/azimuth pairs to sample for every point
	for( apt_real e = Settings::ElevationAngleMin; e <= Settings::ElevationAngleMax; e += Settings::ElevationAngleIncr )
		for ( apt_real a = Settings::AzimuthAngleMin; a <= Settings::AzimuthAngleMax; a += Settings::AzimuthAngleIncr )
			elevazimpairs.push_back( p2d(e,a) );

	cout << "APTCrystallography Vic's method elevation/azimuth pairs to test defined" << endl;
	cout << "Number of pairs to test per materialpoint " << elevazimpairs.size() << endl;
}


void aptcrystHdl::compute_crystallography_vicsmethod()
{
	double tic = MPI_Wtime();

	compute_crystallography_definegrid();

	//##MK::change so that threads work primarily on points in their own region
	#pragma omp parallel
	{
		double mytic = omp_get_wtime();

		int mt = omp_get_thread_num();
		int nt = omp_get_num_threads();

		#pragma omp master
		{
			for( int thr = MASTER; thr < nt; thr++ ) {
				workers.push_back( NULL);
			}
		}
		//necessary as otherwise write back conflict on workers from threads
		#pragma omp barrier

		//setup threadlocal worker which use threadlocal memory
		crystindexer* thrindexer = NULL;
		try {
			thrindexer = new crystindexer;
		}
		catch (bad_alloc &mecroak) {
			//##MK::#####
		}

		//no workers backreference holder can safely be modified by current thread
		//MK::pointer will point to thread-locally allocated memory
		#pragma omp critical
		{
			workers.at(mt) = thrindexer;
		}

		thrindexer->configure();

		apt_real R = thrindexer->info.R;
		apt_real RSQR = SQR(R);

		//get access to a trihull skip samplingpoints if there are not deep enough in tip volume
		//to get full spherical environment/ROI of radius R
		Tree* mylinkedbvh = owner->surf->bvh;

		//execute cooperative processing of samplinggridpoints
		#pragma omp for schedule(dynamic,1) nowait
		for( auto it = samplingpoints.begin(); it != samplingpoints.end(); ++it ) {

			AABB e_aabb( 	trpl(it->x-R, it->y-R, it->z-R),
							trpl(it->x+R, it->y+R, it->z+R) );

			vector<unsigned int> triangle_cand = mylinkedbvh->query( e_aabb ); //MK::assure that this is read only ##MK future improve for better locality multiple Rtrees

			if ( triangle_cand.size() == 0 ) {

				//no candidates so far enough from tip surface to compute unbiased ROI
				//perform analysis per single point
				//define a measurement grid
				p3dm1 probesite = p3dm1(it->x, it->y, it->z, UNKNOWNTYPE);

				//check if point is sufficiently deep enough in tip volume use owner->surf triangle tree for this
				//test if bounding box about spherical inspection volume cuts any triangle or not
				//MK::strictly speaking this is conversatively strong, all triangles might protrude into aabb3d but not cut the sphere
				//however, setting yourself stricter aims is not a bad thing

				//find candidate atoms in spherical environment
				//each member in the team reads details of specific region and initializes threadlocal comparator variables
				//find in which thread region the point is and candidate points
				vector<p3dm1> neighbors3dm1;
				//scan region about this point to find neighboring candidates using KDtree O(lgN)
				for( int thr = MASTER; thr < nt; thr++ ) {
					threadmemory* currentregiondata = owner->sp->db.at(thr);
					vector<p3dm1> const & theseions = currentregiondata->ionpp3_kdtree;
					kd_tree* curr_kauri = currentregiondata->threadtree;
					curr_kauri->range_rball_noclear_nosort_p3d( probesite, theseions, RSQR, neighbors3dm1 );
				}

				if ( neighbors3dm1.size() > 0 ) { //if any, go to next point, if some make angular space sweeping + histogramming + FFT

					thrindexer->res.push_back( vics_materialpoint_result( p3d(it->x, it->y, it->z) ) );

					//precompute difference vectors in preparation for d = dot(n,(p-po)) were po is here and p the elements in neighbors3dm1
					vector<d3d> distvec;
					for(auto jt = neighbors3dm1.begin(); jt != neighbors3dm1.end(); ++jt) { //O(n)
						distvec.push_back( d3d(jt->x - it->x, jt->y - it->y, jt->z - it->z) );
					}
					/*//##MK::bSIMD
					vector<float> distvec;
					distvec.reserve( 3*neighbors3dm1.size() ); //three coordinate values per point
					for(auto jt = neighbors3dm1.begin(); jt != neighbors3dm1.end(); ++jt) {
						distvec.push_back(static_cast<float>(jt->x - it->x));
						distvec.push_back(static_cast<float>(jt->y - it->y));
						distvec.push_back(static_cast<float>(jt->z - it->z));
					}*/

					//perform analysis sample elevation/azimuth space for one point in reconstruction space
					thrindexer->vicsmethod_core_incrementalfft2( *it, distvec ); //N_elev*N_azi*O(NlgN) time complexity
					//thrindexer->vicsmethod_core_batchfft_simd( *it, distvec );
				}
				//do not report for voids in dataset

			} //next samplingpoint
			else {
				/*#pragma omp critical
				{
					cout << "Samplingpoint " << it->x << ";" << it->y << ";" << it->z << " is too close to dataset boundary for being unbiased" << "\n";
				}*/
			}

		} //next sampling point
		#pragma omp barrier

		double mytoc = omp_get_wtime();
		#pragma omp critical
		{
			cout << "Thread " << mt << " finished atom probe crystallography took " << (mytoc-mytic) << " seconds" << endl;
		}
	} //end of parallel region


	double toc = MPI_Wtime();
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ExtractCrystallographyAraullo", APT_UTL, APT_IS_PAR, mm, tic, toc);
	cout << "Conducted multithreaded AtomProbeCrystallography analysis in " << (toc-tic) << " seconds" << endl;
}


void aptcrystHdl::report_crystallography_results()
{
	double tic = MPI_Wtime();

	string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".APTCrystalloPointWiseResults.csv";
	ofstream tlog;
	tlog.open( fn.c_str() );

//	tlog << "ThreadID;X;Y;Z;Elevation;Azimuth;FFTAllocYes;FFTPlacementYes;FFTInitYes;FFTCompYes;FFTFreeYes;Bin1;SignalStrength1;Bin2;SignalStrength2;Bin3;SignalStrength3\n";
//	tlog << ";nm;nm;nm;rad;rad;bool;bool;bool;bool;bool;1;##;1;##;1;##\n";
//	tlog << "ThreadID;X;Y;Z;Elevation;Azimuth;FFTAllocYes;FFTPlacementYes;FFTInitYes;FFTCompYes;FFTFreeYes;Bin1;SignalStrength1;Bin2;SignalStrength2;Bin3;SignalStrength3\n";

	tlog << "ThreadID;X;Y;Z;Elevation;Azimuth;Bin1;SignalStrength1;Bin2;SignalStrength2;Bin3;SignalStrength3\n";
	tlog << "1;nm;nm;nm;rad;rad;1;a.u.;1;a.u.;1;a.u.\n";
	tlog << "ThreadID;X;Y;Z;Elevation;Azimuth;Bin1;SignalStrength1;Bin2;SignalStrength2;Bin3;SignalStrength3\n";

	tlog << setprecision(8);
	for ( size_t thr = 0; thr < workers.size(); ++thr ) {
		for ( auto it = workers.at(thr)->res.begin(); it != workers.at(thr)->res.end(); ++it ) {
			for ( auto jt = it->ElevAzimTblSummary.begin(); jt != it->ElevAzimTblSummary.end(); ++jt ) {
				tlog << thr << ";" << it->MatPointPos.x << ";" << it->MatPointPos.y << ";" << it->MatPointPos.z;
				tlog << ";" << jt->elevation << ";" << jt->azimuth;
				//tlog << ";" << jt->fft_allc_success << ";" << jt->fft_plac_success;
				//tlog << ";" << jt->fft_init_success << ";" << jt->fft_comp_success << ";" << jt->fft_free_success;
				tlog << ";" << jt->max1.first << ";" << jt->max1.second;
				tlog << ";" << jt->max2.first << ";" << jt->max2.second;
				tlog << ";" << jt->max3.first << ";" << jt->max3.second << "\n";
			}
			//free resources
			//it->ElevAzimTblSummary = vector<vicsresult>();
		}
	}

	tlog.flush();
	tlog.close();

	fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".APTCrystalloGridSummary.csv";
	ofstream slog;
	slog.open( fn.c_str() );

	slog << "ThreadID;X;Y;Z;SuccessfulFFTs;FailedFFTs;NumberOfIons\n";
	slog << ";nm;nm;nm;1;1;1\n";
	slog << "ThreadID;X;Y;Z;SuccessfulFFTs;FailedFFTs;NumberOfIons\n";

	for ( size_t thr = 0; thr < workers.size(); ++thr ) {
		for ( auto it = workers.at(thr)->res.begin(); it != workers.at(thr)->res.end(); ++it ) {
			slog << thr << ";" << it->MatPointPos.x << ";" << it->MatPointPos.y << ";" << it->MatPointPos.z;
			slog << ";" << it->FFTSummary.nFFTsSuccess << ";" << it->FFTSummary.nFFTsFailed;
			slog << ";" << it->FFTSummary.nIons << "\n";
		}
	}

	slog.flush();
	slog.close();

	double toc = MPI_Wtime();
	//memory is still high because all intermediate crystallography results still there
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ReportCrystallographyAraullo", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	cout << "Wrote AtomProbeCrystallography results file in " << (toc-tic) << " seconds" << endl;
}


void aptcrystHdl::report_crystallography_results2()
{
	double tic = MPI_Wtime();

	int status = 0;
	//report global cloud of elevation/azimuth pairs PARAPROBE_CRYSTALLO_ELEVAZIMMETA
	vector<float> wfbuf;
	wfbuf.resize(2*elevazimpairs.size()); //total number of elements 1E5-1E6 a 2*4B so <=8MB
	size_t i = 0;
	for( i = 0; i < elevazimpairs.size();  ) {
		wfbuf[2*i+0] = elevazimpairs[i].x;
		wfbuf[2*i+1] = elevazimpairs[i].y; //xy columns each pair a row
		i++;
	}
	h5iometa ifo = h5iometa( PARAPROBE_CRYSTALLO_ELEVAZIMMETA, elevazimpairs.size(), 2 );
	status = owner->owner->crysth5Hdl.create_contiguous_matrix_f32le( ifo );
	h5offsets offs = h5offsets( 0, elevazimpairs.size(), 0, 2, elevazimpairs.size(), 2);
	status = owner->owner->crysth5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, wfbuf );

	//report sampling material point cloud PARAPROBE_CRYSTALLO_MATPOINTMETA
	size_t NumberOfPointsWithResults = 0;
	for ( size_t thr = 0; thr < workers.size(); ++thr ) {
		NumberOfPointsWithResults += workers.at(thr)->res.size();
	}
	wfbuf.resize(3*NumberOfPointsWithResults); //total number of elements 5E6 a 3*4B <= 60MB
	i = 0;
	for( size_t thr = 0; thr < workers.size(); ++thr ) {
		for( auto it = workers.at(thr)->res.begin(); it != workers.at(thr)->res.end(); ++it ) {
			wfbuf[i+0] = it->MatPointPos.x;
			wfbuf[i+1] = it->MatPointPos.y;
			wfbuf[i+2] = it->MatPointPos.z;
			i = i + 3;
		}
	}
	ifo = h5iometa( PARAPROBE_CRYSTALLO_MATPOINTXYZ, NumberOfPointsWithResults, 3 );
	status = owner->owner->crysth5Hdl.create_contiguous_matrix_f32le( ifo );
	offs = h5offsets( 0, NumberOfPointsWithResults, 0, 3, NumberOfPointsWithResults, 3 );
	status = owner->owner->crysth5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, wfbuf );
	wfbuf.clear();

	vector<unsigned int> wuibuf;
	wuibuf.resize(4*NumberOfPointsWithResults); //total number of elements 5E6 a 4*4B so <= 80MB
	i = 0;
	for( size_t thr = 0; thr < workers.size(); ++thr ) {
		for( auto it = workers.at(thr)->res.begin(); it != workers.at(thr)->res.end(); ++it ) {
			wuibuf[i+0] = it->FFTSummary.threadID;
			wuibuf[i+1] = it->FFTSummary.nFFTsFailed;
			wuibuf[i+2] = it->FFTSummary.nFFTsSuccess;
			wuibuf[i+3] = it->FFTSummary.nIons;
			i = i + 4;
		}
	}
	ifo = h5iometa( PARAPROBE_CRYSTALLO_MATPOINTMETA, NumberOfPointsWithResults, 4 );
	status = owner->owner->crysth5Hdl.create_contiguous_matrix_u32le( ifo );
	offs = h5offsets( 0, NumberOfPointsWithResults, 0, 4, NumberOfPointsWithResults, 4 );
	status = owner->owner->crysth5Hdl.write_contiguous_matrix_u32le_hyperslab( ifo, offs, wuibuf );
	wuibuf.clear();

	//report individual material point results
	if ( Settings::IOCrystallography == true ) {
		size_t pid = 0;
		for( size_t thr = 0; thr < workers.size(); ++thr ) {
			for( auto it = workers.at(thr)->res.begin(); it != workers.at(thr)->res.end(); ++it ) {
				//fill buffer with results for a single material point dump to file next one
				i = 0;
				wfbuf.resize(6*it->ElevAzimTblSummary.size());
				for ( auto jt = it->ElevAzimTblSummary.begin(); jt != it->ElevAzimTblSummary.end(); ++jt ) {
					wfbuf[i+0] = jt->max1.first;
					wfbuf[i+1] = jt->max1.second;
					wfbuf[i+2] = jt->max2.first;
					wfbuf[i+3] = jt->max2.second;
					wfbuf[i+4] = jt->max3.first;
					wfbuf[i+5] = jt->max3.second;
					i = i + 6;
				}
				string fwslash = "/";
				string dsfn = PARAPROBE_CRYSTALLO_THREESTRONGEST + fwslash + to_string(pid);
				ifo = h5iometa( dsfn, it->ElevAzimTblSummary.size(), 6 );
				status = owner->owner->crysth5Hdl.create_contiguous_matrix_f32le( ifo );
				offs = h5offsets( 0, it->ElevAzimTblSummary.size(), 0, 6, it->ElevAzimTblSummary.size(), 6 );
				status = owner->owner->crysth5Hdl.write_contiguous_matrix_f32le_hyperslab( ifo, offs, wfbuf );
				//bookkeeping for next dataset
				pid++;
			}
		}
		wfbuf.clear();
	}

	double toc = MPI_Wtime();
	//memory is still high because all intermediate crystallography results still there
	memsnapshot mm = owner->owner->tictoc.get_memoryconsumption();
	owner->owner->tictoc.prof_elpsdtime_and_mem( "ReportHDF5CrystallographyAraullo", APT_IO, APT_IS_SEQ, mm, tic, toc);
	cout << "Wrote AtomProbeCrystallography results to H5 file in " << (toc-tic) << " seconds" << endl;
}


void aptcrystHdl::delete_crystallography_results()
{
	for(size_t i = 0; i < workers.size(); i++) {
		delete workers.at(i);
	}
}


solverHdl::solverHdl()
{
	nevt = 0;

	myRank = MASTER;
	nRanks = SINGLEPROCESS;
}


solverHdl::~solverHdl()
{
	delete_rawdata();
}


bool solverHdl::generate_synthetic_tip()
{
	//##MK::so far single crystalline tip only
	//##MK::implement range check and consistence of planned tip geometry
	double tic, toc;

	tic = MPI_Wtime();

	//##MK::so far generate axis-aligned tips with tip axis z parallel to c axis of fcc lattice and z coordinate

	//get planned tip geometry
	geomodel tipgeometry = geomodel(
			Settings::RelBottomRadius, Settings::RelTopRadius,
			Settings::RelBottomCapHeight, Settings::RelTopCapHeight,
			Settings::LatticeConstant, Settings::NumberOfAtoms);

	//knowing the tip geometry we can now allocate a O(n) time access structure
	rawdata_alf.initcontainer( tipgeometry.mybox );

	//cout << "rawdata_alf/imi/imx/xsz/nx\t\t" << rawdata_alf.mdat.box.xmi << "\t\t" << rawdata_alf.mdat.box.xmx << "\t\t" << rawdata_alf.mdat.box.xsz << "\t\t" << rawdata_alf.mdat.nx << endl;
	//cout << "rawdata_alf/imi/imx/ysz/ny\t\t" << rawdata_alf.mdat.box.ymi << "\t\t" << rawdata_alf.mdat.box.ymx << "\t\t" << rawdata_alf.mdat.box.ysz << "\t\t" << rawdata_alf.mdat.ny << endl;
	//cout << "rawdata_alf/imi/imx/zsz/nz\t\t" << rawdata_alf.mdat.box.zmi << "\t\t" << rawdata_alf.mdat.box.zmx << "\t\t" << rawdata_alf.mdat.box.zsz << "\t\t" << rawdata_alf.mdat.nz << endl;

	//get a solid solution model
	solutemodel dilute = solutemodel();

	//get uvw lattice plane limits to fully enclose tip AABB, ##MK::exemplarily fcc crystal lattice only here
	unitcellaggr rve = unitcellaggr( Settings::LatticeConstant, tipgeometry.mybox, FCC );

	//account for limited detection efficiency and finite spatial resolution
	mt19937 dice;
	dice.seed( Settings::RndDescrStatsPRNGSeed );
	uniform_real_distribution<apt_real> unifrnd(0.f,1.f);
	for (unsigned int i = 0; i < Settings::RndDescrStatsPRNGDiscard; ++i) //warm up MersenneTwister
		apt_real discard = unifrnd(dice);

	//get normally distributed value
	normal_distribution<apt_real> gaussian_x(0.f, Settings::SpatResolutionSigmaX ); //mean is zero
	normal_distribution<apt_real> gaussian_y(0.f, Settings::SpatResolutionSigmaY );
	normal_distribution<apt_real> gaussian_z(0.f, Settings::SpatResolutionSigmaZ );


	//MK::on the fly addition of ions from crystallographic motif to prevent building unnecessary copies
	size_t placed = 0;
	//apt_real placed = 0.f;
	//apt_real progress = 1.f; //##MK::1% increment
	//apt_real total = 1.f / static_cast<apt_real>(Settings::NumberOfAtoms) * 100.f;
	for(size_t b = 0; b < rve.base.size(); ++b) {
		for(int w = rve.wmin; w <= rve.wmax; ++w) {
			for(int v = rve.vmin; v <= rve.vmax; ++v) {
				for(int u = rve.umin; u <= rve.umax; ++u) {
					p3d ap = rve.get_atom(b, u, v, w);
					if ( tipgeometry.is_inside( ap ) == true ) {
						//account for limited detection efficiency
						if ( unifrnd(dice) > Settings::SimDetEfficiency ) //every eta-th gets through so every 1-eta not, if eta = 0.0, almost every dice result > 0.0 discard, if eta --> 1 almost no dice result is > 1-small number
							continue;

//cout << ap << endl;
						//account for finite and anisotropic spatial resolution
						ap.x += gaussian_x(dice); //############MK::limit position...
						ap.y += gaussian_y(dice);
						ap.z += gaussian_z(dice);

//cout << ap << endl;

						//mimick a dilute solid solution alloy without spatial correlations or clustering
						//unsigned int typid = dilute.get_random_speci_typid(); //evaluate a PRNG with the solutemodel to get a random lattice ion
						apt_real mass2charge = dilute.get_random_speci_mq();

						//change to adding into a constant time query container
						rawdata_alf.add_atom( pos(ap.x, ap.y, ap.z, mass2charge) );

						//placed += 1.f;
						placed++;
						if ( placed % 1000000 != 0 )
							continue;
						else
							cout << placed << endl;
					}
				}
			}
		}
	}

/*
	//##MK::naive random point process carving
	size_t placed = 0;
	for( size_t i = 0; i < Settings::NumberOfAtoms; ++i ) {
		p3d ap = p3d(	unifrnd(dice) * tipgeometry.mybox.xsz + tipgeometry.mybox.xmi,
						unifrnd(dice) * tipgeometry.mybox.ysz + tipgeometry.mybox.ymi,
						unifrnd(dice) * tipgeometry.mybox.zsz + tipgeometry.mybox.zmi );

		if ( tipgeometry.is_inside( ap ) == true ) {
			//no-accounting for limited detection efficiency
			//no-accounting for finite and anisotropic spatial resolution
			apt_real mass2charge = dilute.get_random_speci_mq();
			rawdata_alf.add_atom( pos(ap.x, ap.y, ap.z, mass2charge) );

			placed++;
			if ( placed % 1000000 != 0 )
				continue;
			else
				cout << placed << endl;
		}
	}
*/

	reporting( "All tip atoms were placed");

	//##MK  ##potentially deleting the data? rawdata_alf.write_occupancy_raw();

	//add spherical precipitates randomly placed, MK::MIND THAT DIMENSIONS IN NANOMETER
	secondphasemodel ballpark = secondphasemodel(
			tipgeometry,
			Settings::NumberOfCluster,
			Settings::ClusterRadiusMean,
			Settings::ClusterRadiusSigmaSqr );

	if ( ballpark.particles.size() > 0 ) {
		ballpark.reportParticles( Settings::SimID, get_rank() );
		ballpark.reportParticlesVTK( Settings::SimID, get_rank() );
	}

	reporting( "Geometry of cluster/precipitates defined");

	vector<pos> newatoms;
	for(auto it = ballpark.particles.begin(); it != ballpark.particles.end(); ++it ) {
		//erode all host atoms falling inside the spherical ball

		//##MK::we do not relax the structure via MD here!
		p3d ballcenter = it->center;
		apt_real radius = it->radius;

//cout << "Ballcenter/Radius\t\t" << ballcenter.x << ";" << ballcenter.y << ";" << ballcenter.z << "\t\t" << radius << endl;

		erase_log status = rawdata_alf.erase_rball( ballcenter, radius );

		//place instead the atoms within the unit cell aggregate of the spherical cluster
		newatoms.clear();
		it->get_atoms( newatoms ); //the precipitate knows its size

		size_t added = 0;
		for( size_t j = 0; j < newatoms.size(); ++j) {
			p3d cand = p3d( newatoms.at(j).x, newatoms.at(j).y, newatoms.at(j).z );
			if ( tipgeometry.is_inside( cand ) == true ) { //also in reality cluster/precipitates are usually cut by tip boundary
				added++;
				if ( unifrnd(dice) > Settings::SimDetEfficiency ) //every eta-th gets through so every 1-eta not, if eta = 0.0, almost every dice result > 0.0 discard, if eta --> 1 almost no dice result is > 1-small number
					continue;

				//spatial detection inaccuracy
				pos ap = newatoms.at(j);

				ap.x += gaussian_x(dice);  //#####check strictly speaking as there is no relaxatrion and gaussian_i free atoms can overlap...
				ap.y += gaussian_y(dice);
				ap.z += gaussian_z(dice);

				rawdata_alf.add_atom( ap );

			} //placement attempt of included atom in current cluster
		} //next atom of current cluster

cout << "Cluster placed cleared/kept/newatoms planned/new atoms placed\t\t" << status.ncleared << ";" << status.nkept << ";" << newatoms.size() << ";" << added << "\n";
	} //next cluster

	//##MK::all more advanced tip synthesis stuff comes as post-processing
	//##MK::i.e. carve out from the matrix spherical volume and place spherical or ellipsoidal precipitates

	//how many ions in total?
	nevt = 0;
	for (size_t b = 0; b < rawdata_alf.buckets.size(); ++b) {
		if ( rawdata_alf.buckets.at(b) != NULL ) {
			nevt = nevt + static_cast<unsigned int>(rawdata_alf.buckets.at(b)->size()); //##MK::increase size
		}
	}

	toc = MPI_Wtime();
	memsnapshot mm = tictoc.get_memoryconsumption();
	tictoc.prof_elpsdtime_and_mem( "GenerateSyntheticTip", APT_UTL, APT_IS_SEQ, mm, tic, toc);
	string mess = "Sequential tip synthesis was successful and took " + to_string(toc-tic) + " seconds!";
	reporting( get_rank(), mess );


#ifdef UTILIZE_HDF5
	if ( Settings::IOReconstruction == true ) { //store generated tip into an HDF5 file
		tic = MPI_Wtime();

		//##MK::currently switched off
		//write_pos_hdf5( rawdata_alf.buckets, Settings::RAWFilenameIn );

		toc = MPI_Wtime();
		memsnapshot mm = memsnapshot();
		tictoc.prof_elpsdtime_and_mem( "WriteSyntheticTipHDF", APT_IO, APT_IS_SEQ, mm, tic, toc);
		string mess = "Writing synthetic tip to HDF5 was successful and took " + to_string(toc-tic) + " seconds!";
		reporting( get_rank(), mess );
	}
#endif

	return true;
}


bool solverHdl::load_pos_sequentially( const string posfn )
{
	//loads pos file
	//B. Gault et al., Atom Probe Tomography, Springer, 2012, Appendix ##MK
	//*.pos is headerless container of binary big endian float structs (IEEE754?) four to each ion
	//x,y,z, mq interleaved format

	double tic = MPI_Wtime();

//##MK::templatize this function in the future
//##MK::put this section into individual function it is a duplicate from the one in load_epos...
	//get total number of ions by reading file size via OS intrinsics
	unsigned long filesize[1] = {0};
	unsigned long nobjects[1] = {0};
	struct stat oscall;
	if ( stat( posfn.c_str() , &oscall ) != -1 ) { //MK::take care of case when filesize exceeds range of unsigned int...
		filesize[0] = static_cast<unsigned long>(oscall.st_size);
		if ( filesize[0] % static_cast<unsigned long>(sizeof(struct pos_be_2_stl_le)) != 0 ) {
			stopping( "Size of file " + posfn + " is not an integer multiple of expected datatype");
			return false;
		}
		else {
			nobjects[0] = filesize[0] / static_cast<unsigned long>(sizeof(struct pos_be_2_stl_le));
			string mess = "Dataset " + posfn + " contains " + to_string(nobjects[0]) + " ions";
			reporting( myRank, mess);
		}

		if ( nobjects[0] > static_cast<unsigned long>(std::numeric_limits<unsigned int>::max()-2) ) {
			string mess = "ERROR::Dataset is larger than what the code currently can process (4.2 billion ions, (not 4.2 billion byte...). Change data type to unsigned long!";
			stopping( get_rank(), mess );
			//##MK::can be changed by going from unsigned int to unsigned long, but then also the HoshenKopelman needs upgrade, new rawdata structure of lists of buckets of epos objects supports larger than 4.2 billion ions...
			return false;
		}
	}
	else {
		string mess = "File " + posfn + " cannot be read!";
		stopping( get_rank(), mess );
		return false;
	}
//##MK::put this section into individual function it is a duplicate from the one in load_epos...

	//assume length of all containers the same, total number of ions
	nevt = nobjects[0]; //MK::implicit contraction from unsigned long and unsigned int is now okay

	//MK::read successively blocks of C structs and store in blocks the rawdata
	//benefit is that such later rawdata can be processes blockwise and deleted successively reducing the total memory footprint of the application for reconstruction by not more than though factor 2
	size_t efirst = 0;
	size_t elast = static_cast<size_t>(nevt);
	size_t etarget = MPIIO_READ_CACHE / sizeof(struct pos_be_2_stl_le);
	size_t ecurrent = 0;

	struct pos_be_2_stl_le* rbuf = NULL;

	//file exists because size query did not return if reaching this point
	unsigned long probe = 0;
	FILE *file = fopen( posfn.c_str(), "rb" );

	//read block by block
	for ( size_t e = efirst; e < elast; 	) { //number of elements, not total bytes
		//update actual read cache length
		ecurrent = ((e + etarget) < elast) ? etarget : elast - e;

		//allocate buffer if none existent, reallocate if of different size
		if ( rbuf != NULL ) {
			if ( ecurrent != etarget ) {
				//reallocate buffer of size ecurrent
				delete [] rbuf; rbuf = NULL;
				try { rbuf = new struct pos_be_2_stl_le[ecurrent]; }
				catch (bad_alloc &ioexc) {
					stopping( get_rank(), "Unable to reallocate rbuf for pos_be_2_stl_le structs in epos_sequential read" );
					fclose(file);
					return false;
				}
			}
			//else nothing to do reutilize rbuf
		}
		else {
			//least likely case first time allocation
			try { rbuf = new struct pos_be_2_stl_le[ecurrent]; }
			catch (bad_alloc &ioexc) {
				stopping( get_rank(), "Unable to allocate rbuf for pos_be_2_stl_le structs in epos_sequential read" );
				fclose(file);
				return false;
			}
		}
		//by now buffer allocated or already out of subroutine

		//read from file
		probe = fread( rbuf, sizeof(struct pos_be_2_stl_le), ecurrent, file ); //file pointer is advanced implicitly/automatically

		if ( probe == ecurrent ) { //read success
			//allocate new buffer to store results
			vector<pos>* wbuf = NULL;
			try { wbuf = new vector<pos>; }
			catch (bad_alloc &ioexc) {
				stopping( get_rank(), "Unable to allocate wbuf for epos vector in pos_sequential read" );
				delete [] rbuf;
				fclose(file);
				return false;
			}
			wbuf->reserve( ecurrent );

			for ( size_t i = 0; i < ecurrent; ++i ) { //swop endianness implicitly within epos object
				wbuf->push_back( pos(rbuf[i]) );
//cout << wbuf->back() << endl;
			}

			//handshake registration data block wbuf in solverHdl
			rawdata_pos.push_back( NULL );
			rawdata_pos.back() = wbuf;

			e = e + ecurrent;
//cout << e << endl;
		} //block read of ecurrent successful
		else {
			string mess = "I/O failed on " + posfn + " at position " + to_string(e);
			stopping( get_rank(), mess );
			delete [] rbuf;
			fclose(file);
			return false;
		}
	} //next block
	fclose(file);

	delete [] rbuf; rbuf = NULL;

	double toc = MPI_Wtime();
	memsnapshot mm = tictoc.get_memoryconsumption();
	tictoc.prof_elpsdtime_and_mem( "POSFileInputWithEndiannessSwop", APT_IO, APT_IS_SEQ, mm, tic, toc);
	string mess = "Sequential binary-based data loading with BE to LE endianness swopping was successful and took " + to_string(toc-tic) + " seconds!";
	reporting( get_rank(), mess );

	return true;
}


bool solverHdl::load_epos_sequentially( const string eposfn )
{
	//loads native epos file
	//D.J. Larson et al., Local Electrode Atom Probe Tomography: A User\'s Guide
	//DOI 10.1007/978-1-4614-8721-0, Springer Science+Business Media New York 2013 Appendix A
	//*.epos floats are big endian, int as well, file has no header it is just a naive block of the same structured data
	//therefore we can read the file via C char structs and swop bytes for int
	//for floats this is usually dangerous, though, but fact that IVAS writes IEEE754 comes to help
	//see details in functionality of epos struct

	double tic = MPI_Wtime();

	//get total number of ions by reading file size via OS intrinsics
	unsigned long filesize[1] = {0};
	unsigned long nobjects[1] = {0};
	struct stat oscall;
	if ( stat( eposfn.c_str() , &oscall ) != -1 ) { //MK::take care of case when filesize exceeds range of unsigned int...
		filesize[0] = static_cast<unsigned long>(oscall.st_size);
		if ( filesize[0] % static_cast<unsigned long>(sizeof(struct epos_be_2_stl_le)) != 0 ) {
			stopping( "Size of file " + eposfn + " is not an integer multiple of expected datatype");
			return false;
		}
		else {
			nobjects[0] = filesize[0] / static_cast<unsigned long>(sizeof(struct epos_be_2_stl_le));
			string mess = "Dataset " + eposfn + " contains " + to_string(nobjects[0]) + " ions";
			reporting( myRank, mess);
		}

		if ( nobjects[0] > static_cast<unsigned long>(std::numeric_limits<unsigned int>::max()-2) ) {
			string mess = "ERROR::Dataset is large than what the code currently can process (4.2 billion ions, (not 4.2 billion byte...). Change data type to unsigned long!";
			stopping( get_rank(), mess );
			//##MK::can be changed by going from unsigned int to unsigned long, but then also the HoshenKopelman potentially needs upgrade, new rawdata structure of lists of buckets of epos objects supports larger than 4.2 billion ions...
			return false;
		}
	}
	else {
		string mess = "File " + eposfn + " cannot be read!";
		stopping( get_rank(), mess );
		return false;
	}

	//assume length of all containers the same, total number of ions
	nevt = nobjects[0]; //MK::implicit contraction from unsigned long and unsigned int is now okay

	//MK::read successively blocks of C structs and store in blocks the rawdata
	//benefit is that such later rawdata can be processes blockwise and deleted successively reducing the total memory footprint of the application for reconstruction by not more than though factor 2
	size_t efirst = 0;
	size_t elast = static_cast<size_t>(nevt);
	size_t etarget = MPIIO_READ_CACHE / sizeof(struct epos_be_2_stl_le);
	size_t ecurrent = 0;

	struct epos_be_2_stl_le* rbuf = NULL;

	//file exists because size query did not return if reaching this point
	unsigned long probe = 0;
	FILE *file = fopen( eposfn.c_str(), "rb" );

	//read block by block
	for ( size_t e = efirst; e < elast; 	) { //number of elements, not total bytes
		//update actual read cache length
		ecurrent = ((e + etarget) < elast) ? etarget : elast - e;

		//allocate buffer if none existent, reallocate if of different size
		if ( rbuf != NULL ) {
			if ( ecurrent != etarget ) {
				//reallocate buffer of size ecurrent
				delete [] rbuf; rbuf = NULL;
				try { rbuf = new struct epos_be_2_stl_le[ecurrent]; }
				catch (bad_alloc &ioexc) {
					stopping( get_rank(), "Unable to reallocate rbuf for epos_be_2_stl_le structs in epos_sequential read" );
					fclose(file);
					return false;
				}
			}
			//else nothing to do reutilize rbuf
		}
		else {
			//least likely case first time allocation
			try { rbuf = new struct epos_be_2_stl_le[ecurrent]; }
			catch (bad_alloc &ioexc) {
				stopping( get_rank(), "Unable to allocate rbuf for epos_be_2_stl_le structs in epos_sequential read" );
				fclose(file);
				return false;
			}
		}
		//by now buffer allocated or already out of subroutine

		//read from file
		probe = fread( rbuf, sizeof(struct epos_be_2_stl_le), ecurrent, file ); //file pointer is advanced implicitly/automatically

		if ( probe == ecurrent ) { //read success
			//allocate new buffer to store results
			vector<epos>* wbuf = NULL;
			try { wbuf = new vector<epos>; }
			catch (bad_alloc &ioexc) {
				stopping( get_rank(), "Unable to allocate wbuf for epos vector in epos_sequential read" );
				delete [] rbuf;
				fclose(file);
				return false;
			}
			wbuf->reserve( ecurrent );

			for ( size_t i = 0; i < ecurrent; ++i ) { //swop endianness implicitly within epos object
				wbuf->push_back( epos(rbuf[i]) );
//cout << wbuf->back() << endl;
			}

			//handshake registration data block wbuf in solverHdl
			rawdata_epos.push_back( NULL );
			rawdata_epos.back() = wbuf;

			e = e + ecurrent;
//cout << e << endl;
		} //block read of ecurrent successful
		else {
			string mess = "I/O failed on " + eposfn + " at position " + to_string(e);
			stopping( get_rank(), mess );
			delete [] rbuf;
			fclose(file);
			return false;
		}
	} //next block
	fclose(file);

	delete [] rbuf; rbuf = NULL;

	double toc = MPI_Wtime();
	memsnapshot mm = tictoc.get_memoryconsumption();
	tictoc.prof_elpsdtime_and_mem( "EPOSFileInputWithEndiannessSwop", APT_IO, APT_IS_SEQ, mm, tic, toc);
	string mess = "Sequential binary-based data loading with BE to LE endianness swopping was successful and took " + to_string(toc-tic) + " seconds!";
	reporting( get_rank(), mess );

	return true;
}


#ifdef UTILIZE_HDF5

bool solverHdl::load_hdf5_sequentially( const string h5fn )
{
	//##MK::loads H5 file, no necessity for byte swopping as type information is included in the H5 file
	//also no necessity for reading other meta because H5 file is self-contained
	//reading the total number of atoms is required because chunking will generate blocks not necessarily
	//an integer divisor of the number of atoms. Instead, last/trailing chunk has some dummy data we filter off
	//Format is interleaved [NChunks*AtomsPerChunk,4] type H5T_IEEE_F32LE floats, x, y, z, mass2charge
	double tic = MPI_Wtime();

	//open H5 file and read number of atoms
	herr_t status;
	unsigned int rdata[1] = { 0 };

	hid_t fileid = H5Fopen( h5fn.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	if ( fileid < 0 ) {
		stopping( get_rank(), "Loading the HDF5 file was not successful"); return false;
	}
	hid_t dsetid = H5Dopen2( fileid, "/NumberOfIons", H5P_DEFAULT);
	status = H5Dread( dsetid, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
	if ( status < 0 ) {
		status = H5Dclose(dsetid);
		status = H5Fclose(fileid);
		stopping( get_rank(), "Loading the /NumberOfIons data entry was not successful"); return false;
	}
	//reading number of atoms success, release resources and remember number
	status = H5Dclose(dsetid);
	nevt = rdata[0];

	//next get data layout of XYZMQ within H5 file
	dsetid = H5Dopen2( fileid, "/XYZMQ", H5P_DEFAULT );
	hid_t fspcid = H5Dget_space(dsetid);
	//int rank = H5Sget_simple_extent_ndims(fspcid);
	hsize_t rdims[2];
	status = H5Sget_simple_extent_dims(fspcid, rdims, NULL);

	cout << "HDF5 XYZMQ loading NumberOfAtoms/rdims[0]/rdims[1]\t\t" << rdata[0] << "\t\t" << rdims[0] << "\t\t" << rdims[1] << endl;

	//read chunking information
	hid_t cparms = H5Dget_create_plist(dsetid);
	hid_t mspcid;
	float* rbuf = NULL;

	if ( H5D_CHUNKED == H5Pget_layout(cparms) ) {

		hsize_t chunkdims[2];
		int rankchunk = H5Pget_chunk(cparms, 2, chunkdims);

		cout << "HDF5 XYZMQ loading chunkdim[0]/chunkdim[1]\t\t" << chunkdims[0] << "\t\t" << chunkdims[1] << endl;

		mspcid = H5Screate_simple(rankchunk, chunkdims, NULL);

		//get rbuf to store all chunks read into rawdata_pos DO USE the pos(x,y,z,mq) constructor NOT the one with byte swopping !
		size_t NumberOfChunks = static_cast<size_t>(rdims[0]) / static_cast<size_t>(chunkdims[0]);
		size_t elementsPerChunk = chunkdims[0] * chunkdims[1];

		cout << "HDF5 XYZMQ loading NumberOfChunks/elementsPerChunk\t\t" << NumberOfChunks << "\t\t" << elementsPerChunk << endl;

		bool allocation_status = true;
		try {
			//allocate buffer to read a single chunk from the H5 file and interpret
			rbuf = new float[elementsPerChunk];
			for( size_t i = 0; i < elementsPerChunk; ++i )
				rbuf[i] = 0.f;
		}
		catch (bad_alloc &h5croak) {
			stopping( get_rank(), "Allocation error while preparing chunk cache read buffer in load_hdf5_sequentially" );
			allocation_status = false;
		}

		//given that we know that the chunk dimensions and we know the total number of chunks we can allocate, rawdata_pos container also
		size_t chk = 0;
		while ( chk < NumberOfChunks && allocation_status == true ) {
			vector<pos>* wbuf = NULL;
			try {
				wbuf = new vector<pos>;
				wbuf->reserve( chunkdims[0] );
			}
			catch (bad_alloc &h5croak) {
				stopping( get_rank(), "Allocation error while preparing rawdata_pos read buffer in load_hdf5_sequentially" );
				allocation_status = false;
				break;
			}
			rawdata_pos.push_back(NULL); //allocation was successful
			rawdata_pos.back() = wbuf;
			chk++;
		}

		if ( allocation_status == false ) {
			stopping( get_rank(), "Stopping and closing HDF5 file after experiencing errors");
			status = H5Sclose(mspcid);
			status = H5Dclose(dsetid);
			status = H5Sclose(fspcid);
			status = H5Fclose(fileid);
			return false;
		}
		else {
			reporting( get_rank(), "rawdata_pos read buffer prepared");
		}

		//read successively the chunks from the HDF5 file and interpret, take care to avoid reading in dummy data from last chunk
		hsize_t offset[2] = {0, 0};
		hsize_t count[2] = {0, 0};
		bool read_status = true;
		size_t AtomsToReadTotal = nevt;
		size_t AtomsToReadCurr = 0;

		for ( chk = 0; chk < NumberOfChunks; ++chk ) {

			offset[0] = static_cast<hsize_t>(chk) * chunkdims[0];
			offset[1] = 0;
			count[0] = chunkdims[0];
			count[1] = chunkdims[1];
			status = H5Sselect_hyperslab( fspcid, H5S_SELECT_SET, offset, NULL, count, NULL);
			if ( status < 0 ) { read_status = false; break; }
			status = H5Dread( dsetid, H5T_IEEE_F32LE, mspcid, fspcid, H5P_DEFAULT, rbuf );
cout << "-->Reading chunk chk/offset[0]/[1]/count[0]/[1]/status\t\t" << chk << "\t\t" << offset[0] << ";" << offset[1] << "\t\t" << count[0] << ";" << count[1] << "\t\t" << status << endl;
			if ( status < 0 ) {
cout << "--->fail" << endl;
				read_status = false; break;
			}
			else {
cout << "--->dump" << endl;
				//data were read into rbuffer successfully, now interpret into pos structs and store in rawdata_pos
				size_t ReadBufferEnd = ( (AtomsToReadTotal - AtomsToReadCurr) > chunkdims[0] ) ? chunkdims[0] : (AtomsToReadTotal - AtomsToReadCurr);
				vector<pos>* storehere = rawdata_pos.at(chk);
				for ( size_t i = 0; i < ReadBufferEnd; ++i ) {
					storehere->push_back( pos( rbuf[4*i+0], rbuf[4*i+1], rbuf[4*i+2], rbuf[4*i+3] ) );
				}
				AtomsToReadCurr = AtomsToReadCurr + ReadBufferEnd;
cout << "--->dump AtomsToReadCurr\t\t" << AtomsToReadCurr << endl;

				//reset allocate buffer, ##MK::strictly speaking not necessary but enable to detect interpretation error
				for( size_t i = 0; i < elementsPerChunk; ++i )
					rbuf[i] = 0.f;
			}
			if ( read_status == true )
				continue;
			else { //a read error occurred, close file
				status = H5Sclose(mspcid);
				status = H5Dclose(dsetid);
				status = H5Sclose(fspcid);
				status = H5Fclose(fileid);
				break;
			}
		} //read and process next chunk
	}

	//read buffer exists
	delete [] rbuf; rbuf = NULL;

	//release HDF5 resources
	status = H5Sclose(mspcid);
	status = H5Dclose(dsetid);
	status = H5Sclose(fspcid);
	status = H5Fclose(fileid);

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	tictoc.prof_elpsdtime_and_mem( "HDF5FileInputIntoPOS", APT_IO, APT_IS_SEQ, mm, tic, toc);
	string mess = "Sequential HDF5-based data loading was successful and took " + to_string(toc-tic) + " seconds!";
	reporting( get_rank(), mess );

	return true;
}

#endif


bool solverHdl::compare_epos_pos( void )
{
	//##MK::debug compares if x,y,z, mq content of an epos and pos is exactly the same numerically in this case true otherwise returns false
	//length of both the same? EPOS and POS both store interleaved data per ion in the order of the evaporation sequence hence length should be the same
	size_t nep = 0;
	for (size_t ep = 0; ep < rawdata_epos.size(); ep++ ) { nep += rawdata_epos.at(ep)->size(); }
	size_t np = 0;
	for(size_t p = 0; p < rawdata_pos.size(); p++) { np += rawdata_pos.at(p)->size(); }
/*
 	if ( nep != np ) {
		cout << "Total ion count of epos and pos is different!" << endl;
		return false;
	}
*/

	//rawdata_epos and rawdata_pos containers store items of different length therefore
	vector<float> epos;
	for( size_t i = 0; i < rawdata_epos.size(); i++) {
		for( size_t j = 0; j < rawdata_epos.at(i)->size(); j++ ) {
			epos.push_back( rawdata_epos.at(i)->at(j).x );
			epos.push_back( rawdata_epos.at(i)->at(j).y );
			epos.push_back( rawdata_epos.at(i)->at(j).z );
			epos.push_back( rawdata_epos.at(i)->at(j).mq );
		}
	}

cout << "Stored epos values..." << endl;
cout.precision(32);
	//compare with pos content
	size_t k = 0;
	float max[4] = {0.f, 0.f, 0.f, 0.f};
    size_t outliner[4] = {0, 0, 0, 0};
	for( size_t i = 0; i < rawdata_pos.size(); i++) {
		for( size_t j = 0; j < rawdata_pos.at(i)->size(); j++ ) {
			if ( epos[4*k+0] != rawdata_pos.at(i)->at(j).x ) {
				float diff = fabs(epos[4*k+0]) - fabs(rawdata_pos.at(i)->at(j).x);
				max[0] = std::max(diff, max[0]);
                if ( diff >= 0.05 ) outliner[0] += 1;
//				cout << k << "\t\t" << diff << endl;
				//cout << "EPOS and POS inconsistent at item " << k << " for x\t\t" << epos[4*k+0] << "\t\t" << rawdata_pos.at(i)->at(j).x << endl; //return false;
			}
			else
				cout << k << endl;
			if ( epos[4*k+1] != rawdata_pos.at(i)->at(j).y ) {
				float diff = fabs(epos[4*k+1]) - fabs(rawdata_pos.at(i)->at(j).y);
				max[1] = std::max(diff, max[1]);
                if ( diff >= 0.05 ) outliner[1] += 1;
//				cout << k << "\t\t" << diff << endl;
				//cout << "EPOS and POS inconsistent at item " << k << " for y\t\t" << epos[4*k+1] << "\t\t" << rawdata_pos.at(i)->at(j).y << endl; //return false;
			}
			else
				cout << k << endl;
			if ( epos[4*k+2] != rawdata_pos.at(i)->at(j).z ) {
				float diff = fabs(epos[4*k+2]) - fabs(rawdata_pos.at(i)->at(j).z);
				max[2] = std::max(diff, max[2]);
                if ( diff >= 0.05 ) outliner[2] += 1;
//				cout << k << "\t\t" << diff << endl;
				//cout << "EPOS and POS inconsistent at item " << k << " for z\t\t" << epos[4*k+2] << "\t\t" << rawdata_pos.at(i)->at(j).z << endl; //return false;
			}
			else
				cout << k << endl;
			if ( epos[4*k+3] != rawdata_pos.at(i)->at(j).mq ) {
				float diff = fabs(epos[4*k+3]) - fabs(rawdata_pos.at(i)->at(j).mq);
				max[3] = std::max(diff, max[3]);
                if ( diff >= 0.05 ) outliner[3] += 1;
//				cout << k << "\t\t" << diff << endl;
				//cout << "EPOS and POS inconsistent at item " << k << " for mq\t\t" << epos[4*k+3] << "\t\t" << rawdata_pos.at(i)->at(j).mq << endl; //return false;
			}
			else
//				cout << k << endl;

//			cout << endl;
			k++;
			if ( k % 100000 == 0 )
				cout << k << endl;
		}
	}

	cout << "Maximum difference in x,y,z,mq, respectively..." << endl;
	cout << max[0] << endl;
	cout << max[1] << endl;
	cout << max[2] << endl;
	cout << max[3] << endl;
    cout << "How many ions with fabs(epos)-fabs(pos) difference >= 0.05 nm..." << endl;
    cout << outliner[0] << endl;
    cout << outliner[1] << endl;
    cout << outliner[2] << endl;
    cout << outliner[3] << endl;
    

	cout << "Last value in pos" << endl;
	cout << rawdata_pos.back()->back().x << endl;
	cout << rawdata_pos.back()->back().y << endl;
	cout << rawdata_pos.back()->back().z << endl;
	cout << rawdata_pos.back()->back().mq << endl;

	cout << "EPOS and POS content is numerically the same" << endl;
	return true;
}


bool solverHdl::identify_ions( void )
{
	double tic = MPI_Wtime();
	//ranges ions based on periodic system fed with range file, if rangefile load failed UNKNOWN type is assigned
	if ( mypse.read_rangefile2(Settings::RRNGFilenameIn) == true ) {
		size_t nb = 0;
		if ( Settings::SyntheticTips == true )
				nb = rawdata_alf.buckets.size();
		else { //epos, pos, and hdf5
			if ( Settings::InputFileformat == E_IO_EPOS )
				nb = rawdata_epos.size();
			else if ( Settings::InputFileformat == E_IO_POS )
				nb = rawdata_pos.size();
			else if ( Settings::InputFileformat == E_IO_HDF5 )
				nb = rawdata_pos.size();
			else
				nb = 0; //nothing to do
		}

		for ( size_t b = 0; b < nb; ++b ) {
			size_t n = 0; //do data in this bucket exist?
			if ( Settings::SyntheticTips == true ) {
				if ( rawdata_alf.buckets.at(b) != NULL )
					n = rawdata_alf.buckets.at(b)->size();
			}
			else {
				if ( Settings::InputFileformat == E_IO_EPOS ) {
					if ( rawdata_epos.at(b) != NULL )
						n = rawdata_epos.at(b)->size();
				}
				else { //Settings::InputFileformat == E_IO_POS ||  Settings::InputFileformat == E_IO_HDF5
					if ( rawdata_pos.at(b) != NULL )
						n = rawdata_pos.at(b)->size();
				}
			}

			vector<unsigned int>* wbucket = NULL;
			if ( n > 0) { //data exist
				try {
					wbucket = new vector<unsigned int>;
				}
				catch (bad_alloc &croak) {
					stopping( myRank, "Unable to allocate memory for iontype identification");
					return false; //allocated blocks will be cleared by solverHdl destructor
				}
				if ( Settings::SyntheticTips == true ) {
					vector<pos>* these = rawdata_alf.buckets.at(b);
					for(size_t i = 0; i < n; ++i)
						wbucket->push_back( mypse.mq2ionname(these->at(i).mq) );
				}
				else {
					if ( Settings::InputFileformat == E_IO_EPOS ) {
						vector<epos>* these = rawdata_epos.at(b);
						for(size_t i = 0; i < n; ++i)
							wbucket->push_back( mypse.mq2ionname(these->at(i).mq) );
					}
					else { //Settings::InputFileformat == E_IO_POS || Settings::InputFileformat == E_IO_HDF5
						vector<pos>* these = rawdata_pos.at(b);
						for(size_t i = 0; i < n; ++i)
							wbucket->push_back( mypse.mq2ionname(these->at(i).mq) );
					}
				} //copy-construct pointer to registration of identification IDs for all ions referenced by the rbucket
			}
			rawdata_iontype.push_back(wbucket);

		} //next bucket of ions

		cout << "Ranging performed, now summarizing..." << endl;

		//report result of the ranging process, i.e. how many ions of which type are unsigned int from [0,MaximumTypeID)
		vector<size_t> iontype_cnt = vector<size_t>( mypse.get_maxtypeid(), 0L );

		//Iontypes
		for(size_t i = 0; i < rawdata_iontype.size(); ++i) {
			vector<unsigned int>* these = rawdata_iontype.at(i);
			if ( these != NULL ) {
				size_t nj = these->size();
				for(size_t j = 0; j < nj; ++j) {
					unsigned int thisone = these->at(j);
					iontype_cnt.at( thisone ) += 1;
				}
			}
		}
		//report
		cout << "Result of the ranging process is as follows..." << endl;
		size_t total = 0;
		for(unsigned int id = 0; id < mypse.get_maxtypeid(); id++) {
			cout << mypse.ion_type2name( id ) << "\t\t\t" << iontype_cnt.at(id) << endl;
			total += iontype_cnt.at(id);
		}
		cout << "Ions in total\t\t" << total << endl;

		double toc = MPI_Wtime();
		memsnapshot mm = tictoc.get_memoryconsumption();
		tictoc.prof_elpsdtime_and_mem( "RRNGFileInputWithRanging", APT_RRR, APT_IS_SEQ, mm, tic, toc);
		return true;
	}
	return false;
}


bool solverHdl::generate_hdf5_resultsfile()
{
	//MK::initialize a simulation ID specific output file
	//MK::the key idea of PARAPROBE analyses is to have one set of measured file
	//rawdata RHIT/HITS, Cameca LEAP processed data APT, and user-postprocessed files bundling in HDF5
	//h5iometa info = h5iometa();
	int status = 0;
	string prefix = "PARAPROBE.SimID." + to_string(Settings::SimID);
	string fn_res = prefix + ".Results.h5";
	string fn_clu = prefix + ".Clustering.h5";
	string fn_cry = prefix + ".Crystallography.h5";

	if ( resultsh5Hdl.create_paraprobe_results_file( fn_res ) != WRAPPED_HDF5_SUCCESS )
		return false;
	if ( Settings::ClusteringTask != E_NOCLUST ) {
		if ( clusth5Hdl.create_paraprobe_clust_file( fn_clu ) != WRAPPED_HDF5_SUCCESS )
			return false;
	}

	if ( Settings::ExtractCrystallographicInfo != E_NOCRYSTALLO ) {
		if ( crysth5Hdl.create_paraprobe_cryst_file( fn_cry ) != WRAPPED_HDF5_SUCCESS )
			return false;
	}

	//fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".VoroTess.h5";
	//if ( voronoih5Hdl.create_paraprobe_voronoi_file( fn ) != WRAPPED_HDF5_SUCCESS )
	//	return false;

	return true;
}


void solverHdl::delete_rawdata( void )
{
	double tic = MPI_Wtime();

	size_t nb = 0;

	nb = rawdata_alf.buckets.size();
	for(size_t b = 0; b < nb; ++b ) {
		delete rawdata_alf.buckets.at(b);
		rawdata_alf.buckets.at(b) = NULL;
	}
	rawdata_alf.buckets.clear();

	nb = rawdata_epos.size();
	for (size_t b = 0; b < nb; ++b ) {
		delete rawdata_epos.at(b); //destructor of vector takes care that content of bucket b is deleted and memory freed
		rawdata_epos.at(b) = NULL;
	}
	rawdata_epos.clear();

	nb = rawdata_pos.size();
	for (size_t b = 0; b < nb; ++b ) {
		delete rawdata_pos.at(b); //destructor of vector takes care that content of bucket b is deleted and memory freed
		rawdata_pos.at(b) = NULL;
	}
	rawdata_pos.clear();

	nb = rawdata_iontype.size();
	for (size_t b = 0; b < nb; ++b ) {
		delete rawdata_iontype.at(b);
		rawdata_iontype.at(b) = NULL;
	}
	rawdata_iontype.clear();

	//##MK::do not reset to zero nevt = 0;

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	tictoc.prof_elpsdtime_and_mem( "RawdataDestruction", APT_UTL, APT_IS_SEQ, mm, tic, toc);
}


void solverHdl::initialize_thread_binding()
{
	//initialize NUMA binding
	double tic = MPI_Wtime();

	cout << "Initialize NUMA binding..." << "\n";
	vector<NUMANodeType> nodes;
	cout << "Bytesize of NUMANodeType " << sizeof(NUMANodeType) << "\n";
	int status = numa_available();
	if ( status < 0 ) {
		cout << "Your system does not support NUMA API\n";
		return;
	}
	cout << "Your system supports NUMA API\n";

	// returns a mask of CPUs on which the current task is allowed to run.
	bitmask* mask = numa_get_run_node_mask();
	bitmask* cpus = numa_allocate_cpumask();
	cout << "RunNodeMaskSize " << mask->size << endl;

	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			cout << "We are allowed to use node " << j << "\n";
			NUMANodeType node;
			memset(&node, 0xFF, sizeof(node));
			//converts a node number to a bitmask of CPUs.
			//The user must pass a bitmask structure with a mask buffer long enough to represent all possible cpu's
			numa_node_to_cpus(j, cpus);
			node.num_cpus = my_numa_bitmask_weight(cpus);
			cout << "node.num_cpus " << node.num_cpus << "\n";
			cout << "cpus size " << cpus->size << "\n";

			int cpuCounter = 0;
			for (unsigned int i = 0; i < cpus->size; i++) {
				if (numa_bitmask_isbitset(cpus, i)
						&& numa_bitmask_isbitset(numa_all_cpus_ptr, i)) {
					node.numa_cpus[cpuCounter] = i;
					cpuCounter++;
					cout << "\t\ti/cpuCounter " << i << ";" << cpuCounter << endl;
				}
			}
			nodes.push_back(node);
		}
		//else
		//	cout << "We are not allowed to use node " << j << "\n";
	}
	numa_free_cpumask(cpus);

	//we know a priori that the MAWS15,16,17 have 18 HT core pairs each on two nodes
	//first package IDs run from 0/36 to 2/38 and so forth ie even internal names
	//second package IDs run from 1/37 to 3/39 and so forth ie odd internal names
	//verify this using lstopo machine topology assessment
	vector<int> NextFreeCPUonNode;
	NextFreeCPUonNode.push_back(0); //first entry on node 0 e.g. P#1
	NextFreeCPUonNode.push_back(0);	//first entry on node 1 e.g. P#2

	#pragma omp parallel //everything shared by default
	{
		int threadID = omp_get_thread_num(); //threadlocal variables
		for( int tid = 0; tid < omp_get_num_threads(); tid++ ) { //enforced ordered placement of the threads
			if ( tid == threadID ) { //each thread only once
				#pragma omp critical
				{
					int UseThisNodeLocalCoreID = 666; //still unknown
					int UseThisNodeID = 666;

					//MAWS15,16,17
					if ( threadID % 2 == 0 )
						UseThisNodeID = 0; //map to first node package P#0
					else
						UseThisNodeID = 1; //map to second node package P#1

					//MAWS30
					/*
					if ( threadID < 10 )
						UseThisNodeID = 0; //map to first node package P#0
					else
						UseThisNodeID = 1; //map to second node package P#1
					*/

					UseThisNodeLocalCoreID = NextFreeCPUonNode.at(UseThisNodeID);
					NextFreeCPUonNode.at(UseThisNodeID)++;

					cout << "Will bind thread " << threadID << " to package node " << UseThisNodeID << " lstopo core id ";
					cout << nodes.at(UseThisNodeID).numa_cpus[UseThisNodeLocalCoreID];
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(UseThisNodeID).numa_cpus[UseThisNodeLocalCoreID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					if ( res == 0 )
						cout << " Done\n";
					else
						cout << " Failed\n";
				}
			}
			#pragma omp barrier
		}
	}

	double toc = MPI_Wtime();
	memsnapshot mm = memsnapshot();
	tictoc.prof_elpsdtime_and_mem( "InitNUMABinding", APT_UTL, APT_IS_SEQ, mm, tic, toc);
}

/*
unsigned int solverHdl::plan_parameter_study_counts( void )
{
	//double timer1 = MPI_Wtime();
	unsigned int jobid = 0;
	for ( apt_real eta = Settings::DetEffMin; eta <= Settings::DetEffMax; eta += Settings::DetEffIncr ) {
		for ( apt_real kf = Settings::KFMin; kf <= Settings::KFMax; kf += Settings::KFIncr ) {
			for ( apt_real icf = Settings::ICFMin; icf <= Settings::ICFMax; icf += Settings::ICFIncr ) {
				jobid++;
			}
		}
	}
	return jobid;
}


unsigned int solverHdl::jobid2rank( const unsigned int jid ) {
	unsigned int r = (1 + jid) % nRanks;
	//round robin for processes
	return r;
}


void solverHdl::plan_parameter_study( void )
{
	//MK::return 1 on success, zero otherwise
	double timer1 = MPI_Wtime();

	unsigned int jobid = 0;
	unsigned int targetrank = 0;
	unsigned int targetthread = 0;
	//unsigned int NOMP = MPI_COMM_WORLD_OMP_GET_NUM_THREADS; //##MK::make dynamic with call to OMP_thread_max or alike!

	if ( get_rank() == MASTER ) {
		std::cout << "jID\t\tMPIRnk\t\tOMPThr\t\tEta\t\tKF\t\tICF" << std::endl;
	}

	for ( apt_real eta = Settings::DetEffMin; eta <= Settings::DetEffMax; eta += Settings::DetEffIncr ) {
		for ( apt_real kf = Settings::KFMin; kf <= Settings::KFMax; kf += Settings::KFIncr ) {
			for ( apt_real icf = Settings::ICFMin; icf <= Settings::ICFMax; icf += Settings::ICFIncr ) {
				//construct parameter configuration
				targetrank = jobid2rank(jobid);
				targetthread = MASTER;
				
				configuration.push_back( runparm(eta, kf, icf, targetrank, targetthread, jobid ) );
				jobid++;

				if ( get_rank() == MASTER ) { 
					std::cout << jobid << "\t\t" << targetrank << "\t\t" << targetthread << "\t\t" << eta << "\t\t" << kf << "\t\t" << icf << std::endl;
				}
			}
		}
	}
	double timer2 = MPI_Wtime();
	std::cout << "...Worker " << get_rank() << " generated " << configuration.size() << " parameter configurations in " << (timer2-timer1) << " seconds." << std::endl;
}
*/


void solverHdl::set_mpidatatypes( void ) 
{
/*
	//MPI_Datatype MPI_EPOS_EventIO_Type; an example how to do a compound type
	int elementCounts[2] = {9, 2};
	MPI_Aint displacements[2] = {0, 9 * MPIIO_OFFSET_INCR_FLOAT32};
	MPI_Datatype oldTypes[2] = {MPI_FLOAT, MPI_INT};
	MPI_Type_create_struct(2, elementCounts, displacements, oldTypes, &MPI_EPOS_EventIO_Type);
	MPI_Type_commit(&MPI_EPOS_EventIO_Type);
	//std::cout << "...Worker " << this->get_rank() << " MPI_EPOS_EventIO_Type commited sizeof " << sizeof(MPI_EPOS_EventIO) << std::endl;
*/
}

