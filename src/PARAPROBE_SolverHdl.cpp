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

	//bvh = NULL;
}


solver::~solver()
{
	//MK::do not clear owner as it is only a backreference to my owner solverHdl

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

		apt_xyz RSQR = SQR(Settings::SpatStatRadiusMax); //##MK::max of spat and clustering
		sp->db.at(mt)->ion2surfdistance_init( RSQR );

		#pragma omp critical
		{
			cout << "Thread " << mt << " initialized distance field for " << sp->db.at(mt)->ionpp3.size() << " ions" << endl;
		}
	}

	double toc = MPI_Wtime();
	owner->tictoc.prof( "Ion2SurfDistInit", APT_GEO, tic, toc);
	cout << "TipSurface distance fields initialized in " << (toc-tic) << " seconds" << endl;

	if ( hullexists == false ) { //done and out
		return;
	}

	tic = MPI_Wtime();

	//computing distances desired and therefore building threadglobal Rtree
	surf->build_rtree();

	toc = MPI_Wtime();
	owner->tictoc.prof( "SurfTriangleHullRTreeConstruct", APT_BVH, tic, toc);
	cout << "Building triangle tree BVH " << (toc-tic) << " seconds" << endl;

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
		const apt_xyz R = Settings::SpatStatRadiusMax; //##MK::change at some point to max(spatstat and clustering)
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
	surf->chop_rtree();

	toc = MPI_Wtime();
	owner->tictoc.prof( "Ion2SurfThreadparallelDistancing", APT_GEO, tic, toc);
	cout << "Computing parallelized distance to surface in " << (toc-tic) << " seconds" << endl;
}


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
		owner->tictoc.prof( "VTKFileOutputIon2Surf", APT_IO, tic, toc);
	}

	cout << "Ion distance to surface distribution characterized " << endl;
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
	owner->tictoc.prof( "ThreadparallelKDTreeConstruct", APT_BVH, tic, toc);

	tic = MPI_Wtime();

	if ( AllKDTreesValid == true ) {
		sp->kdtree_success = true;
		reporting( owner->get_rank(), "KDTree construction was globally successful!");

		sp->reportpartitioning();
	}
	else {
		sp->kdtree_success = false;
		stopping( owner->get_rank(), "KDTree construction was not successful globally!");
	}

	toc = MPI_Wtime();
	owner->tictoc.prof( "ThreadparallelKDTreeReportResults", APT_IO, tic, toc);
	cout << "Threadlocal KDTrees constructed in parallel in " << (toc-tic) << " seconds" << endl;
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
		hometrics->compute_generic_spatstat2();

		if ( Settings::SpatStatAddLabelRnd == true ) {
			//MK::now randomize iontype assignment over entire tip to get randomize spatial statistics
			rndmizer->reset();
			rndmizer->learn_original_types_global();
			rndmizer->shuffle_types_mt19937();
			rndmizer->apply_shuffled_types_global();

			//so characterize again but this time with randomization
			hometrics->compute_generic_spatstat2();

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

/*
bool reconstructor::reconstruction_accept_synthetic( void )
{
	double tic = MPI_Wtime();

	//feed through rawdata_pos/_epos, where are they?
	solverHdl* here = owner->owner;

	//MK::if process pragma omp parallel here potential to comply with first touch policy that OpenMP threads get portions of the rawdata in threadlocal memory
	//MK::HOWEVER, at this stage the rawdata in the solverHdl object are not necessarily reconstructed yet, therefore we have only a very approximate
	//MK>>idea of a good spatial partitioning of these, in particular if we seek to reconstruct no all ions, therefore, we better build the OpenMP
	//thread/local stuff after the reconstruction...

	vector<p3d>* wpbucket = NULL;
	vector<unsigned int>* wlbucket = NULL;
	try {
		wpbucket = new vector<p3d>;
		wlbucket = new vector<unsigned int>;
	}
	catch (bad_alloc &reconexc) {
		complaining( "Unable to allocate memory in reconstructor for storing synthetic ion locations");
		return false;
	}

	size_t nb = here->rawdata_alf.buckets.size();
	for ( size_t b = 0; b < nb; ++b ) {
		vector<p3dm1>* rbucketsynth = NULL;
		rbucketsynth = here->rawdata_alf.buckets.at(b);
		if ( rbucketsynth != NULL ) {
			for ( size_t i = 0; i < rbucketsynth->size(); ++i ) {
				wpbucket->push_back( p3d(
						rbucketsynth->at(i).x,
						rbucketsynth->at(i).y,
						rbucketsynth->at(i).z ) );
				wlbucket->push_back( rbucketsynth->at(i).m );
			}
		}
	}

	//register in pp3
	pp3.push_back(NULL);
	pp3.back() = wpbucket;
	lbls.push_back(NULL);
	lbls.back() = wlbucket;

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof( "SynthReconstruction", APT_REC, tic, toc);
	reporting( "Reconstruction synthesized in " + to_string(toc-tic) + " seconds");

	if ( Settings::IOReconstruction == true ) {
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".SyntheticRecon.vtk";
		reconstruction_vtk( pp3, lbls, runparm( 0.f, 0.f, 0.f, 0, 0, 0), fn );

		toc = MPI_Wtime();
		owner->owner->tictoc.prof( "VTKFileOutputReconstruction", APT_IO, tic, toc);
	}

	return true;
}
*/


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

		pp3.push_back(NULL);
		lbls.push_back(NULL);

		if ( n > 0) { //data exist
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
			pp3.back() = wpbucket;
			lbls.back() = wlbucket;
		}
	} //process next bin

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof( "AcceptReconstruction", APT_REC, tic, toc);
	reporting( "Reconstruction accepted in " + to_string(toc-tic) + " seconds");

	if ( Settings::IOReconstruction == true ) {
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".AcceptedRecon.vtk";
		reconstruction_vtk( pp3, lbls, runparm( 0.f, 0.f, 0.f, 0, 0, 0), fn );

		toc = MPI_Wtime();
		owner->owner->tictoc.prof( "VTKFileOutputReconstruction", APT_IO, tic, toc);
	}

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
	owner->owner->tictoc.prof( "BasEtAlReconstruction", APT_REC, tic, toc);
	string mess = "Reconstruction including memory setup took " +  to_string(toc - tic) + " seconds!";
	reporting( here->get_rank(), mess);

	if ( Settings::IOReconstruction == true ) {
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".BasEtAlRecon.vtk";
		reconstruction_vtk( pp3, lbls, runparm( ETA, KF, ICF, 0, 0, 0), fn );

		toc = MPI_Wtime();
		owner->owner->tictoc.prof( "VTKFileOutputReconstruction", APT_REC, tic, toc);
	}

	return true;
}


threadmemory::threadmemory()
{
	owner = NULL;
	threadtree = NULL;
	zmi = F32MX;
	zmi = F32MI;
	melast = false;
}


threadmemory::~threadmemory()
{
	//MK::do not delete owner only backreference to decompositor who owns me!
	if ( threadtree != NULL ) {
		delete threadtree; threadtree = NULL;
	}
}


bool threadmemory::init_localmemory( const apt_xyz zmin, const apt_xyz zmax, const bool mlast )
{
	//MK::CALLED FROM WITHIN PARALLEL REGION
	zmi = zmin;
	zmx = zmax;
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

	for(size_t b = 0; b < nb; ++b) {
		vector<p3d>* thesepoints = pointshere->pp3.at(b);
		vector<unsigned int>* theselabels = labelshere->lbls.at(b);
		if ( thesepoints != NULL && theselabels != NULL ) {
			size_t pn = thesepoints->size();
			size_t ln = theselabels->size();
			if ( pn == ln ) {
				for(size_t i = 0; i < pn; ++i) { //z < zmii ? not included : (z < zmxx) ? included : not included
					if ( thesepoints->at(i).z < zmi )
						continue;

					if ( thesepoints->at(i).z < zmx )
						ionpp3.push_back( p3dm1(thesepoints->at(i).x,  thesepoints->at(i).y,  thesepoints->at(i).z, theselabels->at(i)) );
					else
						continue;
				}
			}
		}
		else if ( thesepoints == NULL && theselabels == NULL ) {
			//cout << "--->Bucket " << b << " empty" << endl;
			//stopping("Inaccessible pn and ln in read thread-local geometry and label bucket information");
			//return false;
		}
		else {
			stopping("Inaccessible pn or ln in read thread-local geometry and label bucket information");
			return false;
		}
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

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof( "IonPositionExtrema", APT_UTL, tic, toc);
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
	vector<apt_xyz> zcoordinates;

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
	//sort once ascendingly
	sort(zcoordinates.begin(), zcoordinates.end() );

	//when there is a stack of N domains we split N-1 times
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

	//build threadlocal partitioning objects
	for(size_t mt = MASTER; mt < PARAPROBE_NUM_THREADS; mt++) {
		db.push_back(NULL);
	}


	#pragma omp parallel shared(qq,qqss, zcoordinates) //shared but only read
	{
		unsigned int nt = static_cast<unsigned int>(omp_get_num_threads()); //##MK::what is the default behavior of parallel region start with all?
		unsigned int mt = static_cast<unsigned int>(omp_get_thread_num());

		apt_real zmii = (mt == MASTER) ? zcoordinates.at(0) : qq.at(mt-1);
		apt_real zmxx = (mt == (nt-1)) ? (zcoordinates.back() + EPSILON) : qq.at(mt);
		bool last = (mt == (nt-1)) ? true : false;

		#pragma omp critical
		{
			cout << "Thread " << omp_get_thread_num() << " working on [" <<  zmii << ";" << zmxx << ") is last? " << last << endl;
		}

		//MK::points belong to mt if their z position in [zmii, zmxx) for which we probe as follows left boundary z < zmii ? not included : (z < zmxx) ? included : not included

		threadmemory* memlocal = NULL;
		try {
			memlocal = new threadmemory; //allocate surplus first-touch
			memlocal->owner = this;

			memlocal->init_localmemory( zmii, zmxx, last );

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
	owner->owner->tictoc.prof( "Loadpartitioning", APT_BVH, tic, toc);
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
	owner->owner->tictoc.prof( "VTKFileInputSurfTrianguleHull", APT_IO, tic, toc);
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

	del_pruning_mem();

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof( "PruningIonsSurfaceTriangulation", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "SurfTriangulationDelaunayConstruct", APT_GEO, tic, toc);
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

	//set shape to automatically determined optimum
	//as.set_alpha(*opt);
	as.set_alpha(alpha_solid);
	cout << "Taking the smallest alpha value to get a solid through data points is " << alpha_solid << " for triangulation" << endl;
	assert(as.number_of_solid_components() == 1); //##MK::

	toc = MPI_Wtime();
	owner->owner->tictoc.prof( "SurfTriangulationAlphaShapeConstruct", APT_GEO, tic, toc);
	cout << "CGAL::Delaunay-based alpha shape generated and assessed in " << (toc-tic) << endl;

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
	owner->owner->tictoc.prof( "SurfTriangulationAlphaShapeTriExtract", APT_GEO, tic, toc);
	cout << "Extracted triangularized tip surface (" << nf << "/" << tri.size() << " triangles in total) in " << (toc-tic) << " seconds" << endl;

	if ( Settings::IOTriangulation == true ) {
		tic = MPI_Wtime();

		string fn = "PARAPROBE.SimID." + to_string(Settings::SimID) + ".TipSurface.vtk";
		triangulation_vtk_naive( tipsurface, fn );

		toc = MPI_Wtime();
		owner->owner->tictoc.prof( "VTKFileOutputSurfaceTriangleHull", APT_IO, tic, toc);
	}

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

	return true;
}


void surfacer::chop_rtree()
{
	delete bvh; bvh = NULL;
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
	owner->owner->tictoc.prof( "RndLabelingLearnOldLabelsGlobal", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "RndLabelingShuffleTypesMT19937Global", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "RndLabelingApplyNewLabelsGlobal", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "RndLabelingResetOldLabelsGlobal", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "DescrStatsInitTasks", APT_UTL, tic, toc);
}


/*
void horderdist::compute_generic_spatstat1()
{
	//get all neighbors within Settings::SpatStatRadius
	//each thread machine off his ions binning in local summary statistics which is afterwards dumped to MASTER thread i critical region
	//threads may require to probe KDTrees of their bottom and top neighbors --- potentially multiple
	//(for very flat z slabs---many threads---) can do so in parallel because KDTrees are queried only in local function call but not written to
	double tic = MPI_Wtime();

	cout << "Threadparallel general higher order spatial distribution function..." << endl;
	vector<histogram> globalres; //multi-layer results for the individual tasks
	for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
		globalres.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );

	#pragma omp parallel
	{
		double mytic, mytoc;
		mytic = omp_get_wtime();
		int jnt = omp_get_num_threads();
		int jmt = omp_get_thread_num();
		unsigned int mt = static_cast<unsigned int>(jmt);

		//all distances squared for efficiency unless nbor objects
		//basic operation is as follows take each ion, if it is sufficiently distant from tip take into consideration
		decompositor* partitioning = owner->sp;
		threadmemory* mydata = partitioning->db.at(mt);
		vector<p3dm1> const & theseions = mydata->ionpp3_kdtree;
		vector<apt_xyz> & thesedistances = mydata->ion2surf_kdtree;
		size_t ni = theseions.size();
		apt_xyz R = Settings::SpatStatRadiusMax;
		apt_xyz RSQR = SQR(R);
		kd_tree* kauri = mydata->threadtree;

		apt_real mydata_zmi = mydata->get_zmi();
		apt_real mydata_zmx = mydata->get_zmx();

		//set up thread-local task-local histograms collecting to improve parallel efficiency and avoid critical regions, these results are fused later into globalres
		vector<histogram> myres;
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			myres.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );

		size_t myworkload = 0;
		size_t mylost = 0;

		for(size_t i = 0; i < ni; ++i) {
			p3dm1 me = theseions[i];

			//pre-screening, as we use one global function to extract this is inefficient when probing only tasks in which
			//alloying elements are to be studied, in this case it is much more efficient to screen first when the target ion is at
			//of type in at least one specific task
			bool considerme = false;
			for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
				if ( me.m == spatstat_tasks.at(tsk).target.second ) {
					considerme = true;
					break;
				}
			}

			if ( thesedistances.at(i) >= RSQR && considerme == true ) { //ion should be considered deep enough inside tip volume
				vector<nbor> neighbors;
				neighbors.clear();

				if ( (me.z - R) > mydata_zmi && (me.z + R) < mydata_zmx ) { //we have to probe only my local tree, most likely case, moderate overhead
					kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors );
				}
				else { //we have to probe my local tree and potentially trees of neighboring threads
					kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors );

					//probe potentially top neighbors
					if ( jmt < (jnt-1) ) { //when i am not the topthread already, in which case I would have checked my ions already, i do climb up
						for( int nb = (jmt+1); nb < jnt; nb++) { //will eventually climb up to the topmost threadregion, also in practice this will never happen

							threadmemory* nbordata = partitioning->db.at(nb);
							vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
							kd_tree* nborkauri = nbordata->threadtree;

							nborkauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors );
						}
					}

					//probe potentially bottom neighbors
					if ( jmt > MASTER ) { //when i am not the bottomthread already, also in which case i would have checked my ions already, i do climb down
						for( int nb = (jmt-1); nb > -1; nb-- ) {
							threadmemory* nbordata = partitioning->db.at(nb);
							vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
							kd_tree* nborkauri = nbordata->threadtree;

							nborkauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors );
						}
					}
				}

				//now that we know all neighbors do something useful with it, sort is not required
				myworkload += neighbors.size();

				if ( Settings::SpatialDistributionTask == E_RIPLEYK || Settings::SpatialDistributionTask == E_RDF ) {
					for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
						if ( me.m == spatstat_tasks.at(tsk).target.second ) {
							//consider me central ion only if of type required by that specific task tsk
							for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
								for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin();
										jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
									if ( neighbors.at(nbb).m != jt->second ) { //type of the only one and candidate is different
										continue;
									}
									else {
										//myres.at(tsk).add_nodump( neighbors.at(nbb).d );
										myres.at(tsk).add( neighbors.at(nbb).d );
										break;
									}
								} //done checking all keywords for tsk
							} //done checking all neighboring ions of me for that specific task
						}
					} //reutilize the extracted spatial environment of me again to get histogram of other tasks
				}
				else if ( Settings::SpatialDistributionTask == E_NEAREST_NEIGHBOR ) {
					for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
						if ( me.m == spatstat_tasks.at(tsk).target.second ) {
							apt_real closest = F32MX;
							size_t nborid = numeric_limits<size_t>::max();

							//consider me central ion only if of type required by that specific task tsk
							for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
								for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
									if ( neighbors.at(nbb).m != jt->second ) { //type of the only one and candidate is different
										continue;
									}
									else {
										if ( neighbors.at(nbb).d > closest ) { //most likely most ions in Settings::SpatStatRadiusMax
											continue;
										} else {
											closest = neighbors.at(nbb).d;
											nborid = nbb;
										}
									}
								} //done checking all keywords for tsk
							} //done checking all neighboring ions of me for that specific task

							if ( nborid != numeric_limits<size_t>::max() ) { //myres.at(tsk).add_nodump( neighbors.at(nbb).d );
								myres.at(tsk).add( closest );
							}
							else {
								mylost++;
							}
						}
					} //reutilize the extracted spatial environment of me again to get also histogram of other tasks
				}
				else if ( Settings::SpatialDistributionTask == E_KNN ) {
					//if at least k elements exist at all do n_th element partial sorting for distances
					for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
						if ( me.m == spatstat_tasks.at(tsk).target.second ) {
							//are there at least neighbors first of all regardless of type
							//to be able to name at all a k-th order neighbor?
							//1NN --> order = 0 if size() == 1 ie 1 > 0 we can feed
							//2NN --> order = 1 if size() == 2 ie 2 > 1 we can feed
							if ( neighbors.size() > Settings::SpatStatKNNOrder ) { //MK::even getting the zeroth element requires neighbors to contain at least one element!
								//filter target envcandidates before sorting
								vector<nbor> tmp;
								for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
									for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
										if ( neighbors.at(nbb).m != jt->second ) { //type of the only one and desired candidate is different
											continue;
										}
										else { tmp.push_back( neighbors.at(nbb) ); }
									} //done checking all keywords for tsk
								}
								//having now our envcandidates we make a partial sort to get the n_th element
								//still enough envcandidates to name a KNNOrder?
								if ( tmp.size() > Settings::SpatStatKNNOrder ) { //okay there is no way around O(n) at least partial sorting the envcandidates

									nth_element( tmp.begin(), tmp.begin() + Settings::SpatStatKNNOrder,
											tmp.end(), SortNeighborsForAscDistance );

									nbor luckyone = tmp.at(Settings::SpatStatKNNOrder);
									myres.at(tsk).add( luckyone.d );
								}
								else { mylost++; } //done not enough of specific type at all to give k-th
							}
							else { mylost++; } //done not enough in search sphere of any type at all to give k-th
						} //reutilize the extracted spatial environment of me again to get also histogram of other tasks
						//else return nothing
					}
				}
				else { //all other specific spatial statistics tasks
					continue;
				}
//cout << i << "\t\t" << neighbors.size() << endl;
			} //done checking a single ion
//cout << jmt << "/" << i << endl;
		} //next ion

		//done checking all my ions

		mytoc = omp_get_wtime();

		#pragma omp critical
		{
			//we a accumulating on the master thread
			for( size_t tsk = 0; tsk < myres.size(); ++tsk) {
				globalres.at(tsk).cnts_lowest += myres.at(tsk).cnts_lowest;
				for( size_t b = 0; b < myres.at(tsk).bincount(); ++b) {
					globalres.at(tsk).cnts.at(b) += myres.at(tsk).cnts.at(b);
				}
				globalres.at(tsk).cnts_highest += myres.at(tsk).cnts_highest;
			}

			cout << "Thread " << omp_get_thread_num() << " finished general spatstat " << myworkload << "/" << mylost << " took " << (mytoc-mytic) << " seconds" << endl;
		}

	} //end of parallel region

	//report results
	for(size_t tsk = 0; tsk < globalres.size(); ++tsk) {
		string what = "";
		if ( Settings::SpatialDistributionTask == E_RDF )	what = "RDF";
		else if ( Settings::SpatialDistributionTask == E_NEAREST_NEIGHBOR ) what = "NN";
		else if ( Settings::SpatialDistributionTask == E_RIPLEYK ) what = "RIPK";
		else if ( Settings::SpatialDistributionTask == E_KNN ) what = to_string(Settings::SpatStatKNNOrder) + "NN";
		else break;

		string whichtarget = spatstat_tasks.at(tsk).target.first;
		string whichcand = "";
		for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++) {
			whichcand += jt->first;
		}
		report_apriori_descrstat1( what, whichtarget, whichcand, globalres.at(tsk) );
	}

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof( "DescrStatsThreadparallelCompute", APT_PPP, tic, toc);
	cout << "Computing general spatial statistics completed took " << (toc-tic) << " seconds" << endl;
}
*/

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

	if ( Settings::SpatStatDoRDF == false && Settings::SpatStatDo1NN == false && Settings::SpatStatDoRIPK == false && Settings::SpatStatDoKNN == false && Settings::SpatStatDoMKNN == false ) {
		complaining( "No spatial statistics task to do" );
		return;
	}

	double tic = MPI_Wtime();
	cout << "Threadparallel general higher order spatial distribution function..." << endl;

	vector<histogram> g_res_rdf; //##MK::all threads know Settings so the implicit order of the histograms in globalres is well defined
	vector<histogram> g_res_1nn;
	vector<histogram> g_res_rpk;
	vector<histogram> g_res_knn;

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
	if ( Settings::SpatStatDoKNN == true )
		for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
			g_res_knn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
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
		vector<histogram> m_res_knn;
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
		if ( Settings::SpatStatDoKNN == true )
			for( size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk )
				m_res_knn.push_back( histogram(Settings::SpatStatRadiusMin, Settings::SpatStatRadiusIncr, Settings::SpatStatRadiusMax) );
		if ( Settings::SpatStatDoMKNN == true ) {
			//make threadlocal copy of, ##MK::potentially threadlocal copies of task list might also be beneficial to increase locality
			for( size_t i = 0; i < DescrStatMKNNCandidates.size(); ++i) {
				m_mknn_cand.push_back( DescrStatMKNNCandidates.at(i) );
			}

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

				if ( thesedistances.at(i) >= RSQR && considerme == true ) {
					MyIonsCurrRegionConsider++;

					//probe local environment up to R, get all neighbors, only for 1NN this is inefficient, for RDF and high k kNN it is required
					//and considering that once computed the neighbors are reutilized for all spatstat tasks and multiple distribution functions this is a superior strategy
					vector<nbor> neighbors; neighbors.clear();
					if ( (me.z - R) > currdata_zmi && (me.z + R) < currdata_zmx ) { //we have to probe only in curr_kauri this is the most likely case
						curr_kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors );
					}
					else {
						//we have to probe the curr_kauri and trees of neighboring regions, potentially multiple to the top and bottom
						curr_kauri->range_rball_noclear_nosort( i, theseions, RSQR, neighbors );
						//probe top neighbors and climb up, when i am not the topmost thread
						if ( thr < (jnt-1) ) {
							for( int nb = (thr+1); nb < jnt; nb++) { //will eventually climb up to the topmost threadregion, also in practice this will never happen
								threadmemory* nbordata = owner->sp->db.at(nb);
								vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
								kd_tree* nbor_kauri = nbordata->threadtree;
								nbor_kauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors );
							}
						}
						//probe bottom neighbors and climb down, when i am not the bottommost thread already
						if ( thr > MASTER ) {
							for( int nb = (thr-1); nb > -1; nb-- ) {
								threadmemory* nbordata = owner->sp->db.at(nb);
								vector<p3dm1> const & nborthreadions = nbordata->ionpp3_kdtree;
								kd_tree* nbor_kauri = nbordata->threadtree;
								nbor_kauri->range_rball_noclear_nosort_external( me, nborthreadions, RSQR, neighbors );
							}
						}
					}

					//use the environment for all descriptive statistics tasks and tasks of this ion
					if ( Settings::SpatStatDoRDF == true ) {
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) { //consider me central ion only if of type required by that specific task tsk but only account for once
									for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) { //scan cache locality efficiently threadlocal neighbor array for envcandidates
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors.at(nbb).m != jt->second ) { continue; }
											else { m_res_rdf.at(tsk).add( neighbors.at(nbb).d ); break; }
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
									for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors.at(nbb).m != jt->second ) { continue; }
											else {
												if ( neighbors.at(nbb).d > closest ) { continue; }//most likely most ions in Settings::SpatStatRadiusMax
												else { closest = neighbors.at(nbb).d; nborid = nbb; }
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
									for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) { //scan cache locality efficiently threadlocal neighbor array for envcandidates
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors.at(nbb).m != jt->second ) { continue; }
											else { m_res_rpk.at(tsk).add( neighbors.at(nbb).d ); break; }
										} //done checking all keywords for tsk
									} //done checking all neighboring ions of me for that specific task
									break;
								} //done performing task tsk
							}
						} //reutilize the extracted spatial environment of me again to get histogram of other tasks
					}
					if ( Settings::SpatStatDoKNN == true ) { //if at least k elements exist at all do n_th element partial sorting for distances
						for(size_t tsk = 0; tsk < spatstat_tasks.size(); ++tsk) {
							for( auto kt = spatstat_tasks.at(tsk).trgcandidates.begin(); kt != spatstat_tasks.at(tsk).trgcandidates.end(); ++kt ) {
								if ( me.m == kt->second ) {
									//are there at all neighbors
									//to be able to name at all a k-th order neighbor?
									//1NN --> order = 0 if size() == 1 ie 1 > 0 we can feed
									//2NN --> order = 1 if size() == 2 ie 2 > 1 we can feed
									if ( neighbors.size() > Settings::SpatStatKNNOrder ) { //MK::even getting the zeroth element requires neighbors to contain at least one element!
										//filter target envcandidates before sorting
										vector<nbor> tmp;
										for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
											for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
												if ( neighbors.at(nbb).m != jt->second ) { continue; }//type of the only one and desired candidate is different
												else { tmp.push_back( neighbors.at(nbb) ); }
											} //done checking all keywords for tsk
										}
										//having now our envcandidates we make a partial sort to get the n_th element
										//still enough envcandidates to name a KNNOrder?
										if ( tmp.size() > Settings::SpatStatKNNOrder ) { //okay there is no way around O(n) at least partial sorting the envcandidates
											nth_element( tmp.begin(), tmp.begin() + Settings::SpatStatKNNOrder, tmp.end(), SortNeighborsForAscDistance );
											nbor luckyone = tmp.at(Settings::SpatStatKNNOrder);
											m_res_knn.at(tsk).add( luckyone.d );
										}
										else
											m_res_knn.at(tsk).add( R + EPSILON );
									}
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
									for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
										for ( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); jt++ ) {
											if ( neighbors.at(nbb).m != jt->second ) { continue; }//type of the only one and desired candidate is different
											else { tmp.push_back( neighbors.at(nbb) ); }
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
				for( size_t tsk = 0; tsk < m_res_rdf.size(); ++tsk) {
					g_res_rdf.at(tsk).cnts_lowest += m_res_rdf.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_rdf.at(tsk).bincount(); ++b)
						g_res_rdf.at(tsk).cnts.at(b) += m_res_rdf.at(tsk).cnts.at(b);
					g_res_rdf.at(tsk).cnts_highest += m_res_rdf.at(tsk).cnts_highest;
				}
			}
			if ( Settings::SpatStatDo1NN == true ) {
				for( size_t tsk = 0; tsk < m_res_1nn.size(); ++tsk) {
					g_res_1nn.at(tsk).cnts_lowest += m_res_1nn.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_1nn.at(tsk).bincount(); ++b)
						g_res_1nn.at(tsk).cnts.at(b) += m_res_1nn.at(tsk).cnts.at(b);
					g_res_1nn.at(tsk).cnts_highest += m_res_1nn.at(tsk).cnts_highest;
				}
			}
			if ( Settings::SpatStatDoRIPK == true ) {
				for( size_t tsk = 0; tsk < m_res_rpk.size(); ++tsk) {
					g_res_rpk.at(tsk).cnts_lowest += m_res_rpk.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_rpk.at(tsk).bincount(); ++b)
						g_res_rpk.at(tsk).cnts.at(b) += m_res_rpk.at(tsk).cnts.at(b);
					g_res_rpk.at(tsk).cnts_highest += m_res_rpk.at(tsk).cnts_highest;
				}
			}
			if ( Settings::SpatStatDoKNN == true ) {
				for( size_t tsk = 0; tsk < m_res_knn.size(); ++tsk) {
					g_res_knn.at(tsk).cnts_lowest += m_res_knn.at(tsk).cnts_lowest;
					for( size_t b = 0; b < m_res_knn.at(tsk).bincount(); ++b)
						g_res_knn.at(tsk).cnts.at(b) += m_res_knn.at(tsk).cnts.at(b);
					g_res_knn.at(tsk).cnts_highest += m_res_knn.at(tsk).cnts_highest;
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
		}
	}
	//explicit barrier end of parallel region

	double toc = MPI_Wtime();
	owner->owner->tictoc.prof( "DescrStatsThreadparallelCompute", APT_PPP, tic, toc);
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
			report_apriori_descrstat2( 1, what, whichtarg, whichcand, g_res_rdf.at(tsk) ); //##MK::more elegant cast from enum to long
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
			report_apriori_descrstat2( 2, what, whichtarg, whichcand, g_res_1nn.at(tsk) ); //##MK
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
			report_apriori_descrstat2( 3, what, whichtarg, whichcand, g_res_rpk.at(tsk) ); //##MK
		}
	}
	if ( Settings::SpatStatDoKNN == true ) {
		for(size_t tsk = 0; tsk < g_res_knn.size(); ++tsk) {
			string what = "kNN";
			string whichtarg = ""; //spatstat_tasks.at(tsk).target.first;
			for( auto jt = spatstat_tasks.at(tsk).trgcandidates.begin(); jt != spatstat_tasks.at(tsk).trgcandidates.end(); ++jt)
				whichtarg += jt->first;
			string whichcand = "";
			for( auto jt = spatstat_tasks.at(tsk).envcandidates.begin(); jt != spatstat_tasks.at(tsk).envcandidates.end(); ++jt)
				whichcand += jt->first;
			report_apriori_descrstat2( 4, what, whichtarg, whichcand, g_res_knn.at(tsk) ); //##MK
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

				report_apriori_descrstat2( 5, what, whichtarg, whichcand, g_res_mknn.at(thishist) ); //g_res_mknn.at(tsk).at(kth) );
			}
		}
	}

	toc = MPI_Wtime();
	owner->owner->tictoc.prof( "ReportingDescrStatsResults", APT_IO, tic, toc);
	cout << "Reporting general spatial statistics completed took " << (toc-tic) << " seconds" << endl;
}




void horderdist::report_apriori_descrstat1( const string whichmetric,
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
	apt_real bend = hist.start() + hist.width();

	if ( Settings::SpatialDistributionTask != E_RDF ) {	//mode-specific file layout
		sslog << "BinEnd(r);Counts\n";
		sslog << "nm;1\n";
		sslog << "BinEnd(r);Counts\n"; //we report binends and accumulated counts

		//below hist.start() lower tail dump
		sslog << hist.start() << ";" << hist.cnts_lowest << "\n";
		for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
			sslog << bend << ";" << hist.report(b) << "\n";
		}
		//above hist.end() upper tail dump
		sslog << "" << ";" << hist.cnts_highest << "\n";
	}
	else { //== E_RDF requires additional post-processing
		sslog << "NumberOfIonsInside" << ";" << owner->binner->metadata.nions_inside << "\n";
		sslog << "VolumeInside" << ";" << owner->binner->metadata.volume_inside << " (nm^3)\n";

		if ( owner->binner->metadata.volume_inside > EPSILON ) { //now division is safnm^3...

			apt_real globaldensity = static_cast<apt_real>(owner->binner->metadata.nions_inside);
			globaldensity /= static_cast<apt_real>(owner->binner->metadata.volume_inside);
			sslog << "InsideAtomicDensityEst" << ";" << globaldensity << "\n";

			apt_real norm = (1.f / globaldensity) * (1.f / ((4.f/3.f)*PI));
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

			//according to A. Baddeley

			apt_real half_dr = 0.5*hist.width();

			for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
				apt_real r =  hist.start() + (0.5 + static_cast<apt_real>(b)) * hist.width();
				apt_real diff = (b > 0) ? static_cast<apt_real>(hist.report(b)-hist.report(b-1)) : 0.f;
				apt_real sphvol = CUBE(r + half_dr) - CUBE(r - half_dr);
				apt_real rdfval = (sphvol > EPSILON) ? norm*diff/sphvol : 0.f;

				sslog << r << ";" << hist.report(b) << ";" << diff << ";" << sphvol << ";" << rdfval << "\n";
			}
			//everything above hist.end() is also not of interest to us just report upper tail dump
			sslog << "" << ";" << hist.cnts_highest << ";;;\n";

		}
		else {
			complaining("Insufficient volume detected during binning in order to estimate RDF!");
		}
	}

	sslog.flush();
	sslog.close();
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
	apt_real bend = hist.start() + hist.width();

	//global meta information

	//##MK::add global status information here

	if ( tsktype == E_RDF ) { //mode-specific file layout
		sslog << "NumberOfIonsInside" << ";" << owner->binner->metadata.nions_inside << "\n";
		sslog << "VolumeInside" << ";" << owner->binner->metadata.volume_inside << " (nm^3)\n";

		if ( owner->binner->metadata.volume_inside > EPSILON ) { //now division is safnm^3...

			apt_real globaldensity = static_cast<apt_real>(owner->binner->metadata.nions_inside);
			globaldensity /= static_cast<apt_real>(owner->binner->metadata.volume_inside);
			sslog << "InsideAtomicDensityEst" << ";" << globaldensity << "\n";

			apt_real norm = (1.f / globaldensity) * (1.f / ((4.f/3.f)*PI));
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

			apt_real half_dr = 0.5*hist.width();

			for( unsigned int b = 0; b < hist.bincount(); ++b, bend += hist.width() ) {
				apt_real r =  hist.start() + (0.5 + static_cast<apt_real>(b)) * hist.width();
				apt_real diff = (b > 0) ? static_cast<apt_real>(hist.report(b)-hist.report(b-1)) : 0.f;
				apt_real sphvol = CUBE(r + half_dr) - CUBE(r - half_dr);
				apt_real rdfval = (sphvol > EPSILON) ? norm*diff/sphvol : 0.f;

				sslog << r << ";" << hist.report(b) << ";" << diff << ";" << sphvol << ";" << rdfval << "\n";
			}
			//everything above hist.end() is also not of interest to us just report upper tail dump
			sslog << "" << ";" << hist.cnts_highest << ";;;\n";
		}
	}
	if ( tsktype == E_NEAREST_NEIGHBOR || tsktype == E_KNN || tsktype == E_MKNN ) {
		if ( tsktype == E_NEAREST_NEIGHBOR )	sslog << "1NN\n";
		if ( tsktype == E_KNN )					sslog << to_string(Settings::SpatStatKNNOrder+1) << "NN\n";
		if ( tsktype == E_MKNN )				sslog << whichmetric << "\n";

		sslog << "BinEnd(r);Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n";
		sslog << "nm;1;1;1\n";
		sslog << "BinEnd(r);Counts;ECDFon[BinMinBinMax];CDFon[BinMinBinMax]\n"; //we report binends and accumulated counts

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
	}
	if ( tsktype == E_RIPLEYK ) {
		sslog << "BinEnd(r);AccumulatedCounts\n";
		sslog << "nm;1\n";
		sslog << "BinEnd(r);AccumulatedCounts\n"; //we report binends and accumulated counts

		//below hist.start() lower tail dump
		sslog << hist.start() << ";" << hist.cnts_lowest << "\n";
		//on[ ) interval
		size_t sum = 0; //MK::must be larger enough
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
										tskid, runid );
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
		int jnt = omp_get_num_threads(); unsigned int nt = static_cast<unsigned int>(jnt);
		int jmt = omp_get_thread_num(); unsigned int mt = static_cast<unsigned int>(jmt);

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
				else if ( Settings::SpatialDistributionTask == E_KNN ) {
					//are there at least neighbors first of all regardless of type
					//to be able to name at all a k-th order neighbor?
					if ( neighbors.size() > Settings::SpatStatKNNOrder ) { //MK::even getting the zeroth element requires neighbors to contain at least one element!
						//filter target candidates before sorting
						vector<nbor> tmp;
						for( size_t nbb = 0; nbb < neighbors.size(); ++nbb) {
							if( is_clustered( neighbors.at(nbb).m, mxtypid ) == true )
								continue;
							else
								tmp.push_back( neighbors.at(nbb) );
						}
						//still enough candidates to name a KNNOrder?
						if ( tmp.size() > Settings::SpatStatKNNOrder ) {

							nth_element( tmp.begin(), tmp.begin() + Settings::SpatStatKNNOrder,
									tmp.end(), SortNeighborsForAscDistance );

							nbor luckyone = tmp.at(Settings::SpatStatKNNOrder);
							myres.add( luckyone.d );
						}
						else { mylost++; } //done not enough of specific type at all to give k-th
					}
					else { mylost++; } //done not enough in search sphere of any type at all to give k-th
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
	else if ( Settings::SpatialDistributionTask == E_KNN ) what = to_string(Settings::SpatStatKNNOrder) + "NN";
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
	owner->owner->tictoc.prof( "TipBinning", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "ClusteringInitTasks", APT_UTL, tic, toc);
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
	owner->owner->tictoc.prof( "ClusteringMaximumSeparationDmaxStudy", APT_CLU, tic, toc);
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

	sslog << "DMaxID;DMax;NMin;Ncluster;Nlinkpoints;Ncorepoints;Nnoisepoints;Nclusteredpoints\n";
	sslog << "1;nm;1;1;1;1;1;1\n";
	sslog << "DMaxID;DMax;NMin;Ncluster;Nlinkpoints;Ncorepoints;Nnoisepoints;Nclusteredpoints\n";
	size_t id = 0;
	for( auto it = results.begin(); it != results.end(); ++it, id++ ) {
		sslog << id << ";" << it->Dmax << ";" << it->Nmin << ";" << it->nClusterFound << ";";
		sslog << (it->nClustered - it->nCore) << ";" << it->nCore << ";" << it->nNoise << ";" << it->nClustered << "\n";
	}
	sslog.flush();
	sslog.close();

	toc = MPI_Wtime();
	cout << "Reporting MaximumSeparation dmax range scan results took " << (toc-tic) << " seconds" << endl;
}



solverHdl::solverHdl()
{
	nevt = 0;

	myRank = MASTER;
	nRanks = SINGLE_PROCESS;
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

						/*if ( placed * total < progress )
							continue;
						else {
							cout << placed * total << " % placed" << endl;
							progress += 1.f;
						}*/
					}
				}
			}
		}
	}

	reporting( "All tip atoms were placed");

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
	tictoc.prof( "GenerateSyntheticTip", APT_UTL, tic, toc);
	string mess = "Sequential tip synthesis was successful and took " + to_string(toc-tic) + " seconds!";
	reporting( get_rank(), mess );


#ifdef UTILIZE_HDF5
	if ( Settings::IOReconstruction == true ) { //store generated tip into an HDF5 file
		tic = MPI_Wtime();

		write_pos_hdf5( rawdata_alf.buckets, Settings::RAWFilenameIn );

		toc = MPI_Wtime();
		tictoc.prof( "WriteSyntheticTipHDF", APT_IO, tic, toc);
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
	tictoc.prof( "POSFileInputWithEndiannessSwop", APT_IO, tic, toc);
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
	tictoc.prof( "EPOSFileInputWithEndiannessSwop", APT_IO, tic, toc);
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
	int rank = H5Sget_simple_extent_ndims(fspcid);
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
	tictoc.prof( "HDF5FileInputIntoPOS", APT_IO, tic, toc);
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
	if ( mypse.read_rangefile(Settings::RRNGFilenameIn) == true ) {
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
		vector<size_t> iontype_cnt;
		for(unsigned int i = 0; i < mypse.get_maxtypeid(); i++)
			iontype_cnt.push_back( 0L );

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
		tictoc.prof( "RRNGFileInputWithRanging", APT_RRR, tic, toc);
		return true;
	}
	return false;
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
	tictoc.prof( "RawdataDestruction", APT_UTL, tic, toc);
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
	tictoc.prof( "InitNUMABinding", APT_UTL, tic, toc);
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

