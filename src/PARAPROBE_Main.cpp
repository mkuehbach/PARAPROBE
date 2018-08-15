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

//parameter handshake
#define SIMID									1
#define CONTROLFILE								2


void helloworld ( int pargc, char** pargv )
{
	cout << "Starting up PARAPROBE v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION;
	if ( VERSION_BETASTAGE == 1 ) 	cout << " beta stage" << endl;
	else 						cout << endl;

	cout << "ERRORs are terminating in every case..." << endl;
	if ( pargc < 2 ) {
		std::cout << "\t\tERROR::We need at least a simulation id <unsigned int> and an XML control file <*.xml> before doing something useful!" << endl;
		return;
	}
}


bool init(  int pargc, char** pargv )
{
	Settings::SimID = stoul( pargv[SIMID] );
	try {
		Settings::readXML(pargv[CONTROLFILE]);
	}
	catch (exception& e) { 
		cout << endl << "\t\tERROR::Unable to parse control file! Details:\n" << e.what() << endl; return false;
	}
	if ( Settings::checkUserInput() == false ) {
		cout << endl << "\t\tERROR::Control file settings failed the validity check!" << endl; return false;
	}
	else {
		cout << endl << "\t\tInput is valid under SimulationID = " << "SimID." <<  Settings::SimID << endl;
	}
	cout << "All console prompts which follow are intended for debugging purposes only..." << endl;
	cout << endl << endl;

	return true;
}


//genering global functions to report state, warnings and erros
void reporting( const int rank, const string what ) {
	cout << "VERBOSE::" << rank << " " << what << endl;
}
void reporting( const string what ) {
	cout << "VERBOSE::" << what << endl;
}


void complaining( const int rank, const string what ) {
	cout << "WARNING::" << rank << " " << what << endl;
}
void complaining( const string what ) {
	cout << "WARNING::" << what << endl;
}


void stopping( const int rank, const string what ) {
	cout << "ERROR::" << rank << " " << what << endl;
}
void stopping( const string what ) {
	cout << "ERROR::" << what << endl;
}


int main(int argc, char** argv)
{    
//DEBUG HDF
//	debug_hdf5();
//	return 0;

//SETUP PROGRAM AND PARAMETER BUT DO NOT YET LOAD MEASUREMENT DATA
	helloworld( argc, argv );
	if ( init( argc, argv ) == false ) {
		return 0;
	}
	
//go MPI process parallel with hybrid OpenMP threading capability, funneled means only main thread will make MPI calls
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread( &argc, &argv, supportlevel_desired, &supportlevel_provided);

	double gtic = MPI_Wtime();
	//by now we are parallel...
	int nr = 1;
	int r = MASTER;
	if ( supportlevel_provided < supportlevel_desired ) {
		stopping( r, "Insufficient threading capabilities of the MPI library!");
		MPI_Finalize(); //required because threading insufficiency does not imply process is incapable to work at all
		return 0;
	}
	else { 
		MPI_Comm_size(MPI_COMM_WORLD, &nr);
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
	}
	reporting( r, "-th MPI process initialized, we are now MPI_COMM_WORLD parallel using MPI_THREAD_FUNNELED");

//INITIALIZE MPI SOLVER INSTANCE AND LOAD DATASET
	//generate MPI rank solverHdl instance managing which of the parameter configuration rank 0 does
	int localhealth = 1;
	solverHdl* hdl = NULL;
	try { hdl = new solverHdl; }
	catch (bad_alloc &exc) {
		localhealth = 0;
		stopping( r, "Unable to allocate a solverInstance" );
	}
	//were all processes, who should have, able to generate a process-local solver class object instance?
	int globalhealth = 0;
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if ( globalhealth != nr ) {
		stopping( r, "Not all processes were able to allocate a solverInstance");
		MPI_Finalize();
		return 0;
	}

//all have a solver instance, well then set it up
	hdl->set_rank(r);
	hdl->set_nranks(nr);
	hdl->set_mpidatatypes();

//load APT measurement file by all in parallel, individually identify ions with rangefile
	localhealth = 1; //set to zero if process is get sick because of lacking data or errors
	if ( Settings::InputFileformat == E_IO_SYNTHETIC ) {
    	if ( hdl->generate_synthetic_tip() == true ) {
    		reporting( r, "A synthetic tip was successfully generated");
    	}
    	else {
    		complaining(r, "Unable to synthesize a tip"); localhealth = 0;
    	}
    }
    else {
    	if ( Settings::InputFileformat == E_IO_POS &&
    		Settings::RAWFilenameIn.substr(Settings::RAWFilenameIn.length()-4) == ".pos" ) {
    		if ( hdl->load_pos_sequentially( Settings::RAWFilenameIn ) == true ) {
    			reporting( r, "POS successfully read sequentially...");
    		}
    		else {
    			complaining( r, "Unable to completely process *.pos file"); localhealth = 0;
    		}
    	}
    	else if ( Settings::InputFileformat == E_IO_EPOS &&
    		Settings::RAWFilenameIn.substr(Settings::RAWFilenameIn.length()-5) == ".epos" ) {
    		if ( hdl->load_epos_sequentially( Settings::RAWFilenameIn ) == true ) {
    			reporting( r, "EPOS successfully read sequentially...");
    		}
    		else {
    			complaining( r, "Unable to completely process *.epos file"); localhealth = 0;
    		}
    	}
#ifdef UTILIZE_HDF5
    	else if ( Settings::InputFileformat == E_IO_HDF5 &&
    		Settings::RAWFilenameIn.substr(Settings::RAWFilenameIn.length()-3) == ".h5" ) {
    		if ( hdl->load_hdf5_sequentially( Settings::RAWFilenameIn ) == true ) {
    			reporting( r, "HDF5 successfully read sequentially...");
    		}
    		else {
    			complaining( r, "Unable to completely process *.h5 file"); localhealth = 0;
    		}
    	}
#endif
    	else {
    		complaining( r, "Invalid input file format chosen"); localhealth = 0;
    	}
    }

//identify ions aka ranging
	if ( localhealth == 1 ) { //healthy i.e. files were loaded successfully
		if ( hdl->identify_ions() == true )
			reporting( r, "Ions were identified successfully...");
		else
			complaining( r, "Ions were identified by default type only. Maybe that is not what was intended...");
	}
	//a barrier is necessary such that everybody has the input and is assured properly setup
	//##MK::having rank MASTER read only while others processes wait and broadcast is very likely premature optimization...
	MPI_Barrier(MPI_COMM_WORLD);
	globalhealth = 0;
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if ( globalhealth != nr ) {
		stopping( r, "Not all processes were able to allocate a solverInstance");
		MPI_Finalize();
		return 0;
	}

//generate a solver to perform reconstruction and analysis task, which is equipped with a reconClass and higherOrderHdl
	localhealth = 1;
	solver* sr = NULL;
	try {
		sr = new solver;
		sr->recon->owner = sr;
		sr->sp->owner = sr;
		sr->binner->owner = sr;
		sr->surf->owner = sr;
		sr->rndmizer->owner = sr;
		sr->hometrics->owner = sr;
		sr->cldetect->owner = sr;
	}
	catch (bad_alloc &exc) {
		complaining( r, "Unable to allocate a solver");
		localhealth = 0;
	}

	if ( localhealth == 1 ) { //solver has been initialized
		sr->owner = hdl;

//generate a reconstruction to transform rawdata_pos/_epos into x,y,z point process in every case
		if ( Settings::ReconstructionAlgo != E_RECON_NOTHING  ) {
			sr->volume_reconstruction();

			//MK::from now on the rawdata on hdl are no longer necessary and could already be deleted to safe memory inplace
			hdl->delete_rawdata();

			sr->spatial_decomposition();

//additional datamining in recon space with this reconstruction?
			if ( Settings::AnalysisMode == E_ANALYSIS_DEFAULT ) {

				sr->volume_binning();

				sr->surface_triangulation();
				sr->surface_distancing2();

				sr->characterize_distances();

				sr->init_spatialindexing();

				if ( Settings::SpatialDistributionTask != E_NOSPATSTAT )
					sr->characterize_spatstat();

				if ( Settings::ClusteringTask != E_NOCLUST )
					sr->characterize_clustering();
			}
		}
	}
//report profiling
	hdl->tictoc.spit_profiling( Settings::SimID, hdl->get_rank() );

//delete solver object if existent
	delete sr; sr = NULL;

	double gtoc = MPI_Wtime();
	string mess = "Elapsed time on process " + to_string(hdl->get_rank()) + " " + to_string((gtoc-gtic)) + " seconds";
	reporting( mess );

//delete solverHdl, deconstruct parallel processes, and exit
	delete hdl; hdl = NULL;

//better have all processes joining to exit program cooperatively
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


//##MK::TO DO
//boundary contact analysis


//##MK::TO DO
//##MK::find isolated functions of member functions with void argument
//##MK::change HoshenKopelman to threading and from unsigned int to size_t
//##MK::potentially still a flaw in the HoshenKopelman part
//##MK::add maximum number of ion < UINT32MAX-1 check
//##MK::input check... in Settings


//##MK::DEBUG overload to check numerical representation of content inside EPOS and POS
/*
	hdl->load_epos_sequentially( "R76_30139-v02.epos" );
	hdl->load_pos_sequentially( "R76_30139-v02.pos" );
	hdl->compare_epos_pos();
	cout << "Controlled ending for debugging!" << endl;
	localhealth = 0;
*/
