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


#ifndef __PARAPROBE_SETTINGS_H__
#define __PARAPROBE_SETTINGS_H__

#include "PARAPROBE_Math.h"


//add thirdparty XML reader header library functionality by M. Kalicinski
#include "thirdparty/RapidXML/rapidxml.hpp"

using namespace rapidxml;


enum E_INPUTFILEFORMAT {
	E_IO_NOTHING,
	E_IO_POS,
	E_IO_EPOS,
	E_IO_HDF5,
	E_IO_SYNTHETIC
};


enum E_ANALYSIS_MODE {
	E_ANALYSIS_DEFAULT
};

enum E_RECONSTRUCTION_ALGORITHM {
	E_RECON_NOTHING,				//nothing
	E_RECON_READFROMSYNTH,			//no reconstruction, take x,y,z values from synthetic dataset
	E_RECON_READFROMPOS,			//no reconstruction, take x,y,z values from pos assuming them to be in recon space...
	E_RECON_READFROMEPOS,			//no reconstruction, take x,y,z values from epos assuming them to be in recon space...
	E_RECON_BAS_ETAL				//use epos input and apply Bas et al. reconstruction algorithm from EPOS raw ion data
};


enum E_TIPSURFMESHING_ALGORITHM {
	E_NOSURFACE,
	E_ALPHASHAPE_CGAL,				//with which algorithm to triangularize the surface of the tip
	E_CONVEXHULL_CGAL,
	E_MARCHINGCUBE_IGL,
	E_LOADEXISTENTVTK
};


enum E_DISTANCE_METRICS {
	E_NOSPATSTAT,					//nothing
	E_RDF,							//radial distribution function
	E_NEAREST_NEIGHBOR,				//nearest neighbor
	E_RIPLEYK,						//Ripley K
	E_KNN,							//k-nearest if existent in SpatStatMax sphere
	E_MKNN,							//collection of k-nearest if existent individually in SpatStatMax sphere
	E_MULTISPATSTAT					//multiple analyses
};


enum E_CLUSTERING {
	E_NOCLUST,						//nothing
	E_DBSCAN,						//DBScan
	E_MAXSEPARATION,				//MaximumSeparationMethod
	E_ISOSURF						//Isosurface-based
};


class Settings {
public:
	
	static E_INPUTFILEFORMAT InputFileformat;
	static E_ANALYSIS_MODE AnalysisMode;				//what at all to work on?
	static E_RECONSTRUCTION_ALGORITHM ReconstructionAlgo;
	static E_TIPSURFMESHING_ALGORITHM SurfaceMeshingAlgo;
	static E_DISTANCE_METRICS SpatialDistributionTask;
	static E_CLUSTERING ClusteringTask;
	static string RAWFilenameIn;
	static string RRNGFilenameIn;
	static string SurfaceFilenameIn;
	static string DescrStatTasksCode;
	static string ClusteringTasksCode;
	static string DescrStatMKNNCode;

	static apt_real AdvIonPruneBinWidthMin;
	static apt_real AdvIonPruneBinWidthIncr;
	static apt_real AdvIonPruneBinWidthMax;
	static apt_real FlightLength;
	static apt_real AtomicDensity;
	static apt_real EvaporationField;
	static apt_real DetEffMin;
	static apt_real DetEffIncr;
	static apt_real DetEffMax;
	static apt_real KFMin;
	static apt_real KFIncr;
	static apt_real KFMax;
	static apt_real ICFMin;
	static apt_real ICFIncr;
	static apt_real ICFMax;
	static apt_real SpatStatRadiusMin;
	static apt_real SpatStatRadiusIncr;
	static apt_real SpatStatRadiusMax;
	static apt_real ClustMSDmaxMin;
	static apt_real ClustMSDmaxIncr;
	static apt_real ClustMSDmaxMax;
	//static apt_real ClustMSDdilation;
	//static apt_real ClustMSDerosion;
	static apt_real RelBottomRadius;
	static apt_real RelTopRadius;
	static apt_real RelBottomCapHeight;
	static apt_real RelTopCapHeight;
	static apt_real LatticeConstant;
	static apt_real SimDetEfficiency;			//to account for non-ideal detection efficiency
	static apt_real SpatResolutionSigmaX;
	static apt_real SpatResolutionSigmaY;
	static apt_real SpatResolutionSigmaZ;		//to account for finite APT resolution
	static apt_real ClusterRadiusMean;
	static apt_real ClusterRadiusSigmaSqr;

	static size_t SpatStatKNNOrder;
	static size_t ClustMSNmin;
	static size_t NumberOfAtoms;
	static size_t NumberOfCluster;

	static unsigned int SimID;							//fixed identifier of this particular analysis run
	static unsigned int RndDescrStatsPRNGSeed;
	static unsigned int RndDescrStatsPRNGDiscard;

	static bool IdentifyIonType;
	static bool IOReconstruction;
	static bool IOTriangulation;
	static bool IOHKFilteredIons;
	static bool IORAWHKClusterID;
	static bool IOIonTipSurfDists;
	static bool SyntheticTips;
	static bool DebugComputeDistance;			//bringing unbiased distributions

	static bool SpatStatDoRDF;
	static bool SpatStatDo1NN;
	static bool SpatStatDoRIPK;
	static bool SpatStatDoKNN;
	static bool SpatStatDoMKNN;
	static bool SpatStatAddLabelRnd;
	static bool ClustPostSpatStat;


//prototypes
	static void readXML(string filename = "");
	static bool checkUserInput( void );
	
	//apt_real read_real( xml_node<>* const src, const string keyword, const apt_real defaultvalue );
};

#endif
