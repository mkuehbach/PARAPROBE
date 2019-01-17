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

#include "PARAPROBE_Settings.h"

//using namespace rapidxml;

//##MK::init dummy
#define REALDUMMY				EPSILON;
#define UCHARDUMMY				(0x00);


E_INPUTFILEFORMAT Settings::InputFileformat = E_IO_NOTHING;
string Settings::RAWFilenameIn = "";
bool Settings::SyntheticTips = false;

E_ANALYSIS_MODE Settings::AnalysisMode = E_ANALYZE_NOTHING;
E_RECONSTRUCTION_ALGORITHM Settings::ReconstructionAlgo = E_RECON_NOTHING;
bool Settings::IdentifyIonType = false;
string Settings::RRNGFilenameIn = "";

E_TIPSURFMESHING_ALGORITHM Settings::SurfaceMeshingAlgo = E_NOSURFACE;
E_ALPHASHAPE_ALPHAVALUE_CHOICE Settings::SurfaceAShapeAValue = E_ASHAPE_SMALLEST_SOLID;
E_VOLUME_TESSELLATION Settings::VolumeTessellation = E_NOTESS;
string Settings::SurfaceFilenameIn = "";

E_CRYSTALLOGRAPHY Settings::ExtractCrystallographicInfo = E_NOCRYSTALLO;
E_WINDOWING Settings::WindowingMethod = E_RECTANGULAR_WINDOW;

E_DISTANCE_METRICS Settings::SpatialDistributionTask = E_NOSPATSTAT;

E_CLUSTERING Settings::ClusteringTask = E_NOCLUST;
E_NUMABINDING Settings::NumaBinding = E_NOBINDING;

unsigned int Settings::SimID = 0;

//detector space and signal mapping parameter
/*apt_real Settings::DetBinWidth = 1.f; //nanometer
apt_real Settings::SDMBinWidth = 0.01;*/

//reconstruction parameter
apt_real Settings::FlightLength = 80.0; //best practice 80-100 nm
apt_real Settings::AtomicDensity = EPSILON; //Aluminium
apt_real Settings::EvaporationField = EPSILON;
apt_real Settings::DetEffMin = 1.f;
apt_real Settings::DetEffIncr = 1.f;
apt_real Settings::DetEffMax = 1.f;
apt_real Settings::KFMin = 2.0;
apt_real Settings::KFIncr = 0.2;
apt_real Settings::KFMax = 3.2;
apt_real Settings::ICFMin = 1.65;
apt_real Settings::ICFIncr = 0.1;
apt_real Settings::ICFMax = 1.65;

//tip boundary/surface reconstruction
apt_real Settings::AdvIonPruneBinWidthMin = 1.f; //1nm best practice guess useful for many cases
apt_real Settings::AdvIonPruneBinWidthIncr = 1.f;
apt_real Settings::AdvIonPruneBinWidthMax = 1.f;
bool Settings::DebugComputeDistance = true; //best practice should be to get unbiased distributions

//extraction of crystallographic information
apt_real Settings::TessellationGuardWidth = 5.f; 	//so many average position distances to include in guard needs to be larger than 1.f + EPSILON however do not change unless you have better estimate of a priori guard demands
apt_real Settings::TessellationPointsPerBlock = 5;	//performance compromise spatial querying Voro++ needs to be at least 1.f + EPSILON will be interpreted into integer
apt_real Settings::CrystalloRadiusMax = 0.f;
apt_real Settings::SamplingGridBinWidthX = 0.f;
apt_real Settings::SamplingGridBinWidthY = 0.f;
apt_real Settings::SamplingGridBinWidthZ = 0.f;
apt_real Settings::ElevationAngleMin = 0.f;
apt_real Settings::ElevationAngleIncr = 0.f;
apt_real Settings::ElevationAngleMax = 0.f;
apt_real Settings::AzimuthAngleMin = 0.f;
apt_real Settings::AzimuthAngleIncr = 0.f;
apt_real Settings::AzimuthAngleMax = 0.f;
apt_real Settings::WindowingAlpha = 0.f;
int Settings::CrystalloHistoM = 9;


//descriptive spatial statistics
string Settings::DescrStatTasksCode = "";
apt_real Settings::SpatStatRadiusMin = 0.f; //[0:0.1:2.0] interval best practice
apt_real Settings::SpatStatRadiusIncr = 0.1;
apt_real Settings::SpatStatRadiusMax = 2.f;
size_t Settings::SpatStatKNNOrder = 1-1; //default first order neighbor only, i.e. Cstyle zero-th...
bool Settings::SpatStatDoRDF = false;
bool Settings::SpatStatDo1NN = false;
bool Settings::SpatStatDoRIPK = false;
bool Settings::SpatStatDoMKNN = false;
bool Settings::SpatStatDoNPCorr = false;
bool Settings::SpatStatDoCountNbors = false;
bool Settings::SpatStatAddLabelRnd = false;
string Settings::DescrStatMKNNCode = "";

//clustering analysis
string Settings::ClusteringTasksCode = "";
apt_real Settings::ClustMSDmaxMin = 0.3;
apt_real Settings::ClustMSDmaxIncr = 0.05;
apt_real Settings::ClustMSDmaxMax = 0.3;
size_t Settings::ClustMSCntMin = 1;
size_t Settings::ClustMSCntIncr = 1;
size_t Settings::ClustMSCntMax = 1;
size_t Settings::ClustMSNmin = 1;
bool Settings::ClustPostSpatStat = false;
//apt_real Settings::ClustMSDdilation = DBL_DUMMY;
//apt_real Settings::ClustMSDerosion = DBL_DUMMY;
apt_real Settings::SurfaceCellsCarvingRadius = 2.f; //2nm default


//plotting and I/O options
bool Settings::IOReconstruction = false;
bool Settings::IOTriangulation = false;
bool Settings::IOTriangulationBVH = false;
bool Settings::IOKDTreePartitioning = false;
bool Settings::IOHKFilteredIons = false;
bool Settings::IORAWHKClusterID = false;
bool Settings::IOIonTipSurfDists = false;
bool Settings::IOVoronoiDescrStats = false;
bool Settings::IOVoronoiCellPositions = false;
bool Settings::IOVoronoiTopoGeom = false;
bool Settings::IOCrystallography = false;


//tip synthesis
apt_real Settings::RelBottomRadius = 0.10;
apt_real Settings::RelTopRadius = 0.05;
apt_real Settings::RelBottomCapHeight = 0.05;
apt_real Settings::RelTopCapHeight = 0.05;
apt_real Settings::LatticeConstant = 0.404; //Al lattice in nm
size_t Settings::NumberOfAtoms = 0;
apt_real Settings::SimDetEfficiency = 1.f;
apt_real Settings::SpatResolutionSigmaX = 0.f; //nm
apt_real Settings::SpatResolutionSigmaY = 0.f; //nm
apt_real Settings::SpatResolutionSigmaZ = 0.f; //nm
size_t Settings::NumberOfCluster = 0;
apt_real Settings::ClusterRadiusMean = EPSILON; //nm
apt_real Settings::ClusterRadiusSigmaSqr = EPSILON; //nm


//predefined values
unsigned int Settings::RndDescrStatsPRNGSeed = -1;
unsigned int Settings::RndDescrStatsPRNGDiscard = 700000;


inline apt_real str2real( const string str )
{
#ifdef EMPLOY_SINGLEPRECISION
	return stof( str );
#else
	return stod( str );
#endif
}



void Settings::readXML(string filename) {
	//find the desired .xml file
	if ( 0 == filename.compare("") )
		filename = string("PARAPROBE.Input.Debug.xml");

	ifstream file( filename );
	if ( file.fail() ) 
		throw runtime_error(string("Unable to locate input file ") + filename);

	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		throw runtime_error("Undefined parameters file!");
	}

	unsigned int arg = 0;
	unsigned int mode = 0;

	//WHAT TYPE OF ANALYSIS TO DO
	mode = E_IO_NOTHING;
	if (0 != rootNode->first_node("InputFileformat"))
		mode = stoul( rootNode->first_node("InputFileformat")->value());
	switch(mode)
	{
		case E_IO_POS:
			InputFileformat = E_IO_POS; break;
		case E_IO_EPOS:
			InputFileformat = E_IO_EPOS; break;
		case E_IO_HDF5:
			InputFileformat = E_IO_HDF5; break;
		case E_IO_SYNTHETIC:
			InputFileformat = E_IO_SYNTHETIC; break;
		default:
			InputFileformat = E_IO_NOTHING;
	}
	SyntheticTips = (InputFileformat == E_IO_SYNTHETIC) ? true : false;
	if (0 != rootNode->first_node("RAWFilenameIn"))
		RAWFilenameIn = rootNode->first_node("RAWFilenameIn")->value();

	mode = E_ANALYZE_NOTHING;
	if (0 != rootNode->first_node("AnalysisMode"))
		mode = stoul( rootNode->first_node("AnalysisMode")->value());
	switch (mode)
	{
		case E_ANALYZE_IN_RECONSTRUCTION_SPACE:
			AnalysisMode = E_ANALYZE_IN_RECONSTRUCTION_SPACE; break;
		default:
			AnalysisMode = E_ANALYZE_NOTHING;
	}

	mode = E_RECON_NOTHING;
	if (0 != rootNode->first_node("ReconstructionAlgorithm"))
		mode = stoul( rootNode->first_node("ReconstructionAlgorithm")->value());
	switch (mode)
	{
		case E_RECON_NOTHING:
			ReconstructionAlgo = E_RECON_NOTHING; break;
		case E_RECON_READFROMSYNTH:
			ReconstructionAlgo = E_RECON_READFROMSYNTH; break;
		case E_RECON_READFROMPOS:
			ReconstructionAlgo = E_RECON_READFROMPOS; break;
		case E_RECON_READFROMEPOS:
			ReconstructionAlgo = E_RECON_READFROMEPOS; break;
		case E_RECON_BAS_ETAL:
			ReconstructionAlgo = E_RECON_BAS_ETAL; break;
		default:
			ReconstructionAlgo = E_RECON_BAS_ETAL;
	}

	IdentifyIonType = false;
	if (0 != rootNode->first_node("IdentifyIonType")) {
		arg = stoul( rootNode->first_node("IdentifyIonType")->value() );
		if ( arg != 0 )
			IdentifyIonType = true;
	}
	if (0 != rootNode->first_node("RRNGFilenameIn"))
		RRNGFilenameIn = rootNode->first_node("RRNGFilenameIn")->value();

	mode = E_NOSURFACE;
	if (0 != rootNode->first_node("SurfaceReconstructionType"))
		mode = stoul( rootNode->first_node("SurfaceReconstructionType")->value());
	switch (mode)
	{
		case E_NOSURFACE:
			SurfaceMeshingAlgo = E_NOSURFACE; break;
		case E_ALPHASHAPE_CGAL:
			SurfaceMeshingAlgo = E_ALPHASHAPE_CGAL; break;
		case E_CONVEXHULL_CGAL:
			SurfaceMeshingAlgo = E_CONVEXHULL_CGAL; break;
		case E_MARCHINGCUBE_IGL:
			SurfaceMeshingAlgo = E_MARCHINGCUBE_IGL; break;
		case E_LOADEXISTENTVTK:
			SurfaceMeshingAlgo = E_LOADEXISTENTVTK; break;
		default:
			SurfaceMeshingAlgo = E_NOSURFACE;
	}
	if (0 != rootNode->first_node("SurfaceFilenameIn"))
		SurfaceFilenameIn = rootNode->first_node("SurfaceFilenameIn")->value();
	if ( Settings::SurfaceMeshingAlgo == E_ALPHASHAPE_CGAL ) {
		mode = E_ASHAPE_SMALLEST_SOLID;
		if (0 != rootNode->first_node("AlphaShapeAlphaValueChoice"))
			mode = stoul( rootNode->first_node("AlphaShapeAlphaValueChoice")->value());
		switch (mode)
		{
			case E_ASHAPE_SMALLEST_SOLID:
				SurfaceAShapeAValue = E_ASHAPE_SMALLEST_SOLID; break;
			case E_ASHAPE_CGAL_OPTIMAL:
				SurfaceAShapeAValue = E_ASHAPE_CGAL_OPTIMAL; break;
			default:
				SurfaceAShapeAValue = E_ASHAPE_SMALLEST_SOLID;
		}
	}

	mode = E_NOCRYSTALLO;
	if (0 != rootNode->first_node("AnalysisCrystallographicInfo"))
		mode = stoul( rootNode->first_node("AnalysisCrystallographicInfo")->value() );
	switch(mode)
	{
		case E_ARAULLO_PETERS_METHOD:
			ExtractCrystallographicInfo = E_ARAULLO_PETERS_METHOD; break;
		default:
			ExtractCrystallographicInfo = E_NOCRYSTALLO;
	}
	if ( ExtractCrystallographicInfo == E_ARAULLO_PETERS_METHOD ) {
		if (0 != rootNode->first_node("CrystalloRadiusMax"))
			CrystalloRadiusMax = str2real( rootNode->first_node("CrystalloRadiusMax")->value() );
		if (0 != rootNode->first_node("SamplingGridBinWidthX"))
			SamplingGridBinWidthX = str2real( rootNode->first_node("SamplingGridBinWidthX")->value() );
		if (0 != rootNode->first_node("SamplingGridBinWidthY"))
			SamplingGridBinWidthY = str2real( rootNode->first_node("SamplingGridBinWidthY")->value() );
		if (0 != rootNode->first_node("SamplingGridBinWidthZ"))
			SamplingGridBinWidthZ = str2real( rootNode->first_node("SamplingGridBinWidthZ")->value() );
		if (0 != rootNode->first_node("ElevationAngleMin"))
			ElevationAngleMin = DEGREE2RADIANT(str2real( rootNode->first_node("ElevationAngleMin")->value() ));
		if (0 != rootNode->first_node("ElevationAngleIncr"))
			ElevationAngleIncr = DEGREE2RADIANT(str2real( rootNode->first_node("ElevationAngleIncr")->value() ));
		if (0 != rootNode->first_node("ElevationAngleMax"))
			ElevationAngleMax = DEGREE2RADIANT(str2real( rootNode->first_node("ElevationAngleMax")->value() ));
		if (0 != rootNode->first_node("AzimuthAngleMin"))
			AzimuthAngleMin = DEGREE2RADIANT(str2real( rootNode->first_node("AzimuthAngleMin")->value() ));
		if (0 != rootNode->first_node("AzimuthAngleIncr"))
			AzimuthAngleIncr = DEGREE2RADIANT(str2real( rootNode->first_node("AzimuthAngleIncr")->value() ));
		if (0 != rootNode->first_node("AzimuthAngleMax"))
			AzimuthAngleMax = DEGREE2RADIANT(str2real( rootNode->first_node("AzimuthAngleMax")->value() ));
		if (0 != rootNode->first_node("WindowingAlpha"))
			WindowingAlpha = str2real( rootNode->first_node("WindowingAlpha")->value() );
		if (0 != rootNode->first_node("CrystalloHistoM"))
			Settings::CrystalloHistoM = stoul( rootNode->first_node("CrystalloHistoM")->value() );
		mode = E_RECTANGULAR_WINDOW;
		if (0 != rootNode->first_node("WindowingMethod"))
			mode = stoul( rootNode->first_node("WindowingMethod")->value() );
		switch(mode)
		{
			case E_KAISER_WINDOW:
				WindowingMethod = E_KAISER_WINDOW; break;
			default:
				WindowingMethod = E_RECTANGULAR_WINDOW;
		}
	}

	mode = E_NOSPATSTAT;
	if (0 != rootNode->first_node("AnalysisSpatDistrType")) {
		mode = stoul( rootNode->first_node("AnalysisSpatDistrType")->value() );
		switch (mode)
		{
			case E_RDF:
				SpatialDistributionTask = E_RDF;
				SpatStatDoRDF = true;
				break;
			case E_NEAREST_NEIGHBOR:
				SpatialDistributionTask = E_NEAREST_NEIGHBOR;
				SpatStatDo1NN = true;
				break;
			case E_RIPLEYK:
				SpatialDistributionTask = E_RIPLEYK;
				SpatStatDoRIPK = true;
				break;
			case E_MKNN:
				SpatialDistributionTask = E_MKNN;
				SpatStatDoMKNN = true;
				break;
			case E_NPOINTCORR:
				SpatialDistributionTask = E_NPOINTCORR;
				SpatStatDoNPCorr = true;
				break;
			case E_COUNTNEIGHBORS:
				SpatialDistributionTask = E_COUNTNEIGHBORS;
				SpatStatDoCountNbors = true;
				break;
			default:
				//either zero or multiple mode, so check for E_MULTISPATSTAT
				string modestr = rootNode->first_node("AnalysisSpatDistrType")->value();
				if ( modestr.find("0") != string::npos || modestr.length() > 4 || modestr.length() < 1 ) {
					//if multitask string key contains a zero, more than four different or less than one it is interpreted as E_NOSPATSTAT
					//there are only four spatstat tasks with key 1,2,3,4 but the modestring contains more or no characters
					SpatialDistributionTask = E_NOSPATSTAT;
					break;
				}
				//##MK::is there a more elegant solution to not hardcode the enum variable value but use E_RDF instead?
				if ( modestr.find("1") != string::npos ) { SpatStatDoRDF = true; SpatialDistributionTask = E_MULTISPATSTAT; }
				if ( modestr.find("2") != string::npos ) { SpatStatDo1NN = true; SpatialDistributionTask = E_MULTISPATSTAT; }
				if ( modestr.find("3") != string::npos ) { SpatStatDoRIPK = true; SpatialDistributionTask = E_MULTISPATSTAT; }
				if ( modestr.find("4") != string::npos ) { SpatStatDoMKNN = true; SpatialDistributionTask = E_MULTISPATSTAT; }
				if ( modestr.find("5") != string::npos ) { SpatStatDoNPCorr = true; SpatialDistributionTask = E_MULTISPATSTAT; }
				if ( modestr.find("6") != string::npos ) { SpatStatDoCountNbors = true; SpatialDistributionTask = E_COUNTNEIGHBORS; }
		}
	}

	mode = E_NOCLUST;
	if (0 != rootNode->first_node("AnalysisClusteringType")) {
		mode = stoul( rootNode->first_node("AnalysisClusteringType")->value() );
		switch (mode)
		{
		case E_DBSCAN:
			ClusteringTask = E_DBSCAN; break;
		case E_MAXSEPARATION:
			ClusteringTask = E_MAXSEPARATION; break;
		case E_ISOSURF:
			ClusteringTask = E_ISOSURF; break;
		default:
			ClusteringTask = E_NOCLUST;
		}
	}


	mode = E_NOTESS;
	if (0 != rootNode->first_node("AnalysisVolumeTessellation"))
		mode = stoul( rootNode->first_node("AnalysisVolumeTessellation")->value());
	switch (mode)
	{
		case E_VORONOI_TESS_LEANSTORE:
			VolumeTessellation = E_VORONOI_TESS_LEANSTORE; break;
		default:
			VolumeTessellation = E_NOTESS;
	}


	/*if ( Settings::AnalysisMode == E_ANALYZE_IN_DETECTOR_SPACE ) {
		if (0 != rootNode->first_node("DetBinWidth"))
			DetBinWidth = str2real( rootNode->first_node("DetBinWidth")->value() );
		if (0 != rootNode->first_node("SDMBinWidth"))
			SDMBinWidth = str2real( rootNode->first_node("SDMBinWidth")->value() );
	}*/

	if ( Settings::ReconstructionAlgo == E_RECON_BAS_ETAL ) { /*|| Settings::AnalysisMode == E_ANALYZE_IN_DETECTOR_SPACE ) {*/
		//Reconstruction relevant
		if (0 != rootNode->first_node("FlightLength"))
			FlightLength = str2real( rootNode->first_node("FlightLength")->value() );
		if (0 != rootNode->first_node("AtomicDensity"))
			AtomicDensity = str2real( rootNode->first_node("AtomicDensity")->value() );
		if (0 != rootNode->first_node("EvaporationField"))
			EvaporationField = str2real( rootNode->first_node("EvaporationField")->value() );
		if (0 != rootNode->first_node("DetEffMin"))
			DetEffMin = str2real( rootNode->first_node("DetEffMin")->value() );
		if (0 != rootNode->first_node("DetEffIncr"))
			DetEffIncr = str2real( rootNode->first_node("DetEffIncr")->value() );
		if (0 != rootNode->first_node("DetEffMax"))
			DetEffMax = str2real( rootNode->first_node("DetEffMax")->value() );
		if (0 != rootNode->first_node("KFMin"))
			KFMin = str2real( rootNode->first_node("KFMin")->value() );
		if (0 != rootNode->first_node("KFIncr"))
			KFIncr = str2real( rootNode->first_node("KFIncr")->value() );
		if (0 != rootNode->first_node("KFMax"))
			KFMax = str2real( rootNode->first_node("KFMax")->value() );
		if (0 != rootNode->first_node("ICFMin"))
			ICFMin = str2real( rootNode->first_node("ICFMin")->value() );
		if (0 != rootNode->first_node("ICFIncr"))
			ICFIncr = str2real( rootNode->first_node("ICFIncr")->value() );
		if (0 != rootNode->first_node("ICFMax"))
			ICFMax = str2real( rootNode->first_node("ICFMax")->value() );
	}
	
	//Surface reconstruction relevant
	if (0 != rootNode->first_node("AdvIonPruneBinWidthMin"))
		AdvIonPruneBinWidthMin = str2real( rootNode->first_node("AdvIonPruneBinWidthMin")->value() );
	if (0 != rootNode->first_node("AdvIonPruneBinWidthIncr"))
		AdvIonPruneBinWidthIncr = str2real( rootNode->first_node("AdvIonPruneBinWidthIncr")->value() );
	if (0 != rootNode->first_node("AdvIonPruneBinWidthMax"))
		AdvIonPruneBinWidthMax = str2real( rootNode->first_node("AdvIonPruneBinWidthMax")->value() );
	DebugComputeDistance = false;
	if (0 != rootNode->first_node("DebugComputeDistance")) {
		mode = stoul( rootNode->first_node("DebugComputeDistance")->value() );
		if ( mode == 1 )
			DebugComputeDistance = true;
	}

	if ( SpatialDistributionTask != E_NOSPATSTAT ) {
		//descriptive spatial statistics
		DescrStatTasksCode = "";
		if (0 != rootNode->first_node("DescrStatTasksCode"))
			DescrStatTasksCode = rootNode->first_node("DescrStatTasksCode")->value();
		if (0 != rootNode->first_node("SpatStatRadiusMin"))
			SpatStatRadiusMin = str2real( rootNode->first_node("SpatStatRadiusMin")->value() );
		if (0 != rootNode->first_node("SpatStatRadiusIncr"))
			SpatStatRadiusIncr = str2real( rootNode->first_node("SpatStatRadiusIncr")->value() );
		if (0 != rootNode->first_node("SpatStatRadiusMax"))
			SpatStatRadiusMax = str2real( rootNode->first_node("SpatStatRadiusMax")->value() );
		if (0 != rootNode->first_node("SpatStatKNNOrder")) {
			size_t val = stoul( rootNode->first_node("SpatStatKNNOrder")->value() );
			if ( val > 0 ) {
				SpatStatKNNOrder = val-1;
			}
		}
		else {
			SpatStatKNNOrder = 1-1;
		}
		if ( SpatStatDoMKNN == true ) {
			if (0 != rootNode->first_node("SpatStatMKNNCode"))
				DescrStatMKNNCode = rootNode->first_node("SpatStatMKNNCode")->value();
		}
		SpatStatAddLabelRnd = false;
		if (0 != rootNode->first_node("SpatStatAdditionalLabelRandomization")) {
			mode = stoul( rootNode->first_node("SpatStatAdditionalLabelRandomization")->value() );
			if ( mode == 1 )
				SpatStatAddLabelRnd = true;
		}
	}

	if ( Settings::ClusteringTask != E_NOCLUST ) {
		//clustering analysis
		ClusteringTasksCode = "";
		if (0 != rootNode->first_node("ClusteringTasksCode"))
			ClusteringTasksCode = rootNode->first_node("ClusteringTasksCode")->value();

		if (0 != rootNode->first_node("ClustMaxSepDmaxMin"))
			ClustMSDmaxMin = str2real(rootNode->first_node("ClustMaxSepDmaxMin")->value());
		if (0 != rootNode->first_node("ClustMaxSepDmaxIncr"))
			ClustMSDmaxIncr = str2real(rootNode->first_node("ClustMaxSepDmaxIncr")->value());
		if (0 != rootNode->first_node("ClustMaxSepDmaxMax"))
			ClustMSDmaxMax = str2real(rootNode->first_node("ClustMaxSepDmaxMax")->value());
		if (0 != rootNode->first_node("ClustMaxSepCntMin"))
			ClustMSCntMin = stoul( rootNode->first_node("ClustMaxSepCntMin")->value() );
		if (0 != rootNode->first_node("ClustMaxSepCntIncr"))
			ClustMSCntIncr = stoul( rootNode->first_node("ClustMaxSepCntIncr")->value() );
		if (0 != rootNode->first_node("ClustMaxSepCntMax"))
			ClustMSCntMax = stoul( rootNode->first_node("ClustMaxSepCntMax")->value() );
		if (0 != rootNode->first_node("ClustMaxSepNmin"))
			ClustMSNmin = str2real(rootNode->first_node("ClustMaxSepNmin")->value());
		ClustPostSpatStat = false;
		if (0 != rootNode->first_node("ClustAPosterioriSpatStat")) {
			mode = stoul( rootNode->first_node("ClustAPosterioriSpatStat")->value() );
			if ( mode == 1 )
				ClustPostSpatStat = true;
		}
	}

	if ( Settings::VolumeTessellation != E_NOTESS ) {
		if (0 != rootNode->first_node("SurfaceCellsCarvingRadius"))
			Settings::SurfaceCellsCarvingRadius = str2real( rootNode->first_node("SurfaceCellsCarvingRadius")->value() );
	}

	//visualization options
	IOReconstruction = false;
	if (0 != rootNode->first_node("IOReconstruction")) {
		mode = stoul( rootNode->first_node("IOReconstruction")->value() );
		if ( mode == 1 )
			IOReconstruction = true;
	}
	IOTriangulation = false;
	if (0 != rootNode->first_node("IOTriangulation")) {
		mode = stoul( rootNode->first_node("IOTriangulation")->value() );
		if ( mode == 1 )
			IOTriangulation = true;
	}
	IOTriangulationBVH = false;
	if (0 != rootNode->first_node("IOTriangulationBVH")) {
		mode = stoul( rootNode->first_node("IOTriangulationBVH")->value() );
		if ( mode == 1 )
			IOTriangulationBVH = true;
	}
	IOKDTreePartitioning = false;
	if (0 != rootNode->first_node("IOKDTreePartitioning")) {
		mode = stoul( rootNode->first_node("IOKDTreePartitioning")->value() );
		if ( mode == 1 )
			IOKDTreePartitioning = true;
	}
	IOHKFilteredIons = false;
	if (0 != rootNode->first_node("IOHKFilteredIons")) {
		mode = stoul( rootNode->first_node("IOHKFilteredIons")->value() );
		if ( mode == 1 )
			IOHKFilteredIons = true;
	}
	IORAWHKClusterID = false;
	if (0 != rootNode->first_node("IOHKRawClusterID")) {
		mode = stoul( rootNode->first_node("IOHKRawClusterID")->value() );
		if ( mode == 1 )
			IORAWHKClusterID = true;
	}
	IOIonTipSurfDists = false;
	if (0 != rootNode->first_node("IOIonTipSurfDistances")) {
		mode = stoul( rootNode->first_node("IOIonTipSurfDistances")->value() );
		if ( mode == 1 )
			IOIonTipSurfDists = true;
	}
	IOVoronoiDescrStats = false;
	if (0 != rootNode->first_node("IOVoronoiDescrStats")) {
		mode = stoul( rootNode->first_node("IOVoronoiDescrStats")->value() );
		if ( mode == 1 )
			IOVoronoiDescrStats = true;
	}
	IOVoronoiCellPositions = false;
	if (0 != rootNode->first_node("IOVoronoiCellPositions")) {
		mode = stoul( rootNode->first_node("IOVoronoiCellPositions")->value() );
		if ( mode == 1 )
			IOVoronoiCellPositions = true;
	}
	IOVoronoiTopoGeom = false;
	if (0 != rootNode->first_node("IOVoronoiTopoGeom")) {
		mode = stoul( rootNode->first_node("IOVoronoiTopoGeom")->value() );
		if ( mode == 1 )
			IOVoronoiTopoGeom = true;
	}
	IOCrystallography = false;
	if (0 != rootNode->first_node("IOCrystallography")) {
		mode = stoul( rootNode->first_node("IOCrystallography")->value() );
		if ( mode == 1 )
			IOCrystallography = true;
	}

	//tip synthesis relevant
	if (0 != rootNode->first_node("SimRelBottomRadius"))
		RelBottomRadius = str2real( rootNode->first_node("SimRelBottomRadius")->value() );
	if (0 != rootNode->first_node("SimRelTopRadius"))
		RelTopRadius = str2real( rootNode->first_node("SimRelTopRadius")->value() );
	if (0 != rootNode->first_node("SimRelBottomCapHeight"))
		RelBottomCapHeight = str2real( rootNode->first_node("SimRelBottomCapHeight")->value() );
	if (0 != rootNode->first_node("SimRelTopCapHeight"))
		RelTopCapHeight = str2real( rootNode->first_node("SimRelTopCapHeight")->value() );
	if (0 != rootNode->first_node("SimMatrixLatticeConstant"))
		LatticeConstant = str2real( rootNode->first_node("SimMatrixLatticeConstant")->value() );
	if (0 != rootNode->first_node("SimNumberOfAtoms"))
		NumberOfAtoms = static_cast<size_t>(str2real( rootNode->first_node("SimNumberOfAtoms")->value() ));
	if (0 != rootNode->first_node("SimDetectionEfficiency"))
		SimDetEfficiency = str2real( rootNode->first_node("SimDetectionEfficiency")->value() );
	if (0 != rootNode->first_node("SimFiniteSpatResolutionX"))
		SpatResolutionSigmaX = str2real( rootNode->first_node("SimFiniteSpatResolutionX")->value() );
	if (0 != rootNode->first_node("SimFiniteSpatResolutionY"))
		SpatResolutionSigmaY = str2real( rootNode->first_node("SimFiniteSpatResolutionY")->value() );
	if (0 != rootNode->first_node("SimFiniteSpatResolutionZ"))
		SpatResolutionSigmaZ = str2real( rootNode->first_node("SimFiniteSpatResolutionZ")->value() );
	if (0 != rootNode->first_node("SimNumberOfCluster"))
		NumberOfCluster = static_cast<size_t>(str2real(rootNode->first_node("SimNumberOfCluster")->value()));
	if (0 != rootNode->first_node("SimClusterRadiusMean"))
		ClusterRadiusMean = str2real(rootNode->first_node("SimClusterRadiusMean")->value());
	if (0 != rootNode->first_node("SimClusterRadiusSigmaSqr"))
		ClusterRadiusSigmaSqr = str2real(rootNode->first_node("SimClusterRadiusSigmaSqr")->value());

	//performance
	if (0 != rootNode->first_node("UseNUMABinding")) {
		mode = stoul( rootNode->first_node("UseNUMABinding")->value() );
		if ( mode == 1 )
			NumaBinding = E_DEFAULTBINDING;
	}

	//##MK::convert units to SI
}


bool Settings::checkUserInput()
{
	//##MK::check user input for validity and good sense

	cout << "PARAPROBE utilizes the following settings..." << endl;

	//##MK::add input detector space and reconstruction space AnalysisMode verbosing

	/*if ( Settings::AnalysisMode == E_ANALYZE_IN_DETECTOR_SPACE ) {
		if ( Settings::DetBinWidth < EPSILON ) {
			cout << "Detector binning width needs to be positive and finite!"; return false;
		}
		if ( Settings::SDMBinWidth < EPSILON ) {
			cout << "SignalMapping binning width needs to be positive and finite!"; return false;
		}

		cout << "\t\t\tDetBinWidth\t\t\t" << Settings::DetBinWidth << "\n";
		cout << "\t\t\tSDMBinWidht\t\t\t" << Settings::SDMBinWidth << "\n";
	}*/

	//is input file format consistent with reconstruction mode?
	cout << "Input/Reconstruction" << endl;
	//if ( ReconstructionAlgo == E_RECON_READFROMSYNTH ) {
	//	if ( InputFileformat != E_IO_HDF5 ) {
	//		cout << "ReconstructionAlgorithm is HDF5 but InputFileformat is not HDF5!"; return false;
	//	}
	//}
	if ( ReconstructionAlgo == E_RECON_READFROMPOS ) {
		if ( InputFileformat != E_IO_POS && InputFileformat != E_IO_HDF5 ) {
			cout << "ReconstructionAlgorithm is POS but InputFileformat is neither POS nor HDF5!"; return false;
		}
	}
	if ( ReconstructionAlgo == E_RECON_READFROMEPOS ) {
		if ( InputFileformat != E_IO_EPOS ) {
			cout << "ReconstructionAlgorithm is EPOS but InputFileformat is not EPOS!"; return false;
		}
	}
	if ( ReconstructionAlgo == E_RECON_BAS_ETAL ) {
		if ( InputFileformat != E_IO_EPOS ) {
			cout << "Reconstruction is Bas et al. but InputFileformat is not EPOS!"; return false;
		}
	}
	if ( SyntheticTips == true ) {
		cout << "SYNTHETIC tip" << endl;
		if ( RelBottomRadius <= 0.f ) {
			cout << "RelBottomRadius must be positive and non-zero!" << endl; return false;
		}
		if ( RelTopRadius <= 0.f || RelTopRadius > RelBottomRadius ) {
			cout << "RelTopRadius must not be larger than bottom, positive and non-zero!" << endl; return false;
		}
		if ( RelBottomCapHeight <= 0.f || RelBottomCapHeight > RelBottomRadius ) { //##MK::carve out not more than a halfsphere from bottom of the frustum
			cout << "RelBottomCapHeight must be positive and non-zero and not larger than RelBottomRadius!" << endl; return false;
		}
		if ( RelTopCapHeight <= 0.f || RelTopCapHeight > RelTopRadius ) {
			cout << "RelTopCapHeight must be positive and non-zero and not larger than RelTopRadius!" << endl; return false; //##MK::add further constraints on cap heights
		}
		if ( LatticeConstant <= 0.f ) {
			cout << "LatticeConstant must be positive and non-zero!" << endl; return false;
		}
		if ( NumberOfAtoms <= 0 ) {
			cout << "The tip must contain at least one atom!" << endl; return false;
		}
		if ( SimDetEfficiency <= 0.f || SimDetEfficiency > 1.f ) {
			cout << "The simulated detection efficiency needs to on interval (0,1]!" << endl; return false;
		}
		if ( SpatResolutionSigmaX < 0.f || SpatResolutionSigmaY < 0.f || SpatResolutionSigmaZ < 0.f ) {
			cout << "The simulated finite spatial resolutions must be positive if not zero!" << endl; return false;
		}
		if ( NumberOfCluster < 0 ) {
			cout << "The number of cluster cannot be negative if cluster are desired!" << endl; return false;
		}
		//##MK::at the moment using distribution mu and sigma^2 instead of mean and variance interpreting values in nanometer
		if ( ClusterRadiusMean < 0.f ) {
			cout << "The mean of the lognormal cluster size distribution must not be negative or non-zero!" << endl; return false;
		}
		if ( ClusterRadiusSigmaSqr <= EPSILON ) {
			cout << "The sigma(variance) of the lognormal cluster size distribution must not be negative or non-zero!" << endl; return false;
		}
		//input checks passed
		cout << "\t\tBuild a synthetic tip as follows\n";
		cout << "\t\t\tRelBottomRadius\t\t\t" << RelBottomRadius << "\n";
		cout << "\t\t\tRelTopRadius\t\t\t" << RelTopRadius << "\n";
		cout << "\t\t\tRelBottomCapHeight\t\t" << RelBottomCapHeight << "\n";
		cout << "\t\t\tRelTopCapHeight\t\t\t" << RelTopCapHeight << "\n";
		cout << "\t\t\tLatticeConstant fcc\t\t" << LatticeConstant << "\n";
		cout << "\t\t\tNumber of atoms\t\t\t" << NumberOfAtoms << "\n";
		cout << "\t\t\tDetectionEfficiency\t\t" << SimDetEfficiency << "\n";
		cout << "\t\t\tFiniteSpatResSigmaX\t\t" << SpatResolutionSigmaX << "\n";
		cout << "\t\t\tFiniteSpatResSigmaY\t\t" << SpatResolutionSigmaY << "\n";
		cout << "\t\t\tFiniteSpatResSigmaZ\t\t" << SpatResolutionSigmaZ << "\n";
		cout << "\t\t\tNumber of cluster\t\t" << NumberOfCluster << "\n";
		cout << "\t\t\tClusterRadiusMean\t\t" << ClusterRadiusMean << "\n";
		cout << "\t\t\tClusterRadiusVar\t\t" << ClusterRadiusSigmaSqr << "\n";
	}
	else {
		if ( ReconstructionAlgo == E_RECON_READFROMSYNTH) {
			cout << "\t\tAccept synthetic dataset\n";
		}
		if ( ReconstructionAlgo == E_RECON_READFROMEPOS ) {
			cout << "\t\tAccept EPOS reconstruction\n";
		}
		else if ( ReconstructionAlgo == E_RECON_READFROMPOS ) {
			cout << "\t\tAccept POS reconstruction\n";
		}
		else if ( ReconstructionAlgo == E_RECON_BAS_ETAL ) {
			cout << "\t\tBarr Reconstruction Algorithm (default)\n";
			//##MK::add SI units
			cout << "\t\t\tFlightLength\t\t\t" << FlightLength << "\n";
			cout << "\t\t\tAtomicDensity\t\t\t" << AtomicDensity << "\n";
			cout << "\t\t\tEvaporaField\t\t\t" << EvaporationField << "\n";
			cout << "\t\t\tDefEffMin\t\t\t" << DetEffMin << "\n";
			cout << "\t\t\tDefEffIncr\t\t\t" << DetEffIncr  << "\n";
			cout << "\t\t\tDefEffMax\t\t\t" << DetEffMax << "\n";
			cout << "\t\t\tKFMin\t\t\t\t" << KFMin << "\n";
			cout << "\t\t\tKFIncr\t\t\t\t" << KFIncr << "\n";
			cout << "\t\t\tKFMax\t\t\t\t" << KFMax << "\n";
			cout << "\t\t\tICFMin\t\t\t\t" << ICFMin << "\n";
			cout << "\t\t\tICFIncr\t\t\t\t" << ICFIncr << "\n";
			cout << "\t\t\tICFMax\t\t\t\t" << ICFMax << "\n";
		}
		else {
			cout << "\t\tUnknown reconstruction algorithm!\n";
			return false;
		}
	}

	cout << "\t\t" << RAWFilenameIn << " I/O for heavy data\n";
	cout << "\t\t" << RRNGFilenameIn << " I/O for ranging\n";
	if ( IdentifyIonType == true )
		cout << "\t\tIdentifying ion type\n";
	else
		cout << "\t\tAssuming all ions of the same UNKNOWNTYPE!\n";

	if ( Settings::SurfaceMeshingAlgo == E_LOADEXISTENTVTK ) {
		//compare if filename Raw and VTKTriangles is consistent
		string ivasdata = Settings::RAWFilenameIn.substr(0,Settings::RAWFilenameIn.length()-4-1); //strip sope .
		string tridata = Settings::SurfaceFilenameIn.substr(0,Settings::SurfaceFilenameIn.length()-3-1-4-1); //strip ktv.sope.
		cout << ivasdata << "\n";
		cout << tridata << "\n";
		if ( ivasdata.compare(tridata) == 0 ) { //strip vtk
			cout << "Input dataset and prescribed triangle hull are consistent" << endl;
		}
		else {
			cout << "Input dataset and prescribed triangle hull are not consistent!" << endl;
			//##MK::return false;
		}
	}
	if ( Settings::SurfaceMeshingAlgo == E_ALPHASHAPE_CGAL ) {
		cout << "Building an alpha shape using the CGAL library" << endl;
		if ( Settings::SurfaceAShapeAValue == E_ASHAPE_SMALLEST_SOLID )
			cout << "Using the smallest alpha value to form a solid for this" << endl;
		else
			cout << "Using the CGAL's optimality criterion algorithm for this" << endl;
	}
	if ( Settings::AdvIonPruneBinWidthMin > Settings::AdvIonPruneBinWidthMax
			|| Settings::AdvIonPruneBinWidthMin > Settings::AdvIonPruneBinWidthIncr
			|| Settings::AdvIonPruneBinWidthMax < Settings::AdvIonPruneBinWidthIncr ) {
		cout << "Invalid choice for AdvIonPruneBinning" << endl;
		return false;
	}
	cout << "PARAPROBE performs advanced ion pruning for surface detection..." << endl;
	cout << "\t\tAdvIonPruneBinWidthMin\t\t" << AdvIonPruneBinWidthMin << "\n";
	cout << "\t\tAdvIonPruneBinWidthIncr\t\t" << AdvIonPruneBinWidthIncr << "\n";
	cout << "\t\tAdvIonPruneBinWidthMax\t\t" << AdvIonPruneBinWidthMax << "\n";
	if ( DebugComputeDistance == true )
		cout << "Distance to tip surface is computed to reduce bias of spatial statistics" << endl;
	else
		cout << "Distance to tip surface is NOT computed, spatial statistics will be biased!" << endl;

	//mining crystallographic information
	if ( ExtractCrystallographicInfo != E_NOCRYSTALLO ) {
		if ( CrystalloRadiusMax < EPSILON ) {
			cout << "\t\tSpatStatRadiusMax must be positive!" << endl; return false;
		}
		if (SamplingGridBinWidthX < EPSILON || SamplingGridBinWidthY < EPSILON || SamplingGridBinWidthZ < EPSILON) {
			cout << "\t\tSamplingGridBinWidths must be positive!" << endl; return false;
		}
		//if ( ElevationAngleMin < 0.f ) {
		//	cout << "\t\tElevationAngleMin must be positive!" << endl; return false;
		//}
		if ( ElevationAngleIncr < EPSILON ) {
			cout << "\t\tElevationAngleIncr must be positive!" << endl; return false;
		}
		//if ( ElevationAngleMin >= ElevationAngleMax ) {
		//	cout << "\t\tElevationAngleMin is larger than ElevationAngleMax!" << endl; return false;
		//}
		//if ( AzimuthAngleMin < 0.f ) {
		//	cout << "\t\tAzimuthAngleMin must be positive!" << endl; return false;
		//}
		if ( AzimuthAngleIncr < EPSILON ) {
			cout << "\t\tAzimuthAngleIncr must be positive!" << endl; return false;
		}
		//if ( AzimuthAngleMin >= AzimuthAngleMax ) {
		//	cout << "\t\tAzimuthAngleMin is larger than AzimuthAngleMax!" << endl; return false;
		//}
		if ( Settings::WindowingMethod == E_KAISER_WINDOW ) {
			if ( WindowingAlpha < EPSILON ) {
				cout << "\t\tKaiser windows with alpha = 0 are rectangular windows, use rectangular window instead!" << endl; return false;
			}
		}
		if ( Settings::CrystalloHistoM < 9 ) {
			cout << "\t\tWarning::Histogram binning was very small, reset to m=9" << endl;
			CrystalloHistoM = 9;
		}
		if ( Settings::CrystalloHistoM > 10) {
			cout << "\t\tWarning::Histogram binning was very large, reset to m=10" << endl;
			CrystalloHistoM = 10;
		}
		cout << "PARAPROBE performing Araullo-Peters et al. method for crystallographic feature extraction" << endl;
		cout << "\t\tCrystalloRadiusMax\t\t" << CrystalloRadiusMax << "\n";
		cout << "\t\tSamplingGridBinWidthX\t\t" << SamplingGridBinWidthX << "\n";
		cout << "\t\tSamplingGridBinWidthY\t\t" << SamplingGridBinWidthY << "\n";
		cout << "\t\tSamplingGridBinWidthZ\t\t" << SamplingGridBinWidthZ << "\n";
		cout << "\t\tElevationAngleMin\t\t" << ElevationAngleMin << "\n";
		cout << "\t\tElevationAngleIncr\t\t" << ElevationAngleIncr << "\n";
		cout << "\t\tElevationAngleMax\t\t" << ElevationAngleMax << "\n";
		cout << "\t\tAzimuthAngleMin\t\t" << AzimuthAngleMin << "\n";
		cout << "\t\tAzimuthAngleIncr\t\t" << AzimuthAngleIncr << "\n";
		cout << "\t\tAzimuthAngleMax\t\t" << AzimuthAngleMax << "\n";
		cout << "\t\tHistogram bin\t\t" << CrystalloHistoM << endl;
		if ( WindowingMethod == E_RECTANGULAR_WINDOW )
			cout << "A rectangular window is used prior discrete Fourier transform" << "\n";
		else
			cout << "A Kaiser window with alpha = " << WindowingAlpha << " is used prior discrete Fourier transform" << "\n";
	}

	//descriptive spatial statistics
	if ( SpatialDistributionTask != E_NOSPATSTAT ) {
		if ( SpatStatRadiusMin < 0.f ) {
			cout << "\t\tSpatStatRadiusMin must be positive or zero!" << endl; return false;
		}
		if ( SpatStatRadiusIncr < EPSILON ) {
			cout << "\t\tSpatStatRadiusIncr must be positive!" << endl; return false;
		}
		if ( SpatStatRadiusMax < EPSILON ) {
			cout << "\t\tSpatStatRadiusMax must be positive!" << endl; return false;
		}
		if ( SpatStatRadiusMin >= SpatStatRadiusMax ) {
			cout << "\t\tSpatStatRadiusMin is larger than SpatStatRadiusMax!" << endl; return false;
		}
		//##MK::through interval test maybe helper function
		/*		
		cout << setprecision(32) << floor(static_cast<double>(SpatStatRadiusMin) / static_cast<double>(SpatStatRadiusIncr)) << endl;
		cout << setprecision(32) << ceil(static_cast<double>(SpatStatRadiusMin) / static_cast<double>(SpatStatRadiusIncr)) << endl;
		if ( static_cast<size_t>( floor(static_cast<double>(SpatStatRadiusMin) / static_cast<double>(SpatStatRadiusIncr)) ) != 
			static_cast<size_t>( ceil(static_cast<double>(SpatStatRadiusMin) /static_cast<double>(SpatStatRadiusIncr)) ) ) {
			cout << "\t\tSpatStatRadiusMin is not an integer multiple of the width!" << endl; return false;
		}
		cout << setprecision(32) << floor(static_cast<double>(SpatStatRadiusMax) / static_cast<double>(SpatStatRadiusIncr)) << endl;
		cout << setprecision(32) << ceil(static_cast<double>(SpatStatRadiusMax) / static_cast<double>(SpatStatRadiusIncr)) << endl;
		if ( static_cast<size_t>( floor(static_cast<double>(SpatStatRadiusMax) / static_cast<double>(SpatStatRadiusIncr)) ) != 
			static_cast<size_t>( ceil(static_cast<double>(SpatStatRadiusMax) /static_cast<double>(SpatStatRadiusIncr)) ) ) {
			cout << "\t\tSpatStatRadiusMax is not an integer multiple of the width!" << endl; return false;
		}
		*/
		cout << "PARAPROBE spatial statistics..." << endl;
		//cout << "\t\tSpatStat taskcode\t\t" << DescrStatTaskCode << "\n";
		cout << "\t\tSpatStatRadiusMin\t\t" << SpatStatRadiusMin << "\n";
		cout << "\t\tSpatStatRadiusIncr\t\t" << SpatStatRadiusIncr << "\n";
		cout << "\t\tSpatStatRadiusMax\t\t" << SpatStatRadiusMax << "\n";
		cout << "\t\tSpatStatKNNOrder\t\t" << SpatStatKNNOrder << "\n";
		if ( SpatStatDoRDF == true )	cout << "\t\tDoing RDF\n";
		if ( SpatStatDo1NN == true )	cout << "\t\tDoing 1NN\n";
		if ( SpatStatDoRIPK == true )	cout << "\t\tDoing RipleyK\n";
		if ( SpatStatDoMKNN == true )	cout << "\t\tDoing MKNN\n";
		if ( SpatStatDoNPCorr == true ) cout << "\t\tDoing nPoint correlations\n";
		if ( SpatStatDoCountNbors == true ) cout << "\t\tDoing count neighbors\n";
		if ( SpatStatAddLabelRnd == true )
			cout << "\t\tAdditional run with randomized labels\n";
	}
	//clustering
	if ( Settings::ClusteringTask != E_NOCLUST ) {
		if ( ClustMSDmaxMin < EPSILON ) {
			cout << "\t\tClustMaxSepDmaxMin must be positive and smallest of range values!" << endl; return false;
		}
		if ( ClustMSDmaxIncr < EPSILON || ClustMSDmaxIncr > ClustMSDmaxMax) {
			cout << "\t\tClustMaxSepDmaxIncr must be positive and within min max range!" << endl; return false;
		}
		if ( ClustMSDmaxMax < EPSILON || ClustMSDmaxMax < ClustMSDmaxMin ) {
			cout << "\t\tClustMaxSepDmaxMax must be positive and largest of range!" << endl; return false;
		}
		//if ( ClustMSCntMin < 2 || ClustMSCntMax < ClustMSCntMin || ClustMSCntIncr < 2 ) {
		//	cout << "\t\tClustMaxSepN parameter need to be at least 2" << endl; return false;
		//}
		if ( ClustMSNmin < 2 ) {
			cout << "\t\tMaxSeparationMethod Nmin must be at least 2 !" << endl; return false;
		}
		if ( ClustPostSpatStat == true )
			cout << "\t\tDoing clustering a posteriori spatstat\n";

		cout << "PARAPROBE clustering analyses..." << endl;
		//cout << "\t\tMaxSeparationMethod Taskcode\t\t" << ClusteringTaskCode << "\n";
		cout << "\t\tMaxSepMethod DmaxMin\t\t" << ClustMSDmaxMin << "\n";
		cout << "\t\tMaxSepMethod DmaxIncr\t\t" << ClustMSDmaxIncr << "\n";
		cout << "\t\tMaxSepMethod DmaxMax\t\t" << ClustMSDmaxMax << "\n";
		//cout << "\t\tMaxSepMethod N   Min\t\t" << ClustMSCntMin << "\n";
		//cout << "\t\tMaxSepMethod N   Incr\t\t" << ClustMSCntIncr << "\n";
		//cout << "\t\tMaxSepMethod N   Max\t\t" << ClustMSCntMax << "\n";
		cout << "\t\tMaxSepMethod Nmin\t\t\t" << ClustMSNmin << "\n";

	}
	if ( Settings::VolumeTessellation != E_NOTESS ) {
		cout << "Chris Rycrofts Voro++ library is wrapped and used threaded to compute a Voronoi tessellation of the tip" << "\n";
		cout << "Tessellation guard zone using width " << Settings::TessellationGuardWidth << " times mean distance" << "\n";
		cout << "Tessellation number of points per block " << Settings::TessellationPointsPerBlock << "\n";
		cout << "Tessellation discarding cells closer than " << Settings::SurfaceCellsCarvingRadius << " nm to surface" << "\n";
		if ( Settings::VolumeTessellation == E_VORONOI_TESS_LEANSTORE)
			cout << "The Voronoi lean storage model is used" << "\n";
	}

	cout << "PARAPROBE I/O options..." << endl;
	if ( IOReconstruction == true )
		cout << "\t\tReconstruction" << "\n";
	if ( IOTriangulation == true )
		cout << "\t\tTriangulation" << "\n";
	if ( IOHKFilteredIons == true )
		cout << "\t\tHoshenKopelmanFiltered" << "\n";
	if ( IORAWHKClusterID == true )
		cout << "\t\tHoshenKopelmanRawClusterID" << "\n";
	if ( IOIonTipSurfDists == true )
		cout << "\t\tVTKIonTipDistance" << "\n";
	if ( IOVoronoiDescrStats == true )
		cout << "\t\tHDF5VoronoiDescrStats" << "\n";
	if ( IOVoronoiCellPositions == true )
		cout << "\t\tHDF5VoronoiCellPositions" << "\n";
	if ( IOVoronoiTopoGeom == true )
		cout << "\t\tHDF5VoronoiTopoGeom" << "\n";
	if ( IOCrystallography == true )
		cout << "\t\tHDF5Crystallography" << "\n";

	cout << "PARAPROBE Performance options..." << endl;
	if ( Settings::NumaBinding == E_NOBINDING )
		cout << "\t\tThreads get not pinned" << "\n";
	if ( Settings::NumaBinding == E_DEFAULTBINDING )
		cout << "\t\tThreads are pinned using NUMA library" << "\n";

	//cout << "All console output which follows is intended for debugging primarily..." << endl;
	cout << endl;
	return true;
}

/*apt_real Settings::read_real( xml_node<>* const src, const string keyword, const apt_real defaultvalue )
{
	if (0 != src->first_node(keyword))
		return str2real( src->first_node(keyword)->value() );
	else
		return defaultvalue;
}*/

