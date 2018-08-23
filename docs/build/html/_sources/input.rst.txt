**XML Control File Settings**
=============================
The entire datamining is controlled by a single control parameter settings XML file.

Input
^^^^^
| **InputFileformat**
|	Determines where the rawdata come from
|	0, debug, nothing
|	1, POS
|	2, EPOS
|	3, HDF5 use only to reread existent synthetic dataset and ReconstructionAlgorithm 2
|	4, synthetic dataset
| **RAWFilenameIn**
|	Specifies the input file filename an ending is necessary. It can be .pos or .epos

Analysis mode
^^^^^^^^^^^^^
| **AnalysisMode**
|	0, default

Reconstruction
^^^^^^^^^^^^^^
| **ReconstructionAlgorithm**
|	Determines the location of the ions in the reconstruction space
|	0, no reconstruction
|	1, accept x,y,z from synthetic dataset demands InputFileformat 1
|	2, accept x,y,z coordinates from the POS file
|	3, accept x,y,z coordinates from the EPOS file
|	4, perform common Bas et al. reconstruction based on EPOS file
| **IdentifyIonType**
|	If other than 0, ranging is performed, if invalid range file or set to 0 all ions are assigned a default type
| **RRNGFilenameIn**
|	Specifies a rrng-format-conformant range file with which to convert mass-to-charge to ion type aka ranging

Tip surface reconstruction
^^^^^^^^^^^^^^^^^^^^^^^^^^
| **SurfaceReconstructionType**
|	Determines whether or not and if so which tip surface reconstruction model to apply
|	0, no surface reconstruction
|	1, alpha shape using the CGAL library and its optimal alpha value to get a solid to report an alpha shape
|	4, read existent triangle hull from SurfaceFilenameIn file
| **SurfaceFilenameIn**
|	Specifies a VTK file containing tip surface triangle hull information when using option **SurfaceReconstructionType** 4


Analysis Tasks
^^^^^^^^^^^^^^
| **AnalysisSpatDistrType**
|	Specifies which descriptive spatial statistics should be computed. 
|	Multiple single character numbers can be provided to instruct multiple analyses, eg 42 will do knearest and nearest neighbor
|	0, nothing
|	1, radial distribution function
|	2, nearest neighbor
|	3, Ripley K
|	4, k nearest neighbors
|	5, multiple k nearest neighbors

| **AnalysisClusteringType**
|	Specifies which clustering algoritm to use
|	0, nothing
|	2, maximum separation method

Bas et al reconstruction parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
|	Optional when using reconstruction mode Bas et al
| **FlightLength**
|	In nanometer, instrument dependent
| **AtomicDensity**
|	In atoms per cubic nanometer
| **EvaporationField**
|	In volt per nanometer
| **DetEffMin**
| **DetEffIncr**
| **DetEffMax**
|	Specifies range of detector efficiency minimum, increment, maximum, respectively, physically restricted on (0,1)
| **KFMin**
| **KFIncr**
| **KFMax**
|	Specifies range of kf factor minimum, increment, maximum, respectively
| **ICFMin**
| **ICFIncr**
| **ICFMax**
|	Specifies range of ICF factor

Smart pruning
^^^^^^^^^^^^^
| A technique to identify ions close to the tip surface and dominant poles to avoid downsampling but speed up surface reconstruction.
| **AdvIonPruningBinWidthMin**
| **AdvIonPruningBinWidthIncr**
| **AdvIonPruningBinWidthMax**
|	Specifies range of rectangular binning bin width in nanometer with which to smartly prune the ion cloud
|	before passing candidate points to CGAL to speed up and perform memory leaner surface reconstruction
| **DebugComputeDistance**
|	0 no distance computation to tip surface triangle hull may introduce bias
|	1 exact distance computation to tip surface triangle hull all ions up to SpatStatRadiusMax from the surface

Descriptive spatial statistics parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **DescrStatTaskCode**
|	A semicolon-separated list of string-based codes which analyses on ion types to conduct. Each analysis task, string, requires at least one central ion type string and at least one environment neighboring ion type string both as specified in the rrng rangefile. Multiple central ions in the central ion string can be combined in the analysis against multiple ions in the environment string. For example Al,Al;Ga,Ga will instruct two tasks. The first Al ions against only Al ions as neighbors, the second Gallium against only Gallium. A more complex multiple ion type example is as such AlMn,AlMn;GaMn,GaMn It will instruct in the first analysis to probe Al or Mn as central ions and takes either Al or Mn in the environment, similarly the second Ga or Mn as central against neighboring ions of type either Ga or Mn.
| **SpatStatRadiusMin**
| **SpatStatRadiusIncr**
| **SpatStatRadiusMax**
|	Specifies in nanometer the range of search sphere radii in which to conduct the analyses. Increment needs to be an integer multiple of SpatStatRadiusMax
| **SpatStatKNNOrder**
|	Specifies k for AnalysisSpatDistrType mode 4, output will use C style reporting i.e. order 1 is reported as 0
| **SpatStatMKNNCode**
|	Specifies a semicolon separated list of only non-negative integer values which kth order nearest neighbors should be computed to the central ion, allows to define combinations like 1;2;5;10;100;1000
| **SpatStatAdditionalLabelRandomization**
|	If set to 1 allows to randomize all ion type labels across point cloud and rerun clustering analysis against, otherwise no randomization done. Applied randomizations are reset after the analysis to not invalidate the data set.

Clustering analyses parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **ClusteringTaskCode**
|	Syntax is the same as **DescrStatTaskCode**, different tasks can be defined. 
| **ClustMaxSepDmaxMin**
| **ClustMaxSepDmaxIncr**
| **ClustMaxSepDmaxMax**
|	Specifies in nanometer the range and stepping of the Dmax parameter used to perform a set of independent maximum separation clustering analyses with different Dmax values but same Nmin.
| **ClustMaxSepNmin**
|	Minimum number of ions **inclusive/exclusive** to consider a cluster.
| **ClustAPosterioriSpatStat**
|	If 1 perform a spatial distribution analysis on the clustered ions after the clustering analysis, this support is in beta stage


Visualization options
^^^^^^^^^^^^^^^^^^^^^
If set to value 1 switched on, if set to 0 switched off

| **IOReconstruction**
|	Write ion positions and ranging information to VTK file
| **IOTriangulation**
|	Write tip surface triangle hull into VTK file
| **IOHKFilteredIons**
|	Write candidate ion positions and ranging information to VTK file
| **IOHKClusterID**
|	Write smart pruning bin information to binary file, unsigned int x+yNX+zNXY implicitly encoded
| **IOIonTipSurfDistances**
|	Write distance of ions to surface to VTK file


Synthetic tip
^^^^^^^^^^^^^
Parameter specifying geometry and size of synthetic tip. Shape is conical frustum with spherical cap on top and spherical cap cut out at the bottom.
So far the support for tip geometry is simplistic defining hardcoded single-crystalline pure Al tip with Al3Sc precipitates.

| **SimRelBottomRadius**
| **SimRelTopRadius**
| **SimRelBottomCapHeight**
|	Only frustum height
| **SimRelTopCapHeight**
|	All four relative to frustum height ie. restricted on 0,1
| **SimMatrixLatticeConstant**
|	Currently aluminium single-crystalline pillar.
| **SimNumberOfAtoms**
|	How many ions assuming full efficiency.
| **SimDetectionEfficiency**
|	Fraction of NumberOfAtoms to place, sampling randomly MersenneTwister
| **SimFiniteSpatResolutionX**
| **SimFiniteSpatResolutionY**
| **SimFiniteSpatResolutionZ**
|	Sigma parameter of normal distribution about lattice position by means of which ion is displaced about ideal position, in nanometer
| **SimNumberOfCluster**
|	How many Al3Sc cluster to place in bounding box about the tip
| **SimClusterRadiusMean**
| **SimClusterRadiusSigmaSqr**
|	Lognormal distribution parameter NOT expectation value and variance for cluster size distribution, currently all cluster same size mean, in nanometer


Performance
^^^^^^^^^^^
Options and setting affecting the way PARAPROBE is executed on the machine.

| **UseNUMABinding**
|	If 1, will use the NUMA_ library to pin the threads on specific cores.
| 	Given that the pinning heuristic is not portable, potential code modification are required.

 .. _NUMA: https://linux.die.net/man/3/numa
