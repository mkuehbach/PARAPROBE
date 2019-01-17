v0.1
^^^^
* **Initial implementation**
* POS, EPOS reading, RRNG range file parsing
* Barr et al. reconstruction, supports for up to 4.2 billion ion tips
* Generation of synthetic single-crystalline tip structures
* MPI/OpenMP thread parallelized spatial range querying and indexing tasks
* Tip surface extraction through alpha shapes to entire datasets
* Surface extraction made efficient through smart ion pruning pre-processing algorithm
* Floating point precision exact distancing of ions to the alpha shape triangle hull
* This allows to reduce bias in descriptive statistics and tessellation by excluding close to the surface ions from the analyses
* Thread parallel radial distribution function (RDF), k nearest neighbor (kNN), Ripley K
* Thread parallel 2-point descriptive spatial statistics
* In-built batch processing capability for fully automatic processing of such statistics
* Allows for arbitrary single and molecular ion type combinations
* Optional ion type label randomization
* Thread parallel maximum separation clustering method with parameter space sweeping capability
* This can also be combined with the batch processing functionality
* Thread parallel implementation of V. Araullo-Peters et al. reconstruction-space-based method for quantifying crystallographic signal through discrete Fourier analysis
* Thread parallel wrapper around C. Rycrofts Voro++ library to enable hitherho impossible computations of volume tessellations to the entire tip
* Characterize the cell volume to obtain atomic scale concentration values and topology through nearest neighbor analysis and p-vectors
* Hierarchical Data Format (HDF5) / eXtensible Data Model and Format (XDMF) powered results reporting

Beta-stage functionality
^^^^^^^^^^^^^^^^^^^^^^^^
* Optional a posteriori relabeling of ions after each maximum separation clustering run to perform descriptive spatial statistics in population of remaining non-clustered ions using guard zones

