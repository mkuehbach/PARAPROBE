.. figure:: ../images/PARAPROBEFront_02.png
..   :scale: 20%
..   :align: left
..   :target: https://github.com/mkuehbach/PARAPROBE
      
| PARAPROBE is an MPI/OpenMP/SIMD-parallelized tool for efficient scalable processing of 
| Atom Probe Tomography (APT) data targeting back-end integration.

| The source code is developed by Markus Kühbach, currently Postdoc with Franz Roters 
| at the Max-Planck-Institut für Eisenforschung GmbH (MPIE_) in Düsseldorf.
| Feel free to utilize the tool, do not hesitate contacting me_ to suggest 
| improvements or desirable features, or simply to report your experiences.

 .. _MPIE: https://www.mpie.de

1. Getting started
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
     
   basics
   setup
   

2. Utilize productively
^^^^^^^^^^^^^^^^^^^^^^^
   
.. toctree::
   :maxdepth: 2
      
   input
   executing

3. Examples
^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   tutorials

Version history
^^^^^^^^^^^^^^^

| **v0.1** --- Initial implementation
|	MPI/OpenMP thread-parallelized spatial range querying and indexing tasks
|	CGAL alpha shape tip surface extraction with smart pruning
|	OMP-parallelized RDF, 1NN, kNN, Ripley K descriptive spatial statistics
|	These can be batch processed using different and multiple type combinations
|	As well as multiple kNN settings at a time for each task, descriptive statistics are
|	computed both randomized and not randomized with full ion typ pairing flexibility
|	OMP-parallelized maximum separation clustering method with dmax parameter sweeping
| Beta functionality
|	A posteriori possibility for each dmax single clustering result to relabel the ions
|	for excluding clustered ions and guardzone about them from the again 
|	OMP-parallel a posteriori clustering descriptive spatial statistics
|	Synthetic tip generation for benchmarking and testing, exemplary HDF5 I/O.


Licence
^^^^^^^

.. toctree::
   :maxdepth: 2
   
   licence

References
^^^^^^^^^^

.. toctree::
   :maxdepth: 2
     
   refs

Funding
^^^^^^^
| The author gratefully acknowledges the support from the DFG in the frame of the project BA 4253/2-1 
| as well as the provisioning of computing resources by the Max-Planck Gesellschaft.

Questions, contributions
^^^^^^^^^^^^^^^^^^^^^^^^
Just let me know and contact me_

 .. _me: https://www.mpie.de/person/51206/2656491
