**In a nutshell**
=================

What is PARAPROBE?
^^^^^^^^^^^^^^^^^^
| A software for data mining Atom Probe Tomography (APT) experiment data. It sets prime focus on 
| applying scalable hierarchical parallelism to spatial range querying, clustering, atom probe 
| crystallography, and computational geometry tasks making use of scalable hierarchical parallelism.

What are the user benefits?
^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Open source software**
|	Therefore no usage restriction, unlimited licences when running in parallel
|	surplus full functional transparency and modifiability.
| **Reduced analysis bias**
|	Enabled by state of the art tip surface reconstruction surplus ion to surface distancing.
| **Scalable performance, large datasets**
| 	Thanks to parallel implementation with rigorous hierarchical spatial data partitioning
| 	strategy to improve the utilization of fast caches and ccNUMA-aware data placement policy.

Which parallelization concepts are employed?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Process data parallelism** via the Message Passing Interface (MPI_) API
PARAPROBE processes each individual measurement through a single process. This enables to either distribute parameter sweeping studies of the same tip on practically hundred thousands of processes or to process trivially in parallel multiple tips using the same automatized analysis protocol. At runtime, MPI invokes library calls to communicate pieces of information between physically disjoint computers if necessary. 
As MPI is a library, it requires installation and linking.

**Shared memory thread data parallelism** via the Open Multi-Processing (OpenMP_) API.
PARAPROBE partitions the point data of each measurement into spatially disjoint chunks. Explicit strategies are applied to map and place the data chunks in thread-local memory to reduce false sharing and performance degradation on resources with multiple ccNUMA layers. OpenMP builds on preprocessor directives through which the corresponding pragmas are translated during compilation. As such, OpenMP needs no installation.

.. **Streaming instruction data parallelism (SIMD)** via portable vector intrinsics template libraries (e.g. bSIMD_ or Vc_) is used.
.. At the level of each thread  some core geometrical and numerical tasks can be accelerated further through vectorization. The key idea is to apply vectorized operation which applies to a packet of multiple data elements of the same kind rather than to process single data elements one after another. Technically, this is implementable through highly operation-, problem-, and-CPU-specific instructions of the CPU, the so-called intrinsics.
.. Unfortunately, this renders the code non-portable. Better portability is achieved through portable vector intrinsics. These wrap the individual intrinsics into more abstract commands and compile time instructions with which the choice for the specific realization is delegated to the compiler.

 .. _MPI: http://www.mcs.anl.gov/research/projects/mpi/
 .. _OpenMP: https://www.openmp.org/
 .. _BSIMD: https://developer.numscale.com/bsimd
 .. _Vc: https://github.com/VcDevel/Vc

 
Solid HPC background literature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| J. L. Hennessy and D. A. Patterson
| Computer Architectures: A Quantitative Approach
| 5th edition, 2012, Morgan Kaufmann, Amsterdam

| T. Rauber and G. Rünger
| Parallel Programming for Multicore and Cluster Systems
| 2nd edition, 2013, Springer Heidelberg
| http://dx.doi.org/10.1007/978-3-642-37801-0

| J. Reinders and J. Jeffers
| High Performance Parallelism Pearls Volume One: 
| Multicore and Many-Core Programming Approaches
| 1st edition, 2014, Morgan Kaufmann

| J. Jeffers and J. Reinders
| High Performance Parallelism Pearls Volume Two: 
| Multicore and Many-Core Programming Approaches
| 1st edition, 2015, Morgan Kaufmann
