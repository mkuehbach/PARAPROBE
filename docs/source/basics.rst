**In a nutshell**
=================

What is PARAPROBE?
^^^^^^^^^^^^^^^^^^
| A software for data mining Atom Probe Tomography (APT) experiment data which sets
| prime focus on scalable spatial range querying and computational geometry tasks
| making use of scalable hierarchical parallelism.

What are the user benefits?
^^^^^^^^^^^^^^^^^^^^^^^^^^^
| **Unbiased descriptive spatial statistics**
|	Enabled by state of the art tip surface reconstruction surplus ion to surface distancing.
| **Scalable performance**
| 	Owing to parallel implementation with rigorous hierarchical spatial data partitioning
| 	strategy to improve the utilization of fast caches and ccNUMA-aware data placement policy.
| **Open source software**
|	Therefore no usage restriction, no limited licences when running in parallel
|	surplus full functional transparency and modifiability.

Which parallelization concepts are employed?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Process data parallelism** --- the Message passing interface (MPI_) is utilized.
It processes each individual measurement via a single process. This enables to either distribute parameter sweeping studies of the same tip on practically hundred thousands of processes or to process trivially in parallel multiple tips using the same automatized analysis protocol. At runtime, MPI invokes library calls to communicate pieces of information between physically disjoint computers. Therefore, it requires installation and linking.

**Shared memory thread data parallelism** --- the Open Multiprocessing (OpenMP_) is used.
It allows  to distribute the point data of each measurement spatially into disjoint logical chunks. 
These chunks are mapped to and stored in thread-local memory when processed in parallel.
For some tasks the threads update explicitly data in the memory of other threads processing spatially adjacent points. 
In such cases, explicit care is taken to prevent data invalidation and reduce false sharing.
OpenMP does not require installation because it builds on preprocessor directives, which get evaluated during the compilation stage.

**Streaming instruction data parallelism aka SIMD** --- the BSIMD_ portable vector intrinsics template library is used.
Some core geometrical and numerical tasks can be accelerated further within each thread using vectorization.
Such vectorization is realized via Single Instruction Multiple Data (SIMD) which makes use of highly problem-and-CPU-specific 
instructions, the so-called intrinsics. Their key idea is to apply processing operations on a packet of multiple data elements 
at once instead of sequentially. PARAPROBE employs BSIMD_ in order to improve code portability by solving the problem that 
intrinsics have usually different names and implementation syntax for different CPUs. Upon compilation, the abstract BSIMD intrinsics 
commands are encoded into the specific CPU command realization available on the target architecture.

 .. _MPI: http://www.mcs.anl.gov/research/projects/mpi/
 .. _OpenMP: https://www.openmp.org/
 .. _BSIMD: https://developer.numscale.com/bsimd/documentation/v1.17.6.0/faq.html

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
| High Performance Parallelism Pearls Volume One: Multicore and Many-Core Programming Approaches
| 1st edition, 2014, Morgan Kaufmann

| J. Jeffers and J. Reinders
| High Performance Parallelism Pearls Volume Two: Multicore and Many-Core Programming Approaches
| 1st edition, 2015, Morgan Kaufmann
