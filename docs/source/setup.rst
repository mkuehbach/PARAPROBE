**Setup**
=========

Which operating system is supported?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PARAPROBE is envisioned as a back-end high performance computing (HPC) solution for processing APT datasets. 
Therefore, it targets workstations and computing clusters running Linux. The compilation on a Windows system should in 
principle be technically possible, has so far, though, not been tested.

How large datasets are supported?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Currently, single APT measurements of at most 4.2 billion ions are supported. As of 2018, such successful tip measurement are
to the best of my knowledge not standard. Please contact me_ if you have larger datasets, I am eager to modify my code to become 
capable as well to handle even such higher ion counts per single measurement.

 .. _me: https://www.mpie.de/person/51206/2656491

What are the minimum hardware requirements?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Memory** --- data mining APT datasets is 3d point data processing. Therefore, hardware minimum requirements depend primarily and necessarily on the total number of ions. Sufficient system main memory is required to hold the point data and temporary partial duplicates of it during processing. Internally, each ion is represented as a structure of three 32-bit floating point numbers surplus one 32-bit unsigned integer, hence requiring 16B per ion.

 .. **Empirical overhead factors** are as follows

 ..   .. role:: red

 .. |	Bas et al. reconstruction
 .. |	Spatial range querying indexing structure
 .. |	Tip surface extraction

 .. Consequently, conducting nearest neighbor analyses requires XXXX memory.

**CPU** --- virtually all modern workstation and cluster computing processors are capable of executing PARAPROBE, yet their cost-benefit-ratio and speed of doing so may differ substantially. Consequently, claiming minimum hardware requirements is pointless. 

 .. | **MORE ADVICE SHOULD BE GIVEN HERE**

**GPU** --- PARAPROBE currently does not utilize GPU parallelism.

Which prerequisites are necessary?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| PARAPROBE depends on third-party open source software and open source Linux tools. 
| Please follow these first steps to assure you have a working system 
| and able to compile and execute PARAPROBE.

* Check for a working installation of a **C/C++ build system** including **cmake** and **make**.
* You need a working installation of the **Boost C++ libraries**. Further details about Boost_ here.
* A default installation of Ubuntu in at least version 17.10 provided me with the above-mentioned prerequisites.
* PARAPROBE has been tested to compile with the GNU and the Intel Parallel Studio 2018 compiler.

* You need a working installation of an MPI_ (Message Passing Interface) API_ library.
* The minimum threading support level of the MPI implementation required is **MPI_THREAD_FUNNELED**.
* MPI libraries are not installed by default. They are available here MPICH_ .

* PARAPROBE uses the **Computational Geometry Algorithm Library** (CGAL_).
* This library requires own prerequisites, at least the two arbitrary precision arithmetic libraries GMP_ and MPFR_ .
* Thus, the first step is to prepare these libraries accordingly. 
* CGAL can be downloaded CGAL here_ 
* Next, it is essential to activate the HEADERONLY library mode within the CGAL. Do so by modifying the line containing
  *enable cgal header only* in the CMakeLists.txt document topmost in the CGAL code folder
* Within the toplevel of the CGAL code folder configure once the library by typing cmake .
 
* The PARAPROBE CMakeLists.txt should now be able to use the CGAL functionalities.

 .. * ADD INSTALLATION DETAILS ABOUT CGAL ADD IMAGE HOW TO CHANGE CGAL CMAKE SCRIPT

 .. _MPI: https://www.mcs.anl.gov/research/projects/mpi/
 .. _API: https://www.mpich.org/downloads/
 .. _MPICH: https://www.mpich.org/downloads
 .. _Intel: https://software.intel.com/en-us/intel-mpi-library
 .. _Boost: https://www.boost.org/
 .. _GMP: https://gmplib.org/
 .. _MPFR: https://www.mpfr.org/
 .. _CGAL: https://doc.cgal.org/latest/Manual/installation.html
 .. _here: https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.12/CGAL-4.12.tar.xz


How to compile?
^^^^^^^^^^^^^^^
Once all prerequisites are met, proceed to configure and compile PARAPROBE.

* Download the source from its git repository **https://github.com/mkuehbach/PARAPROBE**
* Unpack the repository such that finally the following ends up in a single folder, from now on referred to as the **root** folder. 
    * a **src** subdirectory with the cpp and the h source code files, 
    * a **build** directory for storing the executable
    * a **XML control file**.
    * Additionally, check that there is a **CMakeLists.txt** file in the root folder.
* You can now rename, if you desire, the root folder to any Linux-conformant name.

* Next, utilize the top section in CMakeList.txt file to choose compiler and specify paths to local Boost and SIMD as well as switch
  on options and give location of HDF5 library if it should be used if in doubt, use OFF. bSIMD and HDF5 are not essential to the program.
* Next, open a console and dive into the **build** directory.
* If you now compile PARAPROBE for the first time type **cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=LOCALPATH ..**. 
* Replace LOCALPATH by the absolute path where the CGAL code folder is on your system.
* Now cmake inspects your system configuration to find compilers and libraries and generates a customized makefile for you.
* Next, or if compiling not for the first time, use this makefile by typing **make** to start the compilation process.
* Warnings about declared but unreferenced variables will appear they can be ignored.
* Upon success, you should now have an **executable** with the name as specified in the CMakeLists.txt within the build.
* Use this executable to perform APT post-processing.

Where to place files?
^^^^^^^^^^^^^^^^^^^^^
| The resulting executable expects the XML control file always in its current location folder!
| Relative indexing is utilized. Other than that restriction, the executable can be renamed and relocated. 
| The latter enables PARAPROBE batch queue processing.

Optimization
^^^^^^^^^^^^
If desired, adjust the level of compiler optimization via the OPTLEVEL variable in the CMakeLists.txt upper section. 
OPTLEVEL "-O0" means no optimization and should be utilized for debugging purposes only, while "-O3" is the maximum and recommended level for production tasks. Improvements between the two extremes vary between a factor of 2 - 5 faster with maximum optimization compared to without.

Troubleshooting?!
^^^^^^^^^^^^^^^^^
If in between the compilation process unrecoverable errors occur, attempt first a **make clean** command. 
If this does not help: Delete everything in the build folder except for the **XML control file** and start over with **cmake -DCMAKE_BUILD_TYPE=Release ..**. 

I faced challenges when attempting to compile the CGAL in header-only mode with the Intel Parallel Suite compilers and tools.
Therefore, I recommend so far to use the GNU compiler and the MPICH installation.


.. Optional prerequisites
.. ^^^^^^^^^^^^^^^^^^^^^^
.. * It is well-known that the general purpose standard malloc(3) memory allocator class implementation does not assure that in particular many small allocations become locality-aware placed. 
.. * Hence, the performance of PARAPROBE can be improved by linking against an alternative memory allocator class, such as **Jason Evans jemalloc**.
.. * A documentation of how to obtain the library and how to compile it can be found online_.
.. * In the future, an OpenMP-parallelization and HDF5 extension of the TopologyTracer is planned. So far however HDF5 is not required to run the TopologyTracer.

..  .. _online: https://www.canonware.com/jemalloc/
