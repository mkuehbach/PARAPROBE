**Setup**
=========

Which operating system is supported?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PARAPROBE is a  high performance computing (HPC) back end solution for processing APT datasets. Therefore, it targets workstations and computing clusters, i.e. Linux-based operation systems. The compilation on a Windows system should in principle be technically possible, has so far, though, not been tested.

How large datasets are supported?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Currently, single APT measurements of at most 4.2 billion ions technically. Dataset sizes of 2.0 billion ions were tested thoroughly. As of 2018, such successful tip measurement are to the best of my knowledge not standard. Please contact me_ if you have larger datasets, I am eager to modify my code to become capable as well to handle even such higher ion counts per single measurement.

 .. _me: https://www.mpie.de/person/51206/2656491

What are the minimum hardware requirements?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Memory** --- data mining APT datasets is 3d point data processing. Therefore, hardware minimum requirements depend primarily and necessarily on the total number of ions. Sufficient system main memory is required to hold the point data and temporary partial duplicates of it during processing. Internally, each ion is represented as a structure of three 32-bit floating point numbers surplus one 32-bit unsigned integer, hence requiring 16B per ion. Quantitative results are detailed in the initial PARAPROBE paper (see Reference section).

**CPU** --- virtually all modern workstation and cluster computing processors are capable of executing PARAPROBE. The cost-benefit-ratio and speed of doing so may differ substantially so. Consequently, claiming minimum hardware requirements is pointless. Quantitative results (see Reference section) document better than 50% strong scaling efficiency for up to 36 threads for all analysis tasks except executing the DBScan algorithm.

**GPU** --- PARAPROBE currently does not utilize GPU parallelism.

Which prerequisites are necessary?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PARAPROBE depends on third-party open source software and open source Linux tools. Please follow these first steps to assure you have a working system 
that is capable of compiling the PARAPROBE source code, link to the libraries required, and execute.

* Check for a working installation of a **C/C++ build system** including **cmake** and **make**.
* You need a working installation of the **Boost C++ libraries**. Further details about Boost_ here.
* A default installation of Ubuntu in at least version 17.10 provided me with the above-mentioned prerequisites.
* For Ubuntu 16.04 LTS an installation of the newest Boost version is necessary.

* PARAPROBE has been tested to compile with the GNU and the Intel Parallel Studio 2018 compiler.
* Given the fact that PARAPROBE uses the IntelMKL library, using it out of the box demands to use the Intel compiler.

* You need a working installation of an **MPI_ (Message Passing Interface)** API_ library.
* The minimum threading support level of the MPI implementation required is **MPI_THREAD_FUNNELED**.
* MPI libraries are not installed by default. They are available here MPICH_ . It is recommended to use the IntelMPI library.

* PARAPROBE uses the **Computational Geometry Algorithms Library (CGAL)**. It has own prerequisites.
* At least the two arbitrary precision arithmetic libraries **GMP_** and **MPFR_**.
* Thus, the first step is to prepare these libraries accordingly. 
* CGAL_ can be downloaded CGAL here_ 
* Next, it suffices to activate the header-only mode of the CGAL library. Do so by modifying the line 
  containing *enable cgal header only* in the CMakeLists.txt file topmost in the CGAL code folder.
* Using the header-only library worked for me with both the CGAL version 4.11.3 and 4.12. It failed so far for 4.13.
* Within the top level CGAL code folder configure once the library by typing::

   cmake -DCGAL_DIR=<CGALLocation> .

* Please note that the angle brackets must not be typed into the command line as they only mark the input parameter!

* The PARAPROBE CMakeLists.txt should now be able to use the CGAL functionalities.
   
* PARAPROBE utilitzes the **Hierarchical Data Format (HDF5)** library and the **eXtensible Data Model and Format** XDMF_
* Personally, I use a local installation of the HDF5_ library. This worked for me with version 1.10.2. 
* I recommend to install the HDFViewer_ a tool for looking into the binary content of an H5 file.

* A local installation of HDF5 worked for me using the following procedure using version 1.10.2
* Starting from the PARAPROBE top level folder the Github repository contains a copy of a HDF5 source code tar archive, unpack it::

   cd src/thirdparty/HDF5/CMake-hdf5-1.10.2.tar.gz
   tar -xvf CMake-hdf5-1.10.2.tar.gz

* Make a local compile folder and make sure your environment has the compiler and MPI you want to use::

   mkdir build
   cd build

* Configure a compilation script which inspects the technical details of your system to configure HDF5::

   cmake -G "Unix Makesfiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=ON -DHDF5_BUILD_TOOLS=ON -DHDF5_BUILD_FORTRAN=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF -DHDF5_ENABLE_SZIP_ENCODING=OFF ../hdf5-1.10.2

* Compile the library from the source code::
   
   cmake --build . --config Release

* Test it on the system::

   ctest . -C Release

* Pack it into an archive and complete the setting up locally::

   cpack -C Release CPackConfig.cmake
   ./HDF5-1.10.2-Linux.sh

* For computing tessellations PARAPROBE builds on **Chris Rycroft's Voro++** (Voro_). Its source code is compilation ready within src/thirdparty/voro.
* Atom probe crystallography analyses demand discrete Fourier transform algorithms. For this PARAPROBE utilizes the Intel Math Kernel Library (IMKL_).
* For many research purposes and students the library is open source (IMKL_). 

 .. _Boost: https://www.boost.org/
 .. _MPI: https://www.mcs.anl.gov/research/projects/mpi/
 .. _API: https://www.mpich.org/downloads/
 .. _MPICH: https://www.mpich.org/downloads
 .. _Intel: https://software.intel.com/en-us/intel-mpi-library
 .. _Boost: https://www.boost.org/
 .. _IMKL: https://software.intel.com/en-us/performance-libraries
 .. _GMP: https://gmplib.org/
 .. _MPFR: https://www.mpfr.org/
 .. _CGAL: https://doc.cgal.org/latest/Manual/installation.html
 .. _here: https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.12/CGAL-4.12.tar.xz
 .. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
 .. _HDFViewer: https://www.hdfgroup.org/downloads/hdfview/
 .. _XDMF: https://www.xdmf.org/index.php/Main_Page
 .. _Voro: https://math.lbl.gov/voro++/

How to compile?
^^^^^^^^^^^^^^^
Once all prerequisites are met, proceed to configure and compile PARAPROBE.

* Download the source from its git repository **https://github.com/mkuehbach/PARAPROBE**
* Unpack the repository such that finally the following ends up in a single folder, from now on referred to as the **root** folder. 
    * a **src** subdirectory with the cpp and the h source code files,
    * a **thirdparty** subdirectory with a compile-ready **RapidXML**, **CGAL**, **Voro**, **HDF5** 
    * a **build** directory for storing the executable
    * a **XML control file**.
    * a **scripts** subdirectory with useful tools for processing PARAPROBE results further.
    * Additionally, check that there is a **CMakeLists.txt** file in the root folder.
* You can now rename, if you desire, the root folder to any Linux-conformant name.
* Next, utilize the top section in CMakeList.txt file to choose compiler and specify the paths as detailed.
* Next, open a console and dive into the **build** directory.
* If you now compile PARAPROBE for the first time type::

   cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=<CGALLocation> ..
   
* Replace <CGALLocation> by the string that specifies the absolute path where the CGAL code folder is on your system.
* Now cmake inspects your system configuration, finds compilers, libraries, which eventually results in a customized **Makefile**.
* Next, or if compiling not for the first time, use this makefile by initiate the compilation process::

   make
   
* **Warnings** will appear but can be ignored.
* Upon success, you should now have a PARAPROBE **executable** with the name as specified in the CMakeLists.txt within the build.
* Use this executable to perform APT post-processing. Always a **XML** control file, a **RRNG** rangefile, and eventually **POS**, **EPOS**, or **APT** measurement raw data file is necessary.

Where to place files?
^^^^^^^^^^^^^^^^^^^^^
The resulting executable expects the XML control file always in its current location folder! Relative indexing is utilized. Other than that restriction, the executable can be renamed and relocated. This enables to script batch queues for PARAPROBE.

Optimization
^^^^^^^^^^^^
If desired, adjust the level of compiler optimization via the OPTLEVEL variable in the CMakeLists.txt upper section. 
OPTLEVEL "-O0" means no optimization and should be utilized for debugging purposes only, while "-O3 -march=native" is the maximum and recommended level for production tasks. Such highly compile time optimized code is not necessarily portable. 
Improvements between the two extremes vary between a factor of 2 - 5 faster with maximum optimization compared to without.

Troubleshooting?!
^^^^^^^^^^^^^^^^^
If unrecoverable errors occur during the compilation process, attempt first to instruct a **make clean** command. This will delete potentially incompletely processed source code files. If this does not help: delete everything in the build folder except for the **XML control file** and start over with **cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=<CGALLocation> ..**. 

.. Optional prerequisites
.. ^^^^^^^^^^^^^^^^^^^^^^
.. * It is well-known that the general purpose standard malloc(3) memory allocator class implementation does not assure that in particular many small allocations become locality-aware placed. 
.. * Hence, the performance of PARAPROBE can be improved by linking against an alternative memory allocator class, such as **Jason Evans jemalloc**.
.. * A documentation of how to obtain the library and how to compile it can be found online_.

..  .. _online: https://www.canonware.com/jemalloc/
