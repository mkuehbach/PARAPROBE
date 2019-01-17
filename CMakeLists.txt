cmake_minimum_required(VERSION 3.0)

################################################################################################################################
##USER INTERACTION##############################################################################################################
##in this section modifications and local paths need to be modified by the user#################################################
################################################################################################################################
#please name your project accordingly
set(MYPROJECTNAME "paraprobe")

#tell the top directory where this local PARAPROBE version is stored
set(MYPROJECTPATH "/home/m.kuehbach/PARAPROBE_Paper14_FINAL/")

#COMPILER, choose which one to use
##GNU 7.3 and Intel 18.1 in combination with MPICH and IntelMPI have been tested
set(EMPLOY_INTELCOMPILER ON) #either or
set(EMPLOY_GNUCOMPILER OFF)

#choose optimization level
##-O0 nothing, debug purposes, -O1 moderate optimization, -O2 -O3 for production level up to aggressive optimization
set(MYOPTLEVEL "-O3 -march=native")

#necessary libraries
#CGAL with internal dependencies gmp mfpr and Boost necessary
#you have to configurate the compilation by typing cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=YOURPATHTOCGAL before
#here the placeholder YOURPATHTOCGAL should be replaced with the absolute path were your local CGAL is
set(EMPLOY_CGAL ON)
set(EMPLOY_MYINTELMKL ON)
set(EMPLOY_MYHDF5 ON)
set(EMPLOY_BOOST ON)

#optional libraries
set(EMPLOY_BSIMDSUPPORT OFF)
set(EMPLOY_JEMALLOC OFF)

################################################################################################################################
#END OF INTERACTION FOR NON PRO USERS###########################################################################################
##in this section advanced users might want/need to make modifications if they use non default places for thirdparty libraries##
################################################################################################################################ 

#INTEL MATH KERNEL LIBRARY threaded numerical library Intel MKL
#manual finding of IntelMKL where do we find the IntelMKL library
set(MYMKLROOT "/opt/intel/compilers_and_libraries_2018/linux/mkl//")
set(MYMKLPATH "${MYMKLROOT}/lib/intel64_lin")
include_directories("${MYMKLROOT}/include")
set(MYMKL_COMPILEFLAGS "")
set(MYMKL_LINKFLAGS "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
#set(MYMKL_LINKFLAGS "-mkl -DMKL_ILP64")


#HDF5 local installation for advanced I/O, collecting metadata and bundle analysis results together
set(MYHDFPATH "${MYPROJECTPATH}/src/thirdparty/HDF5/CMake-hdf5-1.10.2/build/HDF5-1.10.2-Linux/HDF_Group/HDF5/1.10.2")
#download Unix tar archive for the HDF5-1.10.2 release from www.hdf5.org
#copy into local directory, do not use sudo to execute
#unpack via tar -xvf tarballname
#change into directory from now on this is referred to as the hdfroot folder
#execute ./build_unix.sh  this will configure the static libraries, sequential no fortran support and use deflate compression library
#when configure build and test has been completed
#take packed archive and extract into the automatically generated build folder within the hdfroot folder
#deep in these folders are the include directories and the static libraries (a files) we to understand the hdf5.cpp commands and 
#the libraries we need to link statically against to execute hdf5 functionalities
include_directories("${MYHDFPATH}/include")
link_directories("${MYHDFPATH}/lib")
##link against static libraries
set(MYHDFLINKFLAGS "-L${MYHDFPATH}/lib/ ${MYHDFPATH}/lib/libhdf5_hl.a ${MYHDFPATH}/lib/libhdf5.a -lz -ldl")

##########################################################################################
#####PERFORMANCE RELATED##################################################################
##here settings for scalable vectorization and scalable allocator classes are made########
##so far these components are optional####################################################


#NumScale bSIMD scalable vectorization support
if (EMPLOY_BSIMDSUPPORT)
	include_directories("${MYPROJECTPATH}/src/thirdparty/bSIMD/home/m.kuehbach/SIMD/master/include")
	set(MYSIMDFLAGS "-msse4.2")
elseif()
	set(MYSIMDFLAGS "")
endif()

#Linux NUMA Non-uniform memory access (NUMA)
#might need installation on cluster
##needs necessarily linking into the program but not necessarily used
set(MYNUMA "-lnuma")
add_definitions(${MYNUMA})

#Scalable allocation using Jason Evans alternative malloc3 implementation##################
##alternatives are for instance tcmmalloc##################################################
##optional: if desired specify location of working jemalloc installation or local compiled version of the library
##use only in combination with pinned threads to reduce fragmentation and improve memory locality during allocation
if (EMPLOY_JEMALLOC)
	set(MYJEMALLOCPATH "${MYPROJECTPATH}/src/thirdparty/Jemalloc/jemalloc-5.1.0")
	set(MYSCALABLEALLOC "-ljemalloc")
	add_definitions(${MYSCALABLEALLOC})
	include_directories("${MYJEMALLOCPATH}/include")
	link_directories("${MYJEMALLOCPATH}/lib")
elseif()
	set(MYSCALABLEALLOC "")
	add_definitions(${MYSCALABLEALLOC})
endif()


#additional hints IntelMKL but currently not utilized
#set(MYMKL "-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
#set(MYMKL "-DMKL_ILP64 -parallel -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl")
#set(MYMKL_COMPILEFLAGS "-mkl -DMKL_ILP64")

#check also https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#works for now
#set(MYMKL_COMPILEFLAGS "")
#set(MYMKL_LINKFLAGS "-parallel -Wl,--start-group ${MYMKLPATH}/libmkl_intel_lp64.a ${MYMKLPATH}/libmkl_sequential.a ${MYMKLPATH}/libmkl_core.a -Wl,--end-group -lpthread")


#############################################################################################################################
##AUTOMATIC SECTION##########################################################################################################
##normally no further changes should be required below unless for development################################################
#############################################################################################################################
set(VERBOSE_VECTORIZATION ON)

#user input sanity checks
if(EMPLOY_CGAL)
	message([STATUS] "We use the CGAL library for computational geometry")
else()
	message([FATAL_ERROR] "It is necessary to compile in at least the header only CGAL library")
endif()
if(EMPLOY_MYINTELMKL)
	message([STATUS] "We use the Intel Math Kernel library (IMKL) for discrete Fourier transform functionality")
else()
	message([FATAL_ERROR] "It is necessary to link against the IntelMKL library")
endif()
if(EMPLOY_MYHDF5)
	message([STATUS] "We use the HDF5 library for advanced I/O")
else()
	message([FATAL_ERROR] "It is necessary to link against the HDF5 library")
endif()

	
if(EMPLOY_INTELCOMPILER AND EMPLOY_GNUCOMPILER)
	message([FATAL_ERROR] "You cannot utilize two compiler at the same time!")
endif()


#automatically assign project name and compiler flags
project(${MYPROJECTNAME})
set(CMAKE_BUILD_DIR "build")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYOPTLEVEL}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYOPTLEVEL}")

#setting up compiler-specifics
#intel path
if(EMPLOY_INTELCOMPILER)
	message([STATUS] "Employing the Intel compiler!")
	if (VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-vec-report=2") #https://software.intel.com/en-us/node/522949/ "-qopt-report-phase=vec")
	else()
		set(MYVERBOSE "")
	endif()
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE}")
	message([WARNING] "Currently not SIMD flags provided for Intel compile options")
	set(MYOMP "-qopenmp -lpthread")
	add_definitions(${MYOMP})
	add_definitions("-std=c++0x")
	add_definitions("-Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare 
		-Wsequence-point -Wtrigraphs -Wunused-variable")
#gcc path
elseif(EMPLOY_GNUCOMPILER)
	message([STATUS] "Employing the GNU compiler!")
	if(VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-fopt-info")
	endif()
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE} ${MYSIMDFLAGS}")
	set(MYOMP "-fopenmp -lpthread")
	add_definitions(${MYOMP})
	add_definitions("-std=c++11")
	add_definitions("-Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point 
		-Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable")
else()
	message([FATAL_ERROR] "You have to utilize a compiler!")
endif()

message([STATUS] "Projectname is ${MYPROJECTNAME}")
message([STATUS] "We utilize optimization level ${MYOPTLEVEL}")

#parallelization - MPI process-level
#query location of MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(MYSRCPATH "${MYPROJECTPATH}/src/")
set(MYVOROSRCPATH "${MYPROJECTPATH}/src/thirdparty/VoroRycroft/voro++-0.4.6/src/")

add_executable(${MYPROJECTNAME}
	${MYSRCPATH}PARAPROBE_IntelMKL.cpp
	${MYSRCPATH}PARAPROBE_NUMABinding.cpp
	${MYSRCPATH}PARAPROBE_ScalableAllocator.cpp
	${MYSRCPATH}PARAPROBE_Datatypes.cpp
	${MYSRCPATH}PARAPROBE_PeriodicTable.cpp
	${MYSRCPATH}PARAPROBE_Histogram.cpp
	${MYSRCPATH}PARAPROBE_DiscreteHistogram.cpp
	${MYSRCPATH}PARAPROBE_NPointCorrelations3D.cpp
	${MYSRCPATH}PARAPROBE_Math.cpp
	${MYSRCPATH}PARAPROBE_Settings.cpp
	${MYSRCPATH}PARAPROBE_AABBTree.cpp
	${MYSRCPATH}PARAPROBE_SpaceBucketing.cpp
	${MYSRCPATH}PARAPROBE_KDTree.cpp
	${MYSRCPATH}PARAPROBE_HoshenKopelman.cpp
	${MYSRCPATH}PARAPROBE_EPOSEndianness.cpp
	${MYSRCPATH}PARAPROBE_VTKIO.cpp
	${MYSRCPATH}PARAPROBE_HDF5.cpp
	${MYSRCPATH}PARAPROBE_XDMF.cpp
	${MYSRCPATH}PARAPROBE_TAPSimHdl.cpp
	${MYSRCPATH}PARAPROBE_Profiling.cpp
	${MYSRCPATH}PARAPROBE_HPDBScan.cpp
	
	${MYVOROSRCPATH}c_loops.cc
	${MYVOROSRCPATH}cell.cc
	#${MYVOROSRCPATH}cmd_line.cc
	${MYVOROSRCPATH}common.cc
	${MYVOROSRCPATH}container.cc
	${MYVOROSRCPATH}container_prd.cc
	${MYVOROSRCPATH}pre_container.cc
	${MYVOROSRCPATH}unitcell.cc
	${MYVOROSRCPATH}v_base.cc
	#${MYVOROSRCPATH}v_base_wl.cc
	${MYVOROSRCPATH}v_compute.cc
	#${MYVOROSRCPATH}voro++.cc
	${MYVOROSRCPATH}wall.cc
	
	${MYSRCPATH}PARAPROBE_VoroXX.cpp
	${MYSRCPATH}PARAPROBE_SolverHdl.cpp
	${MYSRCPATH}PARAPROBE_Main.cpp
)

#tell us where to find the header only CGAL library paths
##GMP and MPFR is included and found alongside the CGAL library
message([STATUS] "We need the header-only CGAL library...")
find_package(CGAL QUIET)
if( CGAL_FOUND )
	include( ${CGAL_USE_FILE} )
endif()


##assure the Boost package is there to allow and ease the reading of files and folders from the filesystem
message([STATUS] "We need Boost as well...")
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
find_package(Boost COMPONENTS system filesystem REQUIRED)


#linking process
##mind the y and dl libraries at the end which are required for gcc support when using the hdf libraries
target_link_libraries(${MYPROJECTNAME} ${MYOMP} -lgmp -lm ${MYNUMA} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${MYHDFLINKFLAGS} ${MYSCALABLEALLOC} ${MYMKL_LINKFLAGS})

#MPI compilation settings
if(MPI_COMPILE_FLAGS)
  set_target_properties(${MYPROJECTNAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${MYPROJECTNAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
