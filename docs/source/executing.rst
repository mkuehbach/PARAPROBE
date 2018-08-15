**Program execution**
=====================

How to execute
^^^^^^^^^^^^^^	
All set? Excellent! Then, you have to prepare your environment to allow thread parallelism.
Specifically set the OMP_NUM_THREADS environment variable. If sequential execution desired use 1, else
use a number as high as the number of physical cores, counting hyper-threading core pairs with one core, here
exemplified for a workstation with 36 physical cores::

   export OMP_NUM_THREADS=20
   
Next, modify further your environment to allocate sufficient stack space for the application::

   export ulimit -s unlimited

Now PARAPROBE can be executed via a single command line call::

   mpiexec -n <nprocesses> <paraprobe> <simid> <Settings.xml>
   
Or equivalently::

   mpirun -np <nprocesses> <paraprobe> <simid> <Settings.xml>
   
**Please note that the angle brackets must not be typed into the command line as they only mark the input parameter!**
   
In its current version the following input arguments are required:

* <nprocesses> How many MPI processes to utilize?
* <paraprobe> Your specific name of the executable.
* <simid> JobID, an unsigned integer to distinguish the results from runs with other settings but the same raw data. 
* <Settings.xml> a properly formatted XML control file. The name can be changed as long as the file remains a properly formatted XML file.

**Be careful: if the <simid> value is set to the same value during subsequent runs in the same folder, data will be overwritten without prompting!**

PARAPROBE runs always threaded, i.e. with at least one OpenMP thread. More threads can be used to activate parallelism to speed up the processing.

Report runtime diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^
In particular for debugging and getting to know further information how PARAPROBE performed it is useful to store the console prompt during execution.
In order to do so execute the program as follows::

   <regular execution command as shown above> 2>&1 | tee PARAPROBE.Console.STDOUTERR.txt
   
This will instruct to redirect the console output and error verbose into a text file surplus show the results as usual on the console.
Instead::

   <regular execution command as shown above> 1>PARAPROBE.Console.STDOUT.txt 2>PARAPROBE.Console.STDERR.txt
   
will redirect all verbose to separate text files. One for usual output STDOUT. The other one for operating system controlled errors returned by the program STDERR.

Benchmarking
^^^^^^^^^^^^

PARAPROBE has internal functionality to monitor its elapsed time expenditures. After running successfully, this is summarized in the MyProfiling.csv file.

Further details of OpenMP multi-threading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Below is a typical example call which executes the program with 1 MPI process spawning 20 OpenMP threads, reads the settings from MySpecialSettings and marks all results with a consistent ID, 1000 in this case::

   mpiexec -n 1 paraprobe 1000 MySpecialSettings.xml
   
Please note, setting the environment variable OMP_NUM_THREADS is required only once per active console session, if not working with a different number of threads is desired.
Please note as well that it is when using Intel CPUs, Hyper-Threading (HT_) cores suggest that twice as many cores could be used as in practice 
should be. The reason is that a hyper-threading core pair is build of two cores which share most of the fastest caches. In effect, the two
hyper-threading cores will fight for resources when assign individual work packages. For this reasoning, only as many threads as 
real physical core **pairs* exist should be instructed using the OMP_NUM_THREADS environment variable.

 .. _HT: https://www.intel.com/content/www/us/en/architecture-and-technology/hyper-threading/hyper-threading-technology.html

Manage thread affinity
^^^^^^^^^^^^^^^^^^^^^^
Thread affinity deals with the question whether a thread, once it has been mapped on a CPU core, should always be executed on this core or may
be relocated to another in an effort to improve global system performance. In fact, if not explicitly instructed otherwise, 
the operating system usually instructs such migration operations during program operation, in particular when multiple user or programs
work on the same system. Such strategy to improve global system performance, however, may not necessarily mean a performance improvement of the
actual running PARAPROBE job. Moreover, frequent thread migration reduces PARAPROBE performance owing to re-initialization costs. Namely, once a thread is instructed to migrate, its already allocated cache content will in most cases not exist on the new core. Instead, costly refurnishing of cache content is required. For these reasons, overwriting the operating systems thread affinity policy can be worthwhile to prevent
the OS from migrating threads across the available cores. One can instruct this explicitly via setting for instance the KMP AFFINITY_ 
environment variables when using the Intel compiler or the GNU_ thread affinity policy.

 .. _AFFINITY: https://software.intel.com/en-us/node/522691
 .. _GNU: https://gcc.gnu.org/onlinedocs/libgomp/GOMP_005fCPU_005fAFFINITY.html
