Basic Installation:
1. Edit the first line of newgui.py to point to your python installation.
2. If possible, run `make` to compile the fortran code to achieve 2-100x speedup (depending on number of cores).

Operation:
3. `Run ./newgui.py` each time you want to start the program.

Tested using gfortran on Ubuntu14.10, ifort on a raijin login node (no OpenMP), and both on CentOS.  Default is gfortran; to use ifort simply edit the makefile.
(FYI, I had some difficulty with multiple F2PY versions on CentOS such that you need to make sure you are using the version of F2PY which matches your version of Python (2.7+), otherwise just loading the Fortran module causes a Segmentation Fault.)


OS X Fortran Installation:
Apple doesn't provide a recent version of gfortran, but you can download one from http://hpc.sourceforge.net
If you also don't have Xcode Tools, follow the directions here:
https://wiki.helsinki.fi/display/HUGG/Installing+the+GNU+compilers+on+Mac+OS+X
Basically:
1. Install XCode Tools from the App Store.
2. Install the Command Line Tools by running xcode-select --install
3. Download the latest stable gfortran version from http://hpc.sourceforge.net

Note that on OS X you need to close all the old plots before you can run things again.  Otherwise, it's otherwise unresponsive for some reason.

Tested on OS X 10.10 "Yosemite".

