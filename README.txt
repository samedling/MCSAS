Basic Installation:
1. Make sure you have the latest version:
   Go to https://github.com/samedling/MCSAS and click Download ZIP in the lower right.
   Or a current zipfile may be at https://drive.google.com/open?id=0B8EbmzXGZtaZV3MxUlhSUjczQm8
2. If possible, run `make` to compile the fortran code to achieve 2-100x speedup (depending on number of cores).  Or, try copying the fastmath-OS_CPU.so file to fastmath.so.

Operation:
Run `python newgui.py`.


Tested using gfortran on Ubuntu14.10, ifort on a raijin login node (no OpenMP), and both on CentOS.  Default is gfortran; to use ifort simply edit the makefile.


OS X Fortran Installation:
Apple doesn't provide a recent version of gfortran, but you can download one from http://hpc.sourceforge.net
If you also don't have Xcode Tools, follow the directions here:
https://wiki.helsinki.fi/display/HUGG/Installing+the+GNU+compilers+on+Mac+OS+X
Basically:
1. Install XCode Tools from the App Store.
2. Install the Command Line Tools by running xcode-select --install
3. Download the latest stable gfortran version from http://hpc.sourceforge.net

Tested on OS X 10.10 "Yosemite".


Some Troubleshooting:
OS X/EPD/Tkinter: Make sure you have Canopy.  EPD might tell you it's updated everything, but it's still not the same as Canopy.
Centos/F2PY: Make sure you are using the version of F2PY which matches your version of Python (2.7+), otherwise just loading the Fortran module causes a Segmentation Fault.
OS X (and Windows?): You need to close all the old plots before you can run things again.  Otherwise, it's otherwise unresponsive for some reason.
