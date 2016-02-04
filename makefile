CC=f2py     #general needed for OSX
#CC=f2py2.7  #specific needed for scucomp

linux_so:
	$(CC) -c --fcompiler=gnu95 --opt="-O3 -fopenmp" -lgomp -m fastmath fastmath.f90
	#f2py -c --fcompiler=intelem --opt="-O3 -openmp" -m fastmath fastmath.f90
	#/usr/local/epd/bin/f2py -c --fcompiler=intelem --opt="-O3 -openmp" -m fastmath fastmath.f90

no_openmp:
	f2py --opt="-O3" -c -m fastmath fastmath.f90

generic:
	f2py --noarch -c -m fastmath fastmath.f90

uncompiled_module:
	f2py -m fastmath fastmath.f90

sig_pyf:
	f2py -m fastmath -h fastmath.pyf fastmath.f90

compile_w_pyf: sig_pyf
	f2py -c fastmath.pyf fastmath.f90

all: linux_so uncompiled_module sig_pyf

clean:
	rm -f fastmath.so fastmath.pyf fastmathmodule.c
