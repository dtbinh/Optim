####################################################################
#  SLICOT make include file for Intel Fortran compiler.            #
#  SLICOT, Release 5.0                          .\slicot\make.inc  #
#  Vasile Sima, KU Leuven                                          #
#  October 31, 1996.                                               #
#  Revised December 23, 1999; January 2005; February 2009.         #
####################################################################
#
#  The root directory for the SLICOT files
#
HOME = \slicot
PLAT = _win32
#
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = ifort
OPTS     = /O2 /fp:source
NOOPT    = -optimize:0
LOADER   = $(FORTRAN)
LOADOPTS = $(SLICOTLIB) $(LPKAUXLIB) $(LAPACKLIB)
#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH     = lib
ARCHFLAGS= -out
#ARCHFLAGS= 
RANLIB   = echo
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
LPKAUXLIB    = ..\lpkaux.lib
LAPACKLIB    = mkl_lapack.lib mkl_intel_c.lib mkl_sequential.lib mkl_core.lib
SLICOTLIB    = $(HOME)\slicot.lib
