####################################################################
#  SLICOT main makefile                                            #
#  Top Level Makefile for generating SLICOT Library object file,   #
#  the auxiliary library file, and linking and running the example #
#  programs.                                                       #
#  SLICOT, Release 5.0                          .\slicot\makefile  #
#  October 31, 1996                                                #
#  Revised December 7, 1999; January 2005.                         #
####################################################################
#
#  This makefile creates/updates the SLICOT Library object file, the
#  auxiliary library, and compiles, links, and runs the example 
#  programs for the SLICOT Library on Windows platforms.
#  To perform all these actions, enter
#       nmake
#
#  To create/update the library for SLICOT and auxiliary library, 
#  enter 
#       nmake lib
#
#  To compile, link, and run the example programs, enter 
#       nmake example
#
#  To remove the object files for SLICOT routines and auxiliary 
#  routines, enter
#       nmake cleanlib
#
#  To remove the files with the computed results (*.exa), enter
#       nmake cleanexample
#
#  To remove the object files for SLICOT routines and auxiliary 
#  routines, as well as the files with the computed results (*.exa),
#  enter
#       nmake clean
#
####################################################################

!include <make.inc>

all: lib example

clean: cleanlib cleanexample

lib:
	cd $(HOME)\src
	$(MAKE)
	cd $(HOME)\src_aux
	$(MAKE)

example:
	cd $(HOME)\examples
	$(MAKE)

cleanlib:
	cd $(HOME)\src
	$(MAKE) clean
	cd $(HOME)\src_aux
	$(MAKE) clean

cleanexample:
	cd $(HOME)\examples 
	$(MAKE) clean
