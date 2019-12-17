#!/bin/sh

#------------------------  edit this block ---------------------------------

# the name of the executable program to be created
PROG_NAME = tst

# hpchem root
hpc.dir := ../../src

# Directory of the current application
app.dir := .

# objects of the current application
app.obj := p.o

hpc.obj := thermo_idealGas.o thermo_PengRobinson.o hpchem.o

# additional flags to be passed to the linker. If your program
LINK_OPTIONS =

#-----------------------------------------------------------------------------
# paths
VPATH := $(app.dir) $(hpc.dir)

# objects
OBJS := $(hpc.obj) $(app.obj)

#---------------------------------------------------------------------------

# the Fortran 90/95 compiler
F90 = ifort

# the Fortran 77 compiler
F77 = ifort

# Fortran compile flags
FORT_FLAGS =  -O3

# Fortran libraries used to link fortran main programs
FLIBS=

# external libraries

%.o : %.f90
	$(F90) -c $< $(FORT_FLAGS)

%.o : %.f
	$(F77) -c $< $(FORT_FLAGS)


PROGRAM = $(PROG_NAME)$(EXE_EXT)

DEPENDS = $(OBJS:.o=.d)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F90) -o $(PROGRAM) $(OBJS) $(LCXXFLAGS) $(LINK_OPTIONS) \
                    $(EXT_LIBS) $(FLIBS)


clean:
	$(RM) $(OBJS) $(PROGRAM) *.mod

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS)

ifeq ($(wildcard .depends), .depends)
include .depends
endif