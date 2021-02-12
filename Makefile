# This makefile is used to compile new HAMISH code.
# The compiler: gfortran compiler
#
#FCFLAGS=  -Wuse-without-only -g
#FC=mpif90
FC=h5pfc

SRCDIR = src
OBJDIR = obj
BINDIR = bin

FCFLAGS= -O3

# OPTIONS1 = -fcheck=all
OPTIONS2 = -J $(OBJDIR)
OPTIONS3 = -DHDF5 -Dcputime
# OMP = -fopenmp
# OPTIONS1 = -D__GFORTRAN__ -ffree-line-length-none
#OPTIONS2 = -fno-range-check -fcheck=all
#OPTIONS3 = -DEULER
#OPTIONS3=-O3

EXE=astrr.exe

LIBS= -lz -lm

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

srs=  constdef.F90 tecio.F90 commvar.F90 commarray.F90 fludyna.F90 parallel.F90 \
      hdf5io.F90 readwrite.F90 commfunc.F90 gridgeneration.F90 solver.F90       \
      initialisation.F90 statistic.F90 mainloop.F90 astrr.F90
OBJS=$(srs:.F90=.o)

%.o:%.F90
	@mkdir -p $(OBJDIR) 
	$(FC) $(FCFLAGS) $(OPTIONS1) $(OPTIONS2) $(OPTIONS3) $(OMP) -c -o $(OBJDIR)/$@  $<

default: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS) $(OMP)

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET)
