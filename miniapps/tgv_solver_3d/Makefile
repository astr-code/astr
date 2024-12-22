# This makefile is used to compile ASTR code.
# The compiler: gfortran compiler
#
#FCFLAGS=  -Wuse-without-only -g
#FC=mpif90
FC=h5pfc

SRCDIR = ./oacc/ ./com/

FCFLAGS= -O3 -acc -ffast-math

EXE=astr_oacc

TARGET = $(BINDIR)/$(EXE)

VPATH = $(SRCDIR):$(OBJDIR)

srs=  comvardef.f90 constdef.f90 utility.f90 io.f90 tecio.f90 xdmwrite.f90 \
      fluidynamcs.f90 bc.f90 field_init.f90 numerics.f90 statistics.f90   \
      solver.f90 /tgvsolver_oacc.f90 
      
OBJS=$(srs:.F90=.o)

%.o:%.F90
	@mkdir -p $(OBJDIR) 
	$(FC) $(FCFLAGS) $(INCL) $(OPTIONS1) $(OPTIONS2) $(OPTIONS3) $(OPTIONS4) $(OMP) -c -o $(OBJDIR)/$@  $<

default: $(OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJDIR)/*.o $(LIBS) $(INCL) $(OMP)

clean:
	rm -fv $(OBJDIR)/*.o $(OBJDIR)/*.mod $(TARGET) $(OBJDIR)/*.mod
