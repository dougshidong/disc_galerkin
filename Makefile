MF=     Makefile

FC=     gfortran
FFLAGS= -O3
FFLAGS= -g -Wall -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan


LFLAGS= $(FFLAGS)
LIBS= -L/usr/lib -llapack -lblas

EXE=    advection

SRC= \
		prec.f90    \
		glob.f90    \
		quadrature.f90	\
		jacobi.f90	\
		poly.f90	\
		matrices.f90	\
		grid.f90	\
		legendreGLNodesWeights.f90	\
		legendreNM.f90	\
		lift1D.f90	\
		advecRHS1D.f90	\
		main.f90

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=    $(SRC:.f90=.o)

.f90.o:
		$(FC) $(FFLAGS) -c $< $(LIBS)

all:    $(EXE)

$(EXE): $(OBJ)
		$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): $(MF)

tar:
		tar cvf $(EXE).tar $(MF) $(SRC)

clean:
		rm -f $(OBJ) $(EXE) *.mod core
