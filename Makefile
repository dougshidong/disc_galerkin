MF=     Makefile

FC=     gfortran
FV=     -std=f2008
#FFLAGS= -O3
#FFLAGS= -g -Wall -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
#FFLAGS= -g -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
#FFLAGS= -pg -fprofile-arcs -ftest-coverage -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=na
#FFLAGS= -pg -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
FFLAGS= -O3 -fimplicit-none -ffree-line-length-0 -fdefault-real-8



LFLAGS= $(FFLAGS)
LIBS= -L/usr/lib -llapack -lblas

EXE=    advection

SRC= \
		gmsh_io.f90	\
		prec.f90    \
		glob.f90    \
    	legendreGLNodesWeights.f90	\
    	legendreNM.f90	\
    	quadrature.f90	\
		cubature2D.f90	\
		ref_ele_mod.f90\
    	face_mod.f90\
    	bspline.f90	\
    	nurbs.f90	\
    	jacobi.f90	\
    	monomial.f90	\
        bezier.f90      \
	lagrange.f90	\
    	poly.f90	\
    	matrices.f90	\
    	element_mod.f90\
    	grid.f90	\
    	advecRHS2D.f90	\
		testcases.f90 	\
		main.f90
### 	advecRHS2D_old.f90 	\
####	lift1D.f90	\
####	stability.f90	\

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=    $(SRC:.f90=.o)

.f90.o:
		$(FC) $(FV) $(FFLAGS) -c $< $(LIBS)

all:    clean $(EXE)

$(EXE): $(OBJ)
		$(FC) $(FV) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): $(MF)

tar:
		tar cvf $(EXE).tar $(MF) $(SRC)

clean:
		rm -f $(OBJ) $(EXE) *.mod core
