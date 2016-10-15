
.Phony. : all
all : MainProgram.f90
	gfortran -c GaussianElimination.f90
	gfortran -c BuildMatrix.f90
	gfortran -o MainProgram.x MainProgram.f90 BuildMatrix.f90 GaussianElimination.f90 -llapack
	./MainProgram.x 



