
.Phony. : all
all : MainProgram.f90
	gfortran -c TimeMod.f90
	gfortran -c SolverMod.f90
	gfortran -o MainProgram.x MainProgram.f90 TimeMod.f90 SolverMod.f90
	./MainProgram.x



