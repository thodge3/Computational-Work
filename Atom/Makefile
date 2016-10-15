

.PHONY : all
all: atomproject.f90
	gfortran -c ForceEnergyMod.f90
	gfortran -c SolverMod.f90
	gfortran -c StatisticsMod.f90
	gfortran -o atomproject.x atomproject.f90 ForceEnergyMod.o SolverMod.o StatisticsMod.o
	./atomproject.x
	




