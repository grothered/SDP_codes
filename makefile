#Another attempt at a makefile. Note that asa066 and asa241 contain libraries for the standard normal cumulative distribution function and its inverse. 
model: SDP_take2.o SDP_driver.f90
	gfortran -O3 -o SDP_run2 SDP_take2.o SDP_driver.f90

#SDP_take2.mod: SDP_take2.o SDP_take2.f90
#	gfortran -c -O3 SDP_take2.f90

parallel_model: SDP_take2.o SDP_driver.f90
	gfortran -O3 -fopenmp -o SDP_run2 SDP_take2.o SDP_driver.f90


SDP_take2.o: SDP_take2.f90
	gfortran -c -O3 SDP_take2.f90
