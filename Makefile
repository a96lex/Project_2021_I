objects= main.o parameters.o init.o pbc.o integraforces.o statvis.o rad_dist.o
dep_objects= parameters.o init.o pbc.o integraforces.o statvis.o rad_dist.o
mods= parameters.mod init.mod pbc.mod integraforces.mod statvis.mod radial_distribution.mod
compiler=gfortran
# opt=-Wall
opt=

main.x : $(objects)
	$(compiler) -o main.x $(opt) $(objects)

$(mods) : $(dep_objects)  

parameters.o : parameters.f90
	$(compiler) -c $(opt) parameters.f90

init.o : init.f90 parameters.o integraforces.o
	$(compiler) -c $(opt) init.f90 

pbc.o : pbc.f90 parameters.o
	$(compiler) -c $(opt) pbc.f90

rad_dist.o : rad_dist.f90 parameters.o pbc.o
	$(compiler) -c $(opt) rad_dist.f90

integraforces.o : integraforces.f90 parameters.o pbc.o
	$(compiler) -c $(opt) integraforces.f90 

statvis.o : statvis.f90 parameters.o
	$(compiler) -c $(opt) statvis.f90 

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90


.PHONY: clean
clean :
	rm -f $(objects)
	rm -f $(mods)

