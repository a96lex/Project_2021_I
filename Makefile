objects= main.o parameters.o init.o pbc.o statvis.o
dep_objects= parameters.o init.o pbc.o statvis.o
mods= parameters.mod init.mod pbc.mod statvis.mod
compiler=gfortran
opt=

main.x : $(objects)
	$(compiler) -o main.x $(opt) $(objects)

$(mods) : $(dep_objects)  

parameters.o : parameters.f90
	$(compiler) -c $(opt) parameters.f90

init.o : init.f90 parameters.o
	$(compiler) -c $(opt) init.f90 parameters.f90

pbc.o : pbc.f90 parameters.o
	$(compiler) -c $(opt) pbc.f90 parameters.f90

statvis.o : statvis.f90 parameters.o
	$(compiler) -c $(opt) statvis.f90 parameters.f90

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90 parameters.f90 init.f90 pbc.f90


.PHONY: clean
clean :
	rm -f $(objects)
	rm -f $(mods)

