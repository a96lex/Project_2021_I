objects= main.o parameters.o init.o pbc.o integraforces.o statvis.o
dep_objects= parameters.o init.o pbc.o integraforces.o statvis.o
mods= parameters.mod init.mod pbc.mod integraforces.mod statvis.mod
source_files= init pbc integraforces statvis
compiler=gfortran
opt=

main.x : $(objects)
	$(compiler) -o main.x $(opt) $(objects)

$(mods) : $(dep_objects)  

parameters.o : parameters.f90
	$(compiler) -c $(opt) parameters.f90

$(addsuffix .o,$(source_files)) : $(addsuffix .f90,$(source_files))  parameters.o
	$(compiler) -c $(opt) $(addsuffix .f90,$(source_files))  parameters.f90

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90 parameters.f90 init.f90 integraforces.f90 pbc.f90


.PHONY: clean
clean :
	rm -f $(objects)
	rm -f $(mods)

