objects= main.o parameters.o init.o integraforces.o statvis.o
dep_objects= parameters.o init.o integraforces.o statvis.o
mods= parameters.mod init.mod integraforces.mod statvis.mod
compiler=mpif90
# opt=-Wall
opt=-g -Wall

main.x : $(objects)
	mkdir -p results
	$(compiler) -o main.x $(opt) $(objects)

$(mods) : $(dep_objects)  

parameters.o : parameters.f90
	$(compiler) -c $(opt) parameters.f90

init.o : init.f90 parameters.o integraforces.o
	$(compiler) -c $(opt) init.f90 

integraforces.o : integraforces.f90 parameters.o
	$(compiler) -c $(opt) integraforces.f90 

statvis.o : statvis.f90
	$(compiler) -c $(opt) statvis.f90 

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90


.PHONY: sim plots trajectory_video clean backup clean_all
## sim: Compiles the program if needed and runs a simulation with the parameters specified with input=your_input.txt.
## Results are generated at the directory "results"
sim : main.x
	./main.x $(input)	

## plots: Generates the plots of the observables found in results/
plots :
	gnuplot plots.p

## trajectory_video: Creates the movie from results/trajectory.xyz
trajectory_video:
	python3 molecule_plotter/main.py $(filename)

## clean: Removes the compilation files
clean :
	rm -f $(objects) $(mods) 

## clean_all: Removes the compilation files and the results/ folder
clean_all :
	rm -rf -v results
	rm -f $(objects) $(mods)

## backup: Saves the contents from the results folder in the backup folder with a timestamp
backup : 
	mkdir -p backups results
	cp -a results "backups/results_$(shell date +"%Y-%m-%d_%H:%M:%S")"

## help: Shows information about the Makefile options
help :
	@grep '^##' Makefile
