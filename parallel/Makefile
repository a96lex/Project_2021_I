objects= main.o parameters.o init.o integraforces.o statvis.o pbc.o rad_dist.o
dep_objects= parameters.o init.o integraforces.o statvis.o pbc.o rad_dist.o
mods= parameters.mod init.mod integraforces.mod statvis.mod pbc.mod rad_dist.mod
compiler=mpif90
opt=-g

main.x : $(objects)
	mkdir -p results
	mkdir -p results/dimensionalized
	$(compiler) -o main.x $(opt) $(objects)

$(mods) : $(dep_objects)  

parameters.o : parameters.f90
	$(compiler) -c $(opt) parameters.f90

init.o : init.f90 parameters.o integraforces.o
	$(compiler) -c $(opt) init.f90 

integraforces.o : integraforces.f90 parameters.o pbc.o rad_dist.o
	$(compiler) -c $(opt) integraforces.f90 

statvis.o : statvis.f90
	$(compiler) -c $(opt) statvis.f90
	
pbc.o : pbc.f90 parameters.o
	$(compiler) -c $(opt) pbc.f90
	
rad_dist.o : rad_dist.f90 parameters.o
	$(compiler) -c $(opt) rad_dist.f90

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90


.PHONY: sim simN plots trajectory_video clean backup clean_all strong_speed_up
## sim: Compiles the program if needed and runs a simulation with the parameters specified with input=your_input.txt. Uses all processors/workers available.
## Results are generated at the directory "results"
sim : main.x
	mpirun main.x $(input)	

## simN: Compiles the program if needed and runs a simulation with the parameters specified with input=your_input.txt. Uses numproc=N processors/workers.
## Results are generated at the directory "results"
simN : main.x
	mpirun -np $(numproc) main.x $(input) > results/speedup/sim_$(numproc)/sim$(numproc).log	

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
	cp input_template.txt results/.
	cp -a results "backups/results_$(shell date +"%Y-%m-%d_%H:%M:%S")"

## strong_speed_up: Executes the script for different processor number
strong_speed_up :
	bash scripts/strong_speed_up_2.sh

## help: Shows information about the Makefile options
help :
	@grep '^##' Makefile
