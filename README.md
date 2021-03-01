# Project_2021_I

##  Team

- Laia Barjuan.
- Arnau Jurado.
- David March.
- Jaume Ojer.
- Alex Párraga.
- Eloi Sanchez.

## Main files

- init.f90 : initialization subroutines. Responsible: Eloi Sanchez.
- parameters.f90 : common parameters for the program. Responsible: Eloi Sanchez.
- pbc.f90 : boundary conditions subroutines. Responsible: David March.
- integraforces.f90 : integration and dynamic (forces and energies) subroutines. Responsibles: Laia Barjuan (integration) and Arnau Jurado (forces).
- statvis.f90 : statistics and visualization subroutines. Responsibles: Jaume Ojer (statistics) and Alex Párrraga (visualization).

## Installation and usage
After cloning the repository with `git clone https://github.com/EIA-Master/Project_2021_I` you can generate an executable file with `make`, which will generate a an executable called main.x.

### Additional Makefile recipes
 - `make clean`: Removes all intermediate files created during the build
 - `make plots`: Generates plots from the result data in the results/plots directory
 - `make trajectoryVideo`: Generates an animated image (in gif format) with the resulting trajectory. The filename is defaulted to "trajectory.gif", but you can override it by passing an optional filename parameter: `make filename="your_filename" trajectoryVideo` 

### Contact
Voice any concerns to the responsible of this github page: arnau-jr. 
