# Project_2021_I

## Project description

The code of this github page performs a molecular dynamics simulation of an atomic gas/liquid from the parameters defined on an input file (see `input_template.txt` for a template) and outputs different observables from the simulation. The results include:

- Time series and averages of kinetic, potential and total energy of the system.
- Time series and averages of the instanteous temperature and pressure.
- Radial distribution function.
- Correlation function of the total energy of the system.
- Diffusion coefficients for each coordinate.

Scripts for generating plots of the results are also included.

The program has a serial version (in `serie` directory) and a parallel version (in `parallel` directory). Most make functions are the same for both versions.

## Main files
The main file structure is the same for both serie and parallel programs.
- init.f90 : initialization subroutines. Responsible: Eloi Sanchez.
- parameters.f90 : common parameters for the program. Responsible: Eloi Sanchez.
- pbc.f90 : boundary conditions subroutines. Responsible: David March.
- rad_dist.f90 : radial distribution function subroutines. Responsible: David March.
- integraforces.f90 : integration and dynamic (forces and energies) subroutines. Responsibles: Laia Barjuan (integration) and Arnau Jurado (forces).
- statvis.f90 : statistics and visualization subroutines. Responsibles: Jaume Ojer (statistics) and Alex Párrraga (visualization).
- plots.p : `gnuplot` script for generating all plots. Responsible: Alex Párraga.

## Installation and usage

After cloning the repository with `git clone https://github.com/EIA-Master/Project_2021_I` you can generate an executable file with `make` while inside either of the directories `serie` or `parallel`. This will generate an executable called main.x.

You can run the program by executing `./main.x input.txt` (in serie) or `mpirun main.x input.txt` (in parallel) with any proper input file (check the file input_template.txt to see the structure of the input file).
Alternatively `make sim input="input.txt"` with any file as input.txt.

After performing a simulation the results will be present in the `results` directory. You can generate plots from the simulation data by running the gnuplot script (`gnuplot plots.p`, requires `/results/plots/dimensionalized` to be created manually) or by the Makefile rule `plots` (all directory paths are automatically created).

Warning: if you run another simulation the results from the previous one will be erased. You can back them up manually or by using the `backup` Makefile rule.

### Additional Makefile recipes

- `make sim input="your_input.txt"`: Compiles the program if needed and runs a simulation with the parameters specified in `your_input.txt`. In the parallel version it will use all available processors.
- `make simN numproc=N input="your_input.txt"` (only in parallel): Compiles the program if needed and runs a simulation with the parameters specified in `your_input.txt` with `N` tasks.
- `make plots`: (see necessary packages) Generates plots from the result data in the results/plots directory.
- `make trajectory_video`: (see necessary packages) Generates an animated image (in gif format) with the resulting trajectory. The filename is defaulted to "trajectory.gif", but you can override it by passing an optional filename parameter: `make filename="your_filename" trajectory_video`.
- `make backup`: Generates a directory named `results_Y-M-D_H:M:S` and copies the contents of the results folder to it. Use this when you want to save the results obtained febore running the executable again.
- `make clean`: Removes all intermediate files created during the build.
- `make clean_all`: Removes all intermediate files created during the build and removes the results directory and all its content.
- `make strong_speed_up`: (diagnostics) Executes a bash script to evaluate the strong speedup (same load, increasing core number) of the program. Needs acces to a Sun N1 Grid Engine cluster setup, not intended for general usage.  
### Necessary packages
- OpenMPI Fortran 90 compilers (tested with gcc compiler version 4.8.5 and 9.3.0) and OpenMPI release 2.1.1 or later. Only required for parallel program. Avaliable [here](https://www.open-mpi.org/software/ompi/v4.1/) or through some linux package managers.
- Gnuplot: it is a requirement to use `make plots`. Available [here](http://www.gnuplot.info/) or through some linux package managers.
- Python version 3: is a requirement to use `make trajectory_video`.
  - The only required python package is ase ([documentation](https://wiki.fysik.dtu.dk/ase/)). It can be installed using pip: - `pip install ase`
    or - `pip install -r molecule_plotter/requirements.txt` from the root directory.

## Contributors

| Laia Barjuan                                                                   | Arnau Jurado                                                             | David March                                                            | Jaume Ojer                                                                 | Alex Párraga                                                         | Eloi Sanchez                                                                   |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------ | ---------------------------------------------------------------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| ![laiabarjuan](https://avatars.githubusercontent.com/u/79266111 "laiabarjuan") | ![arnau-jr](https://avatars.githubusercontent.com/u/48213666 "arnau-jr") | ![dmarchp](https://avatars.githubusercontent.com/u/79266176 "dmarchp") | ![jaumeojer](https://avatars.githubusercontent.com/u/79266127 "jaumeojer") | ![a96lex](https://avatars.githubusercontent.com/u/62766970 "a96lex") | ![EloiSanchez](https://avatars.githubusercontent.com/u/79266117 "EloiSanchez") |
| [laiabarjuan](https://github.com/laiabarjuan)                                  | [arnau-jr](https://github.com/arnau-jr)                                  | [dmarchp](https://github.com/dmarchp)                                  | [jaumeojer](https://github.com/jaumeojer)                                  | [a96lex](https://github.com/a96lex)                                  | [EloiSanchez](https://github.com/EloiSanchez)                                  |

### Contact

Voice any concerns to the responsible of this github page: [arnau-jr](https://github.com/arnau-jr).
