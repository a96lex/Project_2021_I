# Project_2021_I

## Main files

- init.f90 : initialization subroutines. Responsible: Eloi Sanchez.
- parameters.f90 : common parameters for the program. Responsible: Eloi Sanchez.
- pbc.f90 : boundary conditions subroutines. Responsible: David March.
- integraforces.f90 : integration and dynamic (forces and energies) subroutines. Responsibles: Laia Barjuan (integration) and Arnau Jurado (forces).
- statvis.f90 : statistics and visualization subroutines. Responsibles: Jaume Ojer (statistics) and Alex Párrraga (visualization).

## Installation and usage

After cloning the repository with `git clone https://github.com/EIA-Master/Project_2021_I` you can generate an executable file with `make`, which will generate an executable called main.x.

You can run the program by executing `./main.x input.txt` with any proper input file (check the file input_template.txt to see the structure of the input file).
Alternatively `make sim input="input.txt"` with any file as input.txt.

### Additional Makefile recipes

- `make sim input="your_input.txt"`: Compiles the program if needed and runs a simulation with the parameters specified in `your_input.txt`
- `make plots`: (see necessary packages) Generates plots from the result data in the results/plots directory
- `make trajectory_video`: (see necessary packages) Generates an animated image (in gif format) with the resulting trajectory. The filename is defaulted to "trajectory.gif", but you can override it by passing an optional filename parameter: `make filename="your_filename" trajectory_video`
- `make backup`: Generates a directory named `results_Y-M-D_H:M:S` and copies the contents of the results folder to it. Use this when you want to save the results obtained febore running the executable again.
- `make clean`: Removes all intermediate files created during the build
- `make clean_all`: Removes all intermediate files created during the build and removes the results directory and all its content

### Necessary packages

- Gnuplot: it is a requirement to use `make plots`. Available [here](http://www.gnuplot.info/)
- Python version 3: is a requirement to use `make trajectory_video`
  - The only required python package is ase ([documentation](https://wiki.fysik.dtu.dk/ase/)). It can be installed using pip:
    - `pip install ase`
or
    - `pip install molecule_plotter/requirements.txt` from the root directory

## Contributors

| Laia Barjuan                                                                   | Arnau Jurado                                                             | David March                                                            | Jaume Ojer                                                                 | Alex Párraga                                                         | Eloi Sanchez                                                                   |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------ | ---------------------------------------------------------------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| ![laiabarjuan](https://avatars.githubusercontent.com/u/79266111 "laiabarjuan") | ![arnau-jr](https://avatars.githubusercontent.com/u/48213666 "arnau-jr") | ![dmarchp](https://avatars.githubusercontent.com/u/79266176 "dmarchp") | ![jaumeojer](https://avatars.githubusercontent.com/u/79266127 "jaumeojer") | ![a96lex](https://avatars.githubusercontent.com/u/62766970 "a96lex") | ![EloiSanchez](https://avatars.githubusercontent.com/u/79266117 "EloiSanchez") |
| [laiabarjuan](https://github.com/laiabarjuan)                                  | [arnau-jr](https://github.com/arnau-jr)                                  | [dmarchp](https://github.com/dmarchp)                                  | [jaumeojer](https://github.com/jaumeojer)                                  | [a96lex](https://github.com/a96lex)                                  | [EloiSanchez](https://github.com/EloiSanchez)                                  |

### Contact

Voice any concerns to the responsible of this github page: arnau-jr.
