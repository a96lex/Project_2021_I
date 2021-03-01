# Project_2021_I

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

## Contributors

| Laia Barjuan                                                                   | Arnau Jurado                                                             | David March                                                            | Jaume Ojer                                                                 | Alex Párraga                                                                                                              | Eloi Sanchez                                                                   |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------ | ---------------------------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| ![laiabarjuan](https://avatars.githubusercontent.com/u/79266111 "laiabarjuan") | ![arnau-jr](https://avatars.githubusercontent.com/u/48213666 "arnau-jr") | ![dmarchp](https://avatars.githubusercontent.com/u/79266176 "dmarchp") | ![jaumeojer](https://avatars.githubusercontent.com/u/79266127 "jaumeojer") | ![a96lex](https://avatars.githubusercontent.com/u/62766970?s=460&u=d10e5a4a565ae28797bf4a1bf118b73cf371b720&v=4 "a96lex") | ![EloiSanchez](https://avatars.githubusercontent.com/u/79266117 "EloiSanchez") |
| [laiabarjuan](https://github.com/laiabarjuan)                                  | [arnau-jr](https://github.com/arnau-jr)                                  | [dmarchp](https://github.com/dmarchp)                                  | [jaumeojer](https://github.com/jaumeojer)                                  | [a96lex](https://github.com/a96lex)                                                                                       | [EloiSanchez](https://github.com/EloiSanchez)                                  |

### Contact

Voice any concerns to the responsible of this github page: arnau-jr.
