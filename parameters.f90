      module parameters
      implicit none

      integer :: N,D
      real*8 :: dt_sim
      integer :: n_meas, n_conf, n_total
      real*8 :: rho, T_ref, L, rc, fact_rc
      real*8 :: sigma, epsilon, mass

      real*8 :: unit_of_time,unit_of_energy,unit_of_length,unit_of_pressure
      real*8,parameter :: boltzmann_k = 8.31446261815324 !J/mol, technically R.

      contains

            subroutine reduced_units()
            !Author: Arnau Jurado
            ! Computes the units of time and energy from the LJ parameters
                  implicit none
                  unit_of_time = sigma*sqrt((mass/1000)/epsilon)/sqrt(boltzmann_k)/100.d0!ps
                  unit_of_energy = epsilon*boltzmann_k !J/mol
                  unit_of_length = sigma !Angstroms
                  unit_of_pressure = unit_of_energy/(unit_of_length)**3*(10.d0**(-10))**3 !Pa
            end subroutine reduced_units
      end module parameters
