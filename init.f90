module init
    implicit none
    contains

        subroutine get_param(unit)
            !
            ! Llegim de l'input els parametres del sistema i es calcula el
            ! nº d'iteracions i la longitud de la cel·la.
            ! Els propers llocs on es faci servir 'use parameters' tindran
            ! les variables acualitzades
            ! Eloi (no se ni si compila 22/02 pq falta el parameters.mod)
            !
            use parameters
            implicit none

            integer, intent(in):: unit !AJ: l'asterisk el fa mes generic
            integer :: errstat

            namelist /input/ N, D, rho, dt, n_meas, n_conf, T_ref, fact_rc, sigma, epsilon

            ! Open and read namelist from input file
            read(unit=unit, nml=input, iostat=errstat)
            
            if (errstat > 0) then
                print *, "ERROR reading namelist from input file in init.f90 (code", errstat, ")"
                stop
            end if
            
            ! Calculem el nº de iteracions i la longitud de la cel·la
            n_total = n_meas * n_conf
            L = (N / rho) ** (1.d0 / D)
            rc = fact_rc * L / 2.d0

        end subroutine get_param

end module init
