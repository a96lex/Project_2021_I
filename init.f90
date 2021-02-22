module init
    implicit none
    contains

        subroutine get_param(input_name)
            !
            ! Llegim de l'input els parametres del sistema i es calcula el
            ! nº d'iteracions i la longitud de la cel·la.
            ! Els propers llocs on es faci servir 'use parameters' tindran
            ! les variables acualitzades
            ! Eloi (no se ni si compila 22/02 pq falta el parameters.mod)
            !
            use parameters
            implicit none

            character (len=50), intent(in):: input_name
            integer :: errstat

            namelist /input/ N, D, rho, dt, n_meas, n_conf

            ! Open and read namelist from input file
            open(10, file=input_name, status="old")
            read(unit=10, nml=input, iostat=errstat)
            close(10)
            
            if (errstat > 0) then
                print *, "ERROR reading namelist from input file in init.f90 (code", errstat, ")"
                stop
            end if
            
            ! Calculem el nº de iteracions i la longitud de la cel·la
            n_total = n_meas * n_conf
            L = (N / rho) ** (1.d0 / D)

        end subroutine get_param

end module init
