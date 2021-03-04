module init
    implicit none
    contains

        subroutine get_param(unit)
            !Author: Eloi Sanchez
            ! Llegim de l'input els parametres del sistema i es calcula el
            ! nº d'iteracions i la longitud de la cel·la.
            ! Els propers llocs on es faci servir 'use parameters' tindran
            ! les variables acualitzades
            
            use parameters
            implicit none

            integer, intent(in):: unit
            integer :: errstat

            namelist /input/ N, D, rho, dt_sim, n_meas, n_conf, n_equil, T_ref, fact_rc, sigma, epsilon, mass

            ! Llegim els parametres del input
            read(unit=unit, nml=input, iostat=errstat)
            
            if (errstat > 0) then
                print *, "ERROR reading namelist from input file in init.f90 (code", errstat, ")"
                stop
            end if
            
            ! Calculem el nº de iteracions i la longitud de la cel·la
            n_total = n_meas * n_conf
            L = (N / rho) ** (1.d0 / D)
            rc = fact_rc * L / 2.d0

            call reduced_units()
        end subroutine get_param

        subroutine init_sc(pos)
            !Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions. Es menys eficient que tenir nested loops pero
            ! no se m'ha acudit fer-ho d'altra forma.
            ! Funciona per N^(1/D) no exactes, però la densitat NO sera la requerida

            use parameters, only : D, N, L
            implicit none
            
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a, r  ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j, aux
            real*8 :: raux
            
            raux = N ** (1.d0 / D)
            M = ceiling(raux)
            if (abs(raux - M) > 1.d-3) print *, "WARNING: The number of atoms per dimension is approximated", raux, "->", M
            a = L / dble(M)

            r = - a  ! Necessari per que comenci pel 0 en el i=1 j=1. Si algu té alguna suggerencia o millora que la faci.
            do i = 1, D
                aux = M ** (D - i)
                do j = 1, N
                    if (mod(j - 1, aux) == 0) r = r + a
                    if (r > L - a) r = 0.d0
                    pos(i,j) = r
                end do
            end do
            pos = pos - (L - a)/2.d0 ! Centrem el sistema al (0,0,0)

        end subroutine

        subroutine init_vel(vel, T)
            !Author: Eloi Sanchez
            ! Torna el array de velocitats (aleatories) vel(D,N) consistent amb la T donada.

            use parameters, only : D, N
            use integraforces, only : energy_kin
            implicit none

            real*8, intent(inout) :: vel(D,N)
            real*8, intent(in) :: T
          
            real*8 :: vel_CM(D)
            real*8 :: aux, kin
            integer :: i, j
          
            ! Inicialitza les velocitats de manera random entre -1 i 1
            vel_CM = 0
            do i = 1, N
                do j = 1, D
                    vel(j,i) = 2.*rand() - 1.
                end do
              vel_CM = vel_CM + vel(:,i)
            end do
            vel_CM = vel_CM/dble(N)
          
            ! Eliminem la velocitat neta del sistema
            do i = 1, N
              vel(:,i) = vel(:,i) - vel_CM
            end do
            
            ! Reescalem les velocitats a la temperatura objectiu
            call energy_kin(vel, kin, aux)  
            vel = vel * sqrt(dble(3*N)*T/(2.d0*kin))

        end subroutine

end module init
