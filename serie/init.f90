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

            namelist /input/ N, D, rho, dt_sim, n_meas, n_conf, n_equil, T_ref, & 
                              fact_rc, sigma, epsilon, mass, seed

            ! Llegim els parametres del input
            read(unit=unit, nml=input, iostat=errstat)
            
            if (errstat > 0) then
                print *, "ERROR reading namelist from input file in init.f90 (code", errstat, ")"
                stop
            end if
            call srand(seed)

            ! Calculem el nº de iteracions i la longitud de la cel·la
            n_total = n_meas * n_conf
            L = (N / rho) ** (1.d0 / D)
            rc = fact_rc * L / 2.d0

            call reduced_units()
        end subroutine get_param

        subroutine init_sc(pos)
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Funciona per N^(1/D) no exactes, però la densitat NO sera la requerida.

            use parameters, only : D, N, L
            implicit none
            
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a, r  ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j, index_reset
            real*8 :: M_check, check
            
            M_check = N ** (1.d0 / D)
            M = ceiling(M_check)
            a = L / dble(M)

            r = - a  ! Necessari per que comenci pel 0 en el i=1 j=1.
            do i = 1, D
                index_reset = M ** (D - i)
                do j = 1, N
                    if (mod(j - 1, index_reset) == 0) r = r + a
                    ! Si la posicio correspon al final de la caixa, reiniciem r
                    if (abs(r - L) < 1.d-6) r = 0.d0
                    pos(i,j) = r
                end do
            end do
            pos = pos - (L - a) / 2.d0 ! Centrem el sistema al (0,0,0)
            
            ! Sanity check. Les coord de l'ultim atom han de ser totes iguals!
            check = 0.d0
            do i = 1, D - 1
                check = check + abs(pos(i,N) - pos(i+1,N))
            end do

            if (check > 10.d-6) print *, "WARNING: Check that the initial conditions are correct!"
            if (abs(M_check - M) > 1.d-3) print *, "The number of atoms per dimension is approximated", M_check, "->", M

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
