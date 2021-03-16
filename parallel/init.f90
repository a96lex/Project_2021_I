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
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Funciona per N^(1/D) no exactes

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
            
        end subroutine init_sc

        subroutine init_sc_paralel(pos)
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Funciona per N^(1/D) no exactes

            use parameters, only : D, N, L, taskid, numproc, master
            use mpi
            implicit none
            
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a, r  ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j, index_reset
            real*8 :: M_check, check

            integer, allocatable :: N_block(:), N_block_local(:)
            integer :: task_dim, ierror

            allocate(N_block(D))
            allocate(N_block_local(D))

            ! Aixo ho fan tots els cores
            M_check = N ** (1.d0 / D)
            M = ceiling(M_check)
            a = L / dble(M)
            
            ! IDEA 1. Crec que pitjor
            ! Trobem el nº de cores per dimensio
            do i = 1, D
                N_block(i) = numproc / D
                if (i <= mod(numproc, D)) N_block(i) = N_block(i) + 1
                if (taskid == master) print *, N_block(i)
            end do

            ! ara tocaria assignar a cada core la dim que li toqui

            ! IDEA 2. Crec que millor. Problema: Que passa si numproc < D ????
            ! Necessari per que si no se li enva la pinça al reduce
            N_block = 0
            N_block_local = 0

            ! Assignem a cada core una dimensio
            task_dim = D * taskid / numproc + 1
            N_block_local(task_dim) = 1
            print *, "taskid = ", taskid, "task_dim = ", task_dim

            ! Mirem quants nodes tenim per cada dimensio
            call MPI_Allreduce(N_block_local, N_block, D, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
            if (taskid == master) print *, N_block

            ! r = - a  ! Necessari per que comenci pel 0 en el i=1 j=1.
            ! do i = 1, D
            !     index_reset = M ** (D - i)
            !     do j = 1, N
            !         if (mod(j - 1, index_reset) == 0) r = r + a
            !         ! Si la posicio correspon al final de la caixa, reiniciem r
            !         if (abs(r - L) < 1.d-6) r = 0.d0
            !         pos(i,j) = r
            !     end do
            ! end do
            ! pos = pos - (L - a) / 2.d0 ! Centrem el sistema al (0,0,0)
            
        end subroutine

        ! subroutine init_vel(vel, T)
        !     !Author: Eloi Sanchez
        !     ! Torna el array de velocitats (aleatories) vel(D,N) consistent amb la T donada.

        !     use parameters, only : D, N
        !     use integraforces, only : energy_kin
        !     implicit none

        !     real*8, intent(inout) :: vel(D,N)
        !     real*8, intent(in) :: T
          
        !     real*8 :: vel_CM(D)
        !     real*8 :: aux, kin
        !     integer :: i, j
          
        !     ! Inicialitza les velocitats de manera random entre -1 i 1
        !     vel_CM = 0
        !     do i = 1, N
        !         do j = 1, D
        !             vel(j,i) = 2.*rand() - 1.
        !         end do
        !       vel_CM = vel_CM + vel(:,i)
        !     end do
        !     vel_CM = vel_CM/dble(N)
          
        !     ! Eliminem la velocitat neta del sistema
        !     do i = 1, N
        !       vel(:,i) = vel(:,i) - vel_CM
        !     end do
            
        !     ! Reescalem les velocitats a la temperatura objectiu
        !     call energy_kin(vel, kin, aux)  
        !     vel = vel * sqrt(dble(3*N)*T/(2.d0*kin))

        ! end subroutine

end module init
