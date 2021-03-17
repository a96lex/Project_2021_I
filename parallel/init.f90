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

        subroutine init_sc_inner(pos)
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Paralelització en el inner loop (les N particules) de la assignacio

            use parameters, only : D, N, L, taskid, numproc, master
            use mpi
            implicit none
            
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a, r  ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j, index_reset, index_control
            real*8 :: M_check, check

            real*8, allocatable:: pos_local(:,:)
            integer :: j_0, j_f, ierror
            integer :: local_size, all_size(numproc), all_position(numproc)

            ! Aixo ho fan totes les tasks
            M_check = N ** (1.d0 / D)
            M = ceiling(M_check)
            a = L / dble(M)
            
            ! Paralelitzem nomes el LOOP INTERN
            j_0 = taskid * N / numproc + 1
            j_f = (taskid + 1) * N / numproc

            ! Creem les variables locals de cada task
            allocate(pos_local(D,j_0:j_f))
            local_size = j_f - j_0 + 1

            ! Cada task dona valor a la posició local
            pos_local = 0.d0
            do i = 1, D
                index_reset = M ** (D - i)
                index_control = M ** (D - i + 1)
                do j = j_0, j_f
                    pos_local(i,j) = a * (mod(j - 1, index_control) / index_reset)
                end do
            end do
            pos_local = pos_local - (L - a) / 2.d0  ! Centrem el sistema al (0,0,0)

            ! El master ha de saber quanta informació rebra de cada core
            call MPI_Gather(local_size, 1, MPI_INTEGER, all_size, 1, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierror)

            if (taskid == master) then
                all_position(1) = 0  ! Ha de començar al 0 per coses del OpenMPI
                do i = 1, numproc - 1
                    all_position(i+1) = all_position(i) + all_size(i)
                end do
            end if

            ! Es guarda al master la posicio total a partir de les locals
            do i = 1, D
                call MPI_Gatherv(pos_local(i,:), local_size, MPI_DOUBLE_PRECISION, pos(i,:), &
                                all_size, all_position, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do

        end subroutine init_sc_inner

        subroutine init_sc_outer(pos)
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Paralelització en el outer loop (les N particules) de la assignacio
            ! Fins a un ordre de magnitud millor que la versio inner

            use parameters, only : D, N, L, taskid, numproc, master
            use mpi
            implicit none
            
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a, r  ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j, index_reset, index_control
            real*8 :: M_check, check

            real*8, allocatable:: pos_local(:,:)
            integer :: j_0, j_f, ierror
            integer :: local_size, all_size(numproc), all_position(numproc)

            ! Aixo ho fan totes les tasks
            M_check = N ** (1.d0 / D)
            M = ceiling(M_check)
            a = L / dble(M)
            
            ! Paralelitzem nomes el LOOP INTERN
            j_0 = taskid * N / numproc + 1
            j_f = (taskid + 1) * N / numproc

            ! Creem les variables locals de cada task
            allocate(pos_local(D,j_0:j_f))
            local_size = j_f - j_0 + 1

            ! Cada task dona valor a la posició local
            pos_local = 0.d0
            do j = j_0, j_f
                do i = 1, D 
                    pos_local(i,j) = a * (mod(j - 1, M ** (D - i + 1)) / M ** (D - i))
                end do
            end do
            pos_local = pos_local - (L - a) / 2.d0  ! Centrem el sistema al (0,0,0)

            ! El master ha de saber quanta informació rebra de cada core
            call MPI_Gather(local_size, 1, MPI_INTEGER, all_size, 1, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierror)

            if (taskid == master) then
                all_position(1) = 0  ! Ha de començar al 0 per coses del OpenMPI
                do i = 1, numproc - 1
                    all_position(i+1) = all_position(i) + all_size(i)
                end do
            end if

            ! Es guarda al master la posicio total a partir de les locals
            do i = 1, D
                call MPI_Gatherv(pos_local(i,:), local_size, MPI_DOUBLE_PRECISION, pos(i,:), &
                                all_size, all_position, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do

        end subroutine init_sc_outer

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
