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
                            fact_rc, sigma, epsilon, mass

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
            call divide_particles()
        end subroutine get_param

        subroutine divide_particles()
        !Author: Arnau Jurado & Eloi Sanchez
        !Divides the work among the processors by assigning each one an "imin" and
        !a "imax", which are the indexes of the first and last particle they have
        !to process e.g. with forces, each processor computes the forces
        !from the imin-th particle to the imax-th particles, both included.
            use parameters
            implicit none
            include 'mpif.h'
            integer :: i, local_size
  
            integer :: ierror
            
            imin = taskid * N / numproc + 1
            imax = (taskid + 1) * N / numproc
            local_size = imax - imin + 1
            ! print*,taskid,imin,imax,imax-imin+1
            
            ! We create aux_size(numproc) and aux_pos(numproc) only in master
            call MPI_Gather(local_size, 1, MPI_INTEGER, aux_size, 1, MPI_INTEGER, &
            master, MPI_COMM_WORLD, ierror)
            
            if (taskid == master) then
                aux_pos(1) = 0  ! Ha de començar al 0 per coses del OpenMPI
                do i = 1, numproc - 1
                    aux_pos(i+1) = aux_pos(i) + aux_size(i)
                end do
            end if

        end subroutine divide_particles

        subroutine init_sc_gather(pos)
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Paralelització en el outer loop (les N particules) de la assignacio
            ! Fins a un ordre de magnitud millor que la versio del reduce

            use parameters, only : D, N, L, taskid, numproc, master
            implicit none
            include 'mpif.h'
            
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

            deallocate(pos_local)
        end subroutine init_sc_gather

        subroutine init_sc_reduce(pos)
            ! Author: Eloi Sanchez
            ! Crea una xarxa cristal·lina ordenada (cuadrada en 2D i cubica en 3D)
            ! Està fet general per D dimensions.
            ! P. ex. N=27 i L=3 tindrem atoms amb coords a 0, 1 i 2.
            ! En la dimensio 1 farem 000000000111111111222222222
            ! En la dimensio 2 farem 000111222000111222000111222
            ! En la dimensio 3 farem 012012012012012012012012012
            ! Així, cada columna indica les 3 coord de un atom. Al final es centra la grid.
            ! Paralelització en el outer loop (les N particules) i reduce enlloc de gather
            ! Del ordre del gatherv pero una mica mes lenta.

            use parameters, only : D, N, L, taskid, numproc, master
            implicit none
            include 'mpif.h'
            
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a, r  ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j, index_reset, index_control
            real*8 :: M_check, check

            real*8 :: pos_local(D,N)
            integer :: j_0, j_f, ierror

            ! Aixo ho fan totes les tasks
            M_check = N ** (1.d0 / D)
            M = ceiling(M_check)
            a = L / dble(M)
            
            ! Paralelitzem nomes el LOOP INTERN
            j_0 = taskid * N / numproc + 1
            j_f = (taskid + 1) * N / numproc

            ! Cada task dona valor a la posició local
            pos_local = 0.d0
            do j = j_0, j_f
                do i = 1, D 
                    pos_local(i,j) = a * (mod(j - 1, M ** (D - i + 1)) / (M ** (D - i)))
                end do
            end do
            
            ! El master rep la info al vector pos
            if (taskid == master) pos = 0.d0
            do i = 1, D
                call MPI_Reduce(pos_local(i,:), pos(i,:), N, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierror)
            end do
            if (taskid == master) pos = pos - (L - a) / 2.d0  ! Centrem el sistema al (0,0,0)

        end subroutine init_sc_reduce

        subroutine init_vel_gather(vel, T)
            ! Author: Eloi Sanchez
            ! Torna el array de velocitats (aleatories) vel(D,N) consistent amb la T donada.
            ! --- VARIABLES ---
            ! vel(D,N) -> Array on es tornaran les velocitats de les part.
            !        T -> Temp. a la que s'inicialitzara la velocitat de les part.
            use parameters, only : D, N, taskid, numproc, master
            use integraforces, only : energy_kin
            implicit none
            include 'mpif.h'

            real*8, intent(inout) :: vel(D,N)
            real*8, intent(in) :: T

            real*8, allocatable :: vel_local(:,:)
            integer :: local_size, all_size(numproc), all_position(numproc)

            integer :: seed
          
            real*8 :: vel_CM_local(D), aux_CM(D), vel_CM(D)
            real*8 :: dummy_T, kin
            integer :: i, j, i_0, i_f, ierror
          
            ! Creem les variables locals de cada task
            i_0 = taskid * N / numproc + 1
            i_f = (taskid + 1) * N / numproc
            allocate(vel_local(D,i_0:i_f))
            local_size = i_f - i_0 + 1

            ! Fem una seed per cada task
            seed = int(MPI_Wtime() * 1000000 * (taskid * 2 + 1))
            ! print*, "taskid", taskid, "has seed", seed
            call srand(seed)

            ! Inicialitza les velocitats de manera random entre -1 i 1
            vel_CM_local = 0.d0
            do i = i_0, i_f
                do j = 1, D
                    vel_local(j,i) = 2.*rand() - 1.
                end do
              vel_CM_local = vel_CM_local + vel_local(:,i)
            end do
            vel_CM_local = vel_CM_local/dble(local_size)
            
            ! Eliminem la velocitat neta del sistema
            do i = i_0, i_f
                vel_local(:,i) = vel_local(:,i) - vel_CM_local
            end do
            
            call MPI_Gather(local_size, 1, MPI_INTEGER, all_size, 1, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierror)

            if (taskid == master) then
                all_position(1) = 0  ! Ha de començar al 0 per coses del OpenMPI
                do i = 1, numproc - 1
                    all_position(i+1) = all_position(i) + all_size(i)
                end do
            end if

            ! Es guarda al master la velocitat total a partir de les locals
            do i = 1, D
                call MPI_Gatherv(vel_local(i,:), local_size, MPI_DOUBLE_PRECISION, vel(i,:), &
                                all_size, all_position, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do
            
            ! Reescalem les velocitats a la temperatura objectiu
            call energy_kin(vel, kin, dummy_T)  
            do i = 1, D
                call MPI_Bcast(vel(i,:), N, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do
            call MPI_Bcast(kin, 1, MPI_DOUBLE_PRECISION, master, &
                            MPI_COMM_WORLD, ierror)
            
            do i = i_0, i_f
                vel_local(:,i) = vel(:,i) * sqrt(dble(3*N)*T/(2.d0*kin))
            end do

            do i = 1, D
                call MPI_Gatherv(vel_local(i,:), local_size, MPI_DOUBLE_PRECISION, vel(i,:), &
                                all_size, all_position, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do

        end subroutine

end module init
