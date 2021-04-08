module init
    use parameters
    implicit none
    contains

        subroutine get_param(unit)
            !Author: Eloi Sanchez
            ! Llegim de l'input els parametres del sistema i es calcula el
            ! nº d'iteracions i la longitud de la cel·la.
            ! Els propers llocs on es faci servir 'use parameters' tindran
            ! les variables acualitzades
            implicit none
            integer, intent(in):: unit
            integer :: errstat
            namelist /input/ N, D, rho, dt_sim, n_meas, n_conf, n_equil, T_ref, &
                            fact_rc, nu, sigma, epsilon, mass, seed

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

            seed = seed + taskid
            call srand(seed)

            call reduced_units()
            call divide_particles()
            call divide_particles_pairs()
            call divide_particles_pairs_improv()
        end subroutine get_param

        subroutine divide_particles()
           !Author: Arnau Jurado & Eloi Sanchez
           !Divides the work among the processors by assigning each one an "imin" and
           !a "imax", which are the indexes of the first and last particle they have
           !to process e.g. with forces, each processor computes the forces
           !from the imin-th particle to the imax-th particles, both included.
           !aux_pos and aux_size are stored in master to be used in Gatherv calls
            implicit none
            include 'mpif.h'
            integer :: i
            integer :: ierror
            
            imin = taskid * N / numproc + 1
            imax = (taskid + 1) * N / numproc
            local_size = imax - imin + 1
            
            ! We create aux_size(numproc) and aux_pos(numproc) only in master
            call MPI_Gather(local_size, 1, MPI_INTEGER, aux_size, 1, MPI_INTEGER, &
                            master, MPI_COMM_WORLD, ierror)
            
            if (taskid == master) then
                aux_pos(1) = 0  ! Must start at 0 for OpenMPI issues
                do i = 1, numproc - 1
                    aux_pos(i+1) = aux_pos(i) + aux_size(i)
                end do
            end if
            
            call MPI_Bcast(aux_pos, numproc, MPI_INTEGER, master, &
                          MPI_COMM_WORLD, ierror)
            call MPI_Bcast(aux_size, numproc, MPI_INTEGER, master, &
                          MPI_COMM_WORLD, ierror)

        end subroutine divide_particles
        
        subroutine divide_particles_pairs()
           ! Author: David March
           ! Distribute particles so they each processor computes an approx. equal number of pairs in a nested loop such as:
           ! do i=imin_p,imax_p
           !    do j=i+1,N
           ! Sets the particles ranges per processor in imin_p, imax_p
           implicit none
           !include 'mpif.h'
           integer :: i,j
           integer, dimension(N) :: num_pairs
           integer, dimension(:,:), allocatable :: ranges_proc
           real(8) total_pairs, pairs_per_proc, sum_pairs
           
           do i=1,N
              num_pairs(i) = N-i
           enddo
           total_pairs = dble(N*(N-1))/dble(2)
           pairs_per_proc = total_pairs/dble(numproc)
           
           allocate(ranges_proc(numproc,2)) ! inferior limit at (:,1), superior at (:,2) for each processor
           ranges_proc(1,1) = 1
           ranges_proc(numproc,2) = N-1
           do i=1,numproc-1
              sum_pairs = 0d0
              limits: do j=ranges_proc(i,1),N
                 sum_pairs = sum_pairs + dble(num_pairs(j))
                 if(sum_pairs.gt.pairs_per_proc) then
                    ranges_proc(i,2) = j
                    ranges_proc(i+1,1) = j+1
                    exit limits
                 endif
              enddo limits
           enddo
          
          ! Finally, assignate the min and max index to the global variables:
          imin_p = ranges_proc(taskid+1,1)
          imax_p = ranges_proc(taskid+1,2)
          !print*, "task ",taskid, " with particle ranges ", imin_p, imax_p
          deallocate(ranges_proc)
       end subroutine divide_particles_pairs
       
       subroutine divide_particles_pairs_improv()
           ! Author: David March
           ! Distribute particles so they each processor computes an approx. equal number of pairs in a nested loop such as:
           ! do i=imin_p,imax_p
           !    do j=jmin_p(i),jmax_p(i)
           ! Sets the particles ranges per processor in imin_p, imax_p; Each processor for will have jmin_p(i),jmax_p(i)
           ! for their working particles
           implicit none
           !include 'mpif.h'
           integer :: i,j,k,count_pairs,processor
           integer, dimension(numproc) :: track_pairs
           integer, dimension(N) :: num_pairs
           integer, dimension(:,:), allocatable :: ranges_proc_i, ranges_proc_j_imin, ranges_proc_j_imax
           real(8) total_pairs, pairs_per_proc, sum_pairs
           
           allocate(ranges_proc_i(numproc,2)) ! (:,1) for min i, (:,2) for max i
           allocate(ranges_proc_j_imin(numproc,2)) ! (:,1) for min j from min i, (:,2) for min j from max i
           allocate(ranges_proc_j_imax(numproc,2)) ! (:,1) for max j from min i, (:,2) for max j from max i
           
           count_pairs = 0
           total_pairs = N*(N-1)/2
           pairs_per_proc = dble(total_pairs)/dble(numproc)
           
           processor = 1
           ranges_proc_i(processor,1) = 1
           ranges_proc_j_imin(processor,1) = 2
           do i=1,N-1
             do j=i+1,N
               if(j.eq.N.and.i.eq.ranges_proc_i(processor,1)) ranges_proc_j_imin(processor,2) = N
               count_pairs = count_pairs + 1
               if(dble(count_pairs).ge.pairs_per_proc.or.(i.eq.(N-1).and.(j.eq.N))) then
                   ! Final indexes from processor
                   ranges_proc_i(processor,2) = i
                   ranges_proc_j_imax(processor,1) = i+1
                   ranges_proc_j_imax(processor,2) = j
                   track_pairs(processor) = count_pairs
                   count_pairs = 0
                   ! Start new processor, saving starting indexes:
                   if(processor.lt.numproc) then
                       processor = processor + 1
                       if(j.lt.N) then
                           ranges_proc_i(processor,1) = i
                           ranges_proc_j_imin(processor,1) = j+1
                       else
                           ranges_proc_i(processor,1) = i+1
                           ranges_proc_j_imin(processor,1) = i+2
                       endif
                       
                   endif
               endif
             enddo
           enddo
          
          allocate(jmin_p(ranges_proc_i(taskid+1,1):ranges_proc_i(taskid+1,2)))
          allocate(jmax_p(ranges_proc_i(taskid+1,1):ranges_proc_i(taskid+1,2)))
          imin_p = ranges_proc_i(taskid+1,1)
          imax_p = ranges_proc_i(taskid+1,2)
          jmin_p(ranges_proc_i(taskid+1,1)) = ranges_proc_j_imin(taskid+1,1)
          jmax_p(ranges_proc_i(taskid+1,1)) = ranges_proc_j_imin(taskid+1,2)
          jmin_p(ranges_proc_i(taskid+1,2)) = ranges_proc_j_imax(taskid+1,1)
          jmax_p(ranges_proc_i(taskid+1,2)) = ranges_proc_j_imax(taskid+1,2)
          ! In between, go normaly from i+1 to N
          do i=ranges_proc_i(taskid+1,1)+1,ranges_proc_i(taskid+1,2)-1
              jmin_p(i) = i+1
              jmax_p(i) = N
          enddo
          
          !if(taskid==master) CALL sleep(5)
          ! CHECK:
          !write(*,*) "Task", taskid
          !do i=1,N
          !    write(*,*) jmin_p(i),jmax_p(i)
          !enddo
          
          ! Prova serie per veure els rangs ben ordenats
          !if(taskid==master) then
          !do k=1,numproc
          !    write(*,'(A17, I2)') "Ranges processor ", k
          !    write(*,'(A8, 2(I4,1X))') "     i: ", ranges_proc_i(k,1),ranges_proc_i(k,2)
          !    write(*,'(A21, I5)') "     Pairs assigned: ", track_pairs(k)
          !    write(*,'(A20, I4, A12, 2(I4,1X))') "          For min i ", ranges_proc_i(k,1), " j range is: ", &
          !                                                   ranges_proc_j_imin(k,1), ranges_proc_j_imin(k,2)
          !    write(*,'(A20, I4, A12, 2(I4,1X))') "          For max i ", ranges_proc_i(k,2), " j range is: ", &
          !                                                   ranges_proc_j_imax(k,1), ranges_proc_j_imax(k,2)
          !enddo
          !endif
          deallocate(ranges_proc_i)
          deallocate(ranges_proc_j_imin)
          deallocate(ranges_proc_j_imax)
       end subroutine divide_particles_pairs_improv

        subroutine init_sc(pos)
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

            implicit none
            include 'mpif.h'
            real*8, intent(out):: pos(D,N)
            integer :: M    ! Nº atoms en cada dimensio.
            real*8 :: a     ! Distancia interatomica i variable per assignar posicions al loop
            integer :: i, j
            real*8 :: M_check
            real*8, allocatable:: pos_local(:,:)
            integer :: ierror

            ! Aixo ho fan totes les tasks
            M_check = N ** (1.d0 / D)
            M = ceiling(M_check)
            a = L / dble(M)
            allocate(pos_local(D,imin:imax))

            ! Cada task dona valor a la posició local
            pos_local = 0.d0
            do i = imin, imax
                do j = 1, D 
                    pos_local(j,i) = a * (mod(i - 1, M ** (D - j + 1)) / M ** (D - j))
                end do
            end do
            pos_local = pos_local - (L - a) / 2.d0  ! Centrem el sistema al (0,0,0)

            ! Es guarda al master la posicio total a partir de les locals
            do i = 1, D
                call MPI_Allgatherv(pos_local(i,:), local_size, MPI_DOUBLE_PRECISION, pos(i,:), &
                                aux_size, aux_pos, MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD, ierror)
            end do

            deallocate(pos_local)           
        end subroutine init_sc

        subroutine init_vel(vel, T)
            ! Author: Eloi Sanchez
            ! Torna el array de velocitats (aleatories) vel(D,N) consistent amb la T donada.
            ! --- VARIABLES ---
            ! vel(D,N) -> Array on es tornaran les velocitats de les part.
            !        T -> Temp. a la que s'inicialitzara la velocitat de les part.
            use integraforces, only : energy_kin
            implicit none
            include 'mpif.h'

            real*8, intent(inout) :: vel(D,N)
            real*8, intent(in) :: T

            real*8, allocatable :: vel_local(:,:)
          
            real*8 :: vel_CM_local(D)
            real*8 :: dummy_T, kin
            integer :: i, j, ierror
          
            allocate(vel_local(D,imin:imax))

            ! Inicialitza les velocitats de manera random entre -1 i 1
            vel_CM_local = 0.d0
            do i = imin, imax
                do j = 1, D
                    vel_local(j,i) = 2.*rand() - 1.
                end do
              vel_CM_local = vel_CM_local + vel_local(:,i)
            end do
            vel_CM_local = vel_CM_local/dble(local_size)
            
            ! Eliminem la velocitat neta del sistema
            do i = imin, imax
                vel_local(:,i) = vel_local(:,i) - vel_CM_local
            end do

            ! Es guarda al master la velocitat total a partir de les locals
            do i = 1, D
                call MPI_Gatherv(vel_local(i,:), local_size, MPI_DOUBLE_PRECISION, vel(i,:), &
                                aux_size, aux_pos, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do
            
            ! Reescalem les velocitats a la temperatura objectiu
            call energy_kin(vel, kin, dummy_T)  
            call MPI_Bcast(kin, 1, MPI_DOUBLE_PRECISION, master, &
                            MPI_COMM_WORLD, ierror)
            
            do i = imin, imax
                vel_local(:,i) = vel_local(:,i) * sqrt(dble(3*N)*T/(2.d0*kin))
            end do

            do i = 1, D
                call MPI_Gatherv(vel_local(i,:), local_size, MPI_DOUBLE_PRECISION, vel(i,:), &
                                aux_size, aux_pos, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do

        end subroutine init_vel

        subroutine init_vel_bad(vel, T)
            ! Author: Eloi Sanchez
            ! Torna el array de velocitats (aleatories) vel(D,N) consistent amb la T donada.
            ! --- VARIABLES ---
            ! vel(D,N) -> Array on es tornaran les velocitats de les part.
            !        T -> Temp. a la que s'inicialitzara la velocitat de les part.
            use integraforces, only : energy_kin
            implicit none
            include 'mpif.h'

            real*8, intent(inout) :: vel(D,N)
            real*8, intent(in) :: T

            real*8, allocatable :: vel_local(:,:)
          
            real*8 :: vel_CM_local(D)
            real*8 :: dummy_T, kin
            integer :: i, j, ierror
          
            allocate(vel_local(D,imin:imax))

            ! Inicialitza les velocitats de manera random entre -1 i 1
            vel_CM_local = 0.d0
            do i = imin, imax
                do j = 1, D
                    vel_local(j,i) = 2.*rand() - 1.
                end do
              vel_CM_local = vel_CM_local + vel_local(:,i)
            end do
            vel_CM_local = vel_CM_local/dble(local_size)
            
            ! Eliminem la velocitat neta del sistema
            do i = imin, imax
                vel_local(:,i) = vel_local(:,i) - vel_CM_local
            end do

            ! Es guarda al master la velocitat total a partir de les locals
            do i = 1, D
                call MPI_Gatherv(vel_local(i,:), local_size, MPI_DOUBLE_PRECISION, vel(i,:), &
                                aux_size, aux_pos, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do
            
            ! Reescalem les velocitats a la temperatura objectiu
            call energy_kin(vel, kin, dummy_T)  
            call MPI_Bcast(kin, 1, MPI_DOUBLE_PRECISION, master, &
                            MPI_COMM_WORLD, ierror)
            
            do i = imin, imax
                vel_local(:,i) = vel_local(:,i) * sqrt(dble(3*N)*T/(2.d0*kin))
            end do

            do i = 1, D
                call MPI_Gatherv(vel_local(i,:), local_size, MPI_DOUBLE_PRECISION, vel(i,:), &
                                aux_size, aux_pos, MPI_DOUBLE_PRECISION, master, &
                                MPI_COMM_WORLD, ierror)
            end do

        end subroutine init_vel_bad
end module init
