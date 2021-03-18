module integraforces
      use parameters
   !    use pbc
      contains

      subroutine compute_force_LJ(r,f,U,P)
      !Author: Arnau Jurado
      ! Computes the force, potential energy and pressure of a system
      ! of N particles. Needs the external parameters M,D,L,rc to work
            implicit none
            include 'mpif.h'
            real*8,intent(in)  :: r(D,N)
            real*8,intent(out) :: f(D,N),U,P
            real*8             :: flocal(D,N),Ulocal,Plocal
            real*8,allocatable :: fi(:,:)
            real*8             :: distv(D),dist
            real*8             :: rlocal(D,N),coord(N)
            integer            :: i,j
            integer            :: imin,imax,particles_per_proc

            integer            :: ierror,request

            !Initialize quantities.
            f = 0.d0
            U = 0.d0
            P = 0.d0

            flocal = 0.d0
            Ulocal = 0.d0
            Plocal = 0.d0
            rlocal = r

            !If divisible everything ok
            if(mod(N,numproc)==0) then
               particles_per_proc = N/numproc
               imin = (particles_per_proc * taskid) + 1
               imax = (particles_per_proc *(taskid+1))
            else 
               !If not last cpu gets less particles, this gets more inefficient
               !the more cpus you have
               particles_per_proc = ceiling(N/real(numproc))
               if(taskid /= numproc-1) then
                  imin = (particles_per_proc * taskid) + 1
                  imax = (particles_per_proc *(taskid+1))
               else
                  imin = (particles_per_proc * taskid) + 1
                  imax = N
               end if
            end if
            ! print*,"taskid:",taskid,particles_per_proc
            ! print*,"taskid:",taskid,imin,imax

            allocate(fi(D,particles_per_proc)) !Not used in first version of this
            !subroutine

            do i=1,D
               coord = rlocal(i,:)
               call MPI_BCAST(coord,N,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
               rlocal(i,:) = coord
            end do


            do i=imin,imax !Loop over assigned particles
                  do j=1,N !We do the full double loop
                     if(i /= j) then
                        !Compute distance and apply minimum image convention.
                        distv = rlocal(:,i)-rlocal(:,j)
                     !    call min_img_2(distv) !Commented until pbc is added again
                        dist = sqrt(sum((distv)**2))


                        if(dist<rc) then !Cutoff
                              !Compute forces and pressure.
                              flocal(:,i) = flocal(:,i)& 
                              + (48.d0/dist**14 - 24.d0/dist**8)*distv

                              Ulocal = Ulocal + 4.d0*((1.d0/dist)**12-(1.d0/dist)**6)-&
                                      4.d0*((1.d0/rc)**12-(1.d0/rc)**6)
                              Plocal = Plocal + sum(distv * flocal(:,i))
                        end if
                     end if
                  end do
            end do
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            print*,taskid,flocal(:,1)

            do i=1,D
               coord = flocal(i,:)
               call MPI_REDUCE(coord,f(i,:),N,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            end do

            call MPI_REDUCE(Ulocal,U,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(Plocal,P,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)

            if(taskid==master) then
               U = U/2.d0
               ! P = P/2.d0
               P = P/(dble(N)*(dble(N)-1.d0))
               !Add 1/3V factor to potential pressure.
               P = 1.d0/(3.d0*L**3)*P
            end if
      end subroutine compute_force_LJ



!      subroutine energy_kin(v,ekin,Tins)
!      !Author: Laia Barjuan
!      ! Computes the kinetic energy and instant temperature of the
!      ! system
!         implicit none
!         real(8), intent(in)   :: v(D,N), vlocal(D,N), vec(N)
!         real(8), intent (out) :: ekin, ekinlocal, Tins
!         integer :: i
!         integer :: part_per_proc, imin, imax
!         integer :: reques, ierror
!         ! --------------------------------------------------
!           input: v  --> velocities of the system
!           output: ekin --> kinetic energy of the system 
!                   Tins --> instant temperature of the system
!         ! --------------------------------------------------
!
!         !Initialization
!         ekin=0.0d0
!         ekinlocal=0.d0
!
!         !Assign number of particles to processor
!         part_per_proc=int(dble(N)/dble(numproc)) 
!             !Warning: leaves behind particles if the division is not exact
!         imin=taskid*part_per_proc + 1
!         imax=(taskid+1)*part_per_proc
!
!         do i=1,D
!            if (taskid==master) vec = v(i,:)
!            call MPI_BCAST(vec,N,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
!            vlocal(i,:) = vec
!         end do
!
!         do i=imin,imax
!            ekinlocal=ekinlocal+sum(vlocal(:,i)**2)/2.0d0  
!         enddo
!
!         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!         call MPI_REDUCE(ekinlocal,ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
!      
!         !Store instant temperature information in master
!         if (taskid==master) Tins=2.0d0*ekin/(3.0d0*N)
!
!         return 
!      end subroutine energy_kin


end module integraforces
