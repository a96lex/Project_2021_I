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
               flocal = 0.d0
               U = 0.d0
               P = 0.d0

               particles_per_proc = int(N/float(numproc))
               imin = (particles_per_proc* taskid)     + 1
               imax = (particles_per_proc*(taskid+1))

               allocate(fi(D,particles_per_proc)) !Not used in first version of this
               !subroutine

               do i=1,D
                  coord = rlocal(i,:)
                  call MPI_BCAST(coord,N,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
                  rlocal(i,:) = coord
               end do


               do i=imin,imax !Loop over assigned particles
                     do j=1,N !We do the full double loop
                           !Compute distance and apply minimum image convention.
                           distv = rlocal(:,i)-rlocal(:,j)
                        !    call min_img_2(distv) !Commented until pbc is added again
                           dist = sqrt(sum((distv)**2))


                           if(dist<rc) then !Cutoff
                                 !Compute forces and pressure.
                                 flocal(:,i) = flocal(:,i)& 
                                 + (48.d0/dist**14 - 24.d0/dist**8)*distv
                                 flocal(:,j) = flocal(:,j)& 
                                 - (48.d0/dist**14 - 24.d0/dist**8)*distv
   
                                 Ulocal = Ulocal + 4.d0*((1.d0/dist)**12-(1.d0/dist)**6)-&
                                         4.d0*((1.d0/rc)**12-(1.d0/rc)**6)
                                 Plocal = Plocal + sum(distv * f(:,i))
                           end if
                     end do
               end do
               call MPI_BARRIER(MPI_COMM_WORLD,ierror)

               do i=1,N
                  do j=1,D
                     call MPI_REDUCE(flocal(j,i),f(j,i),1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
                  end do
               end do

               call MPI_REDUCE(Ulocal,U,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
               call MPI_REDUCE(Plocal,P,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)

               if(taskid==master) then
                  P = P/(dble(N)*(dble(N)-1.d0)/2.d0)
                  !Add 1/3V factor to potential pressure.
                  P = 1.d0/(3.d0*L**3)*P
               end if
         end subroutine compute_force_LJ

         
         subroutine andersen_therm(velocities,temperature)
            !Author: Arnau Jurado
            !     Applies the Andersen thermostat to a system of N particles. 
            !           Input
            !           -----
            !                 v : real*8,dimension(D,N)
            !                       Velocities of the system.
            !                 T : real*8
            !                       Target temperature of the thermostat.
            !           Output
            !           ------
            !                 v : real*8,dimension(D,N)
            !                       Velocities of the system after the thermostat application.
               implicit none
               include 'mpif.h'
               real*8,intent(inout) :: velocities(D,N)
               real*8,intent(in) :: temperature
               real*8 :: std,nu,x1,x2,pi
               integer :: i,j

               if (taskid.eq.master) then
                  do i=1,N
                     do j=1,D
                        call MPI_SEND(velocities(i,j),1,MPI_DOUBLE_PRECISION,mod(i,numproc),master,MPI_COMM_WORLD,ierror)
                     enddo
                  enddo
               endif
               
               pi = 4d0*datan(1d0)
               std = sqrt(temperature) !Standard deviation of the gaussian.
               nu = 0.1


               do i=1,numproc
                  do i=1,N
                        if (rand()<nu) then ! Check if a collision happens.
                              do j=1,D
                                    call MPI_RECEIVE(velocities(i,j), 1, MPI_DOUBLE_PRECISION, mod(i,numproc), 1, MPI_COMM_WORLD, stat, ierror)
                                    x1 = rand()
                                    x2 = rand()
                                    velocities(i,j) = std*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*pi*x2)
                              enddo
                        endif
                  enddo
               enddo

               call MPI_BARRIER(MPI_COMM_WORLD,ierror)

            end subroutine andersen_therm  
                     
   end module integraforces
