   module integraforces
         use parameters
         use pbc
         contains

         subroutine compute_force_LJ(r,f,U,P)
         ! Computes the force, potential energy and pressure of a system
         ! of N particles. Needs the external parameters M,D,L,rc to work
               implicit none
               real*8,intent(in) :: r(D,N)
               real*8,intent(out) :: f(D,N),U,P
               real*8 :: distv(D),dist
               integer :: i,j
               f = 0.d0
               U = 0.d0
               P = 0.d0
               do i=1,N
                     do j=i+1,N
                           distv = r(:,i)-r(:,j)
                           call min_img_2(distv)
                           dist = sqrt(sum((distv)**2))
                           if(dist<rc) then
                                 f(:,i) = f(:,i)& 
                                 + (48.d0/dist**14 - 24.d0/dist**8)*distv
                                 f(:,j) = f(:,j)& 
                                 - (48.d0/dist**14 - 24.d0/dist**8)*distv
   
                                 U = U + 4.d0*((1.d0/dist)**12-(1.d0/dist)**6)-&
                                         4.d0*((1.d0/rc)**12-(1.d0/rc)**6)
                                 P = P + sum(distv * f(:,i))
                           endif
                     enddo
               enddo
               ! P = P/(dble(N)*(dble(N)-1.d0)/2.d0)

               P = 1.d0/(3.d0*L**3)*P
         end subroutine compute_force_LJ

         subroutine andersen_therm(v,dt,Temp)
         !     Applies the Andersen thermostat to a system of N particles. Requires
         !     boxmuller subroutine.
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
            real*8,intent(inout) :: v(D,N)
            real*8,intent(in) :: dt,Temp
            real*8 :: std,nu,x1,x2
            integer :: i,j
            std = sqrt(Temp) !Standard deviation of the gaussian.
            nu = 0.1*dt
            do i=1,N
                  if (rand()<nu) then ! Check if collision happens.
                        do j=1,D
                              x1 = rand()
                              x2 = rand()
                              call box_muller(std,x1,x2) !Modify the velocity.
                              v(i,j) = x1
                        enddo
                  endif
            enddo
         end subroutine andersen_therm  
                  

      subroutine box_muller(std, X1,X2)
            implicit none
            real*8,intent(in) :: std
            real*8,intent(out) :: x1,x2
            real*8 :: PI
            PI = 4d0*datan(1d0)
            x1 = std*dsqrt(-2d0*dlog(1.d0-x1))*dcos(2d0*PI*x2)
            x2 = std*dsqrt(-2d0*dlog(1.d0-x1))*dsin(2d0*PI*x2)
      end subroutine box_muller


         subroutine energy_kin(v,ekin,Tins)
         ! Computes the kinetic energy and instant temperature of the
         ! system
            implicit none
            real(8), intent(in) :: v(D,N)
            real(8), intent (out) :: ekin,Tins
            integer :: i
            
            ekin=0.0d0
            do i=1,N
               ekin=ekin+sum(v(:,i)**2)/2.0d0  
            enddo

            Tins=2.0d0*ekin/(3.0d0*N)
            return 
         end subroutine energy_kin


         subroutine verlet_v_step(r,v,t,dt)
         ! Computes one time step with velocity verlet algorithm 
            implicit none 
            integer :: i
            real(8) :: r(D,N), v(D,N), U, f(D,N), t, dt, P
            ! --------------------------------------------------
            ! input: r and v
            ! output: modified r and v
            ! --------------------------------------------------

            !forces at t
            call compute_force_LJ(r,f,U,P)

            !change of positions
            do i=1,N
               r(:,i)=r(:,i)+v(:,i)*dt +f(:,i)*dt**2/2.0d0
               !Apply PBC
               call min_img_2(r(:,i)) 

               v(:,i)=v(:,i)+f(:,i)*dt*0.5d0
            enddo

            !forces at t+dt
            call compute_force_LJ(r,f,U,P)

            do i=1,N
               v(:,i)=v(:,i)+f(:,i)*dt*0.5d0
            enddo 

            t=t+dt !Update time
            !AJ : this is not a good way to update time, prone to carry over errors.
            return
         end subroutine verlet_v_step


         subroutine vvel_solver(Nt,dt,r,v,Temp,eunit,eunit_g,flag_g) !Nt -- n_total?
   !     Performs Nt steps of the velocity verlet algorithm while computing
   !     different observables and writing to file.
   !     Set flag to different to a non-zero int to write files.
            use radial_distribution
            implicit none
            integer, intent(in) :: Nt, eunit, eunit_g
            real(8) :: dt, Temp
            real(8) :: r(D,N), v(D,N), f(D,N)
            real(8) :: ekin, U, t, Tins, Ppot, Ptot
            integer :: i,j
            ! Flags for writing:
            integer, intent(in) :: flag_g
            integer :: Nshells
            
            ! Initialization of the g(r) calculation:
            if(flag_g.ne.0) then
              Nshells = 100
              call prepare_shells(Nshells)
            endif

            t = 0.d0
            call compute_force_LJ(r,f,U,Ppot) !Initial force, energy and pressure
            call energy_kin(v,ekin,Tins) !Compute initial kinetic energy
            Ptot = rho*Tins + Ppot !AJ: modifed to use the rho in parameters.

            !Write intial results.
            write(eunit,*)"t","K","U","E","T","v_tot"
            write(eunit,*) t, ekin, U, ekin+U, Tins, sum(v,2)
            !do i=1,N   !Writes initial configuration of particles
              ! write(runit,*) r(:,i)
            !enddo
      
            do i=1,Nt !Main time loop.
               call verlet_v_step(r,v,t,dt) !Perform Verlet step.
               call andersen_therm(v,dt,Temp) !Apply thermostat
               call compute_force_LJ(r,f,U,Ppot)
               call energy_kin(v,ekin,Tins)
               Ptot = rho*Tins + Ppot
            !    call rad_distr_fun(r)

               !Write to file.
               write(eunit,*) t, ekin, U, ekin+U, Tins, sum(v,2)
               
               !Write snapshot of g(r)
               if(flag_g.ne.0) then
                 do j=1,Nshells
                     write(eunit_g,*) (j-1)*grid_shells, g(j)
                 enddo
                 write(eunit_g,*) ! separation line
               endif
               
            enddo

            return
         end subroutine vvel_solver


   end module integraforces
