   module integraforces
         use parameters
         use pbc
         contains

         subroutine compute_force_LJ(r,f,U,P)
         !Author: Arnau Jurado
         ! Computes the force, potential energy and pressure of a system
         ! of N particles. Needs the external parameters M,D,L,rc to work
               implicit none
               real*8,intent(in) :: r(D,N)
               real*8,intent(out) :: f(D,N),U,P
               real*8 :: distv(D),dist,fij(D)
               integer :: i,j
               !Initialize quantities.
               f = 0.d0
               U = 0.d0
               P = 0.d0
               do i=1,N !loop over i<j particles
                     do j=i+1,N
                           !Compute distance and apply minimum image convention.
                           distv = r(:,i)-r(:,j)
                           ! call min_img_2(distv)
                           dist = sqrt(sum((distv)**2))


                           if(dist<rc) then !Cutoff
                                 !Compute forces and pressure.
                                 fij = (48.d0/dist**14 - 24.d0/dist**8)*distv
                                 f(:,i) = f(:,i) + fij
                                 f(:,j) = f(:,j) - fij
   
                                 U = U + 4.d0*((1.d0/dist)**12-(1.d0/dist)**6)-&
                                         4.d0*((1.d0/rc)**12-(1.d0/rc)**6)
                                 P = P + sum(distv * fij)
                           endif
                     enddo
               enddo
               !Add 1/3V factor to potential pressure.
               P = 1.d0/(3.d0*L**3)*P
         end subroutine compute_force_LJ

         subroutine andersen_therm(v,dt,Temp)
         !Author: Arnau Jurado
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
            nu = 0.1
            do i=1,N
                  if (rand()<nu) then ! Check if collision happens.
                        do j=1,D
                              x1 = rand()
                              x2 = rand()
                              call box_muller(std,x1,x2) !Modify the velocity.
                              v(i,j) = x1
                              !Here we are effectively throwing away one of the
                              !two gaussian numbers given by box_muller
                        enddo
                  endif
            enddo
         end subroutine andersen_therm  
                  

      subroutine box_muller(std, X1,X2)
         !Author: Arnau Jurado
         !Generates two gaussian distributed numbers from two uniform random numbers
         !using the box muller method.
            implicit none
            real*8,intent(in) :: std
            real*8,intent(inout) :: x1,x2
            real*8 :: x1aux,x2aux
            real*8 :: PI
            PI = 4d0*datan(1d0)
            x1aux = x1
            x2aux = x2
            x1 = std*dsqrt(-2d0*dlog(1.d0-x1aux))*dcos(2d0*PI*x2aux)
            x2 = std*dsqrt(-2d0*dlog(1.d0-x1aux))*dsin(2d0*PI*x2aux)
      end subroutine box_muller


         subroutine energy_kin(v,ekin,Tins)
         !Author: Laia Barjuan
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


         subroutine verlet_v_step(r,v,t,time_i,dt,U,P)
         !Author: Laia Barjuan
         ! Computes one time step with velocity verlet algorithm 
            implicit none 
            integer :: i, time_i
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

            t=(time_i-1)*dt !Update time
            return
         end subroutine verlet_v_step


         subroutine vvel_solver(Nt,dt,r,v,Temp,eunit,eunit_g,flag_g)
         !Author: Laia Barjuan
         !Co-Authors: David March (radial distribution), Arnau Jurado (interface with forces)
   !     Performs Nt steps of the velocity verlet algorithm while computing
   !     different observables and writing to file.
   !     Set flag to different to a non-zero int to write files.
            use rad_dist
            implicit none
            integer, intent(in) :: Nt, eunit, eunit_g
            real(8) :: dt, Temp
            real(8) :: r(D,N), v(D,N), f(D,N)
            real(8) :: ekin, U, t, Tins, Ppot, Ptot
            integer :: i,j
            ! Flags for writing g:
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
      
            do i=1,Nt !Main time loop.
               call verlet_v_step(r,v,t,i,dt,U,Ppot) !Perform Verlet step.
               call andersen_therm(v,dt,Temp) !Apply thermostat
               call compute_force_LJ(r,f,U,Ppot)
               call energy_kin(v,ekin,Tins)
               Ptot = rho*Tins + Ppot

               !Write to file.
               write(eunit,*) t, ekin, U, ekin+U, Tins, sum(v,2)
               
               !Write snapshot of g(r)
               if(flag_g.ne.0) then
                 call rad_distr_fun(r,Nshells)
                 do j=1,Nshells
                     write(eunit_g,*) (j-1)*grid_shells, g(j)
                 enddo
                 write(eunit_g,*) ! separation line
               endif
               
            enddo

            return
         end subroutine vvel_solver


   end module integraforces
