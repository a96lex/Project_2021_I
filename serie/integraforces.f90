   module integraforces
         use parameters
         use pbc
         contains

         subroutine compute_force_LJ(r,f,U,P)
         !Author: Arnau Jurado
         ! Computes the force, potential energy and pressure of a system
         ! of N particles. Needs the external parameters D,rc to work.
         ! The pair-wise interaction is the LJ interaction potential.
         !           Input
         !           -----
         !                 r : real*8,dimension(D,N)
         !                       Positions of the particles of the system.
         !           Output
         !           ------
         !                 f : real*8,dimension(D,N)
         !                       Force on the particles. F(:,i) is the force
         !                       vector of the forces acting on the i*th particle.
         !                 U : real*8
         !                       Potential energy of the system.
         !                 P : real*8
         !                       Potential contribution of the pressure.
         !                       Factor 1/(3V) is already included on output
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
                           call min_img(distv)
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

         subroutine andersen_therm(v,Temp)
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
            real*8,intent(in) :: Temp
            real*8 :: PI,std,x1,x2
            integer :: i,j
            std = sqrt(Temp) !Standard deviation of the gaussian.
            PI = 4d0*datan(1d0)
            do i=1,N
                  if (rand()<nu) then ! Check if collision happens.
                        do j=1,D
                              x1 = rand()
                              x2 = rand()
                              v(j,i) = std*dsqrt(-2d0*dlog(1.d0-x1))*dcos(2d0*PI*x2)
                              !Here we are effectively throwing away one of the
                              !two gaussian numbers given by box_muller
                        enddo
                  endif
            enddo
         end subroutine andersen_therm  
                  


         subroutine energy_kin(v,ekin,Tins)
         !Author: Laia Barjuan
         ! Computes the kinetic energy and instant temperature of the
         ! system
            ! --------------------------------------------------
            !  input: 
            !          v  --> velocities of the system
            !  output: 
            !          ekin --> kinetic energy of the system 
            !          Tins --> instant temperature of the system
            ! --------------------------------------------------
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


         subroutine verlet_v_step(r,v,fold,t,time_i,dt,U,P)
         !Author: Laia Barjuan
         ! Computes one time step with velocity verlet algorithm 
            ! --------------------------------------------------
            ! input: 
            !        r --> positions of the particles
            !        v  --> velocities of the system
            !        fold --> force at t
            !        time_i --> number of step
            !        dt --> time step
            ! output: 
            !         modified r and v
            !         t --> final time
            !         U --> potential energy 
            !         P --> potential pressure
            ! --------------------------------------------------
            implicit none 
            integer :: i, time_i
            real(8) :: r(D,N), v(D,N), U, f(D,N), fold(D,N), t, dt, P

            !change of positions
            do i=1,N
               r(:,i)=r(:,i)+v(:,i)*dt +fold(:,i)*dt**2/2.0d0
               !Apply PBC
               call min_img(r(:,i)) 
            enddo

            !forces at t+dt
            call compute_force_LJ(r,f,U,P)

            !update of velocities
            do i=1,N
               v(:,i)=v(:,i)+(fold(:,i)+f(:,i))*dt*0.5d0
            enddo

            fold=f

            t=time_i*dt !Update time
            return
         end subroutine verlet_v_step


         subroutine vvel_solver(Nt,dt,r,v,Temp,eunit,eunit_dim,eunit_g,eunit_g_dim,flag_g)
         !Author: Laia Barjuan
         !Co-Authors: David March (radial distribution), Arnau Jurado (interface with forces)
   !     Performs Nt steps of the velocity verlet algorithm while computing
   !     different observables and writing to file.
            ! --------------------------------------------------
            ! input: 
            !        Nt --> number of time steps
            !        dt  --> time step
            !        r --> positions of the particles
            !        v  --> velocities of the system
            !        T --> temperature of the system
            !        eunit --> unit of file to write energies and temperature
            !        eunit_dim --> "" with physical units
            !        eunit_g --> unit of file to write g(r)
            !        flag_g --> different to a non-zero int to write files
            ! output: 
            !         modified r and v
            ! --------------------------------------------------
            use rad_dist
            implicit none
            integer, intent(in) :: Nt, eunit, eunit_dim, eunit_g, eunit_g_dim
            real(8) :: dt, Temp
            real(8) :: r(D,N), v(D,N), f(D,N), fold(D,N)
            real(8) :: ekin, U, t, Tins, Ppot, Ptot
            integer :: i,j
            ! Flags for writing g:
            integer, intent(in) :: flag_g
            integer :: Nshells
            real(8), dimension(:), allocatable :: g_avg, g_squared_avg
            
            ! Initialization of the g(r) calculation:
            if(flag_g.ne.0) then
              Nshells = 100
              call prepare_shells(Nshells)
              allocate(g_avg(Nshells))
              allocate(g_squared_avg(Nshells))
              g_avg = 0d0
              g_squared_avg = 0d0
            endif

            t = 0.d0
            call compute_force_LJ(r,f,U,Ppot) !Initial force, energy and pressure
            call energy_kin(v,ekin,Tins) !Compute initial kinetic energy
            Ptot = rho*Tins + Ppot !AJ: modifed to use the rho in parameters.

            !Write intial results.
            write(eunit,*)"#t,   K,   U,  E,  T,  v_tot,  Ptot"
            write(eunit,*) t, ekin, U, ekin+U, Tins, dsqrt(sum(sum(v,2)**2)), Ptot

            write(eunit_dim,*)"#t,   K,   U,  E,  T,  Ptot"
            write(eunit_dim,*) t*unit_of_time,&
                        ekin*unit_of_energy, U*unit_of_energy, (ekin+U)*unit_of_energy,&
                        Tins*epsilon, Ptot*unit_of_pressure

            fold = f
            do i=1,Nt !Main time loop.
               call verlet_v_step(r,v,fold,t,i,dt,U,Ppot) !Perform Verlet step.
               call andersen_therm(v,Temp) !Apply thermostat
               call energy_kin(v,ekin,Tins)
               Ptot = rho*Tins + Ppot

               !Write to files.
               write(eunit,*) t, ekin, U, ekin+U, Tins, dsqrt(sum(sum(v,2)**2)), Ptot
               write(eunit_dim,*) t*unit_of_time,&
                        ekin*unit_of_energy, U*unit_of_energy, (ekin+U)*unit_of_energy,&
                        Tins*epsilon, Ptot*unit_of_pressure
               !Save snapshot of g(r) to average
               if(flag_g.ne.0) then
                 call rad_distr_fun(r,Nshells)
                 g_avg = g_avg + g
                 g_squared_avg = g_squared_avg + g**2
               endif

               if(mod(i,int(0.001*Nt))==0) then
                  write (*,"(A,F5.1,A)",advance="no") "Progress: ",i/dble(Nt)*100.,"%"
                  if (i<Nt) call execute_command_line('echo "\033[A"')
                endif            
            enddo
            
            if(flag_g.ne.0) then
              do j=1,Nshells
                write(eunit_g,*) (j-1)*grid_shells+grid_shells/2d0, g(j), dsqrt(g_squared_avg(i) - g_avg(i)**2)
                write(eunit_g_dim,*) ((j-1)*grid_shells+grid_shells/2d0)*sigma, g(j), dsqrt(g_squared_avg(i) - g_avg(i)**2)
              enddo
            endif

            return
         end subroutine vvel_solver


   end module integraforces
