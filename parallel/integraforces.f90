module integraforces
      use parameters
      use pbc
      contains

      subroutine compute_force_LJ(r,f,U,P)
      !Author: Arnau Jurado
      ! Computes the force, potential energy and pressure of a system
      ! of N particles. Needs the external parameters D,rc to work.
      ! The pair-wise interaction is the LJ interaction potential.
      ! Parallel version, on calling the subroutine all tasks should know all
      ! the updated positions, r, of the particles.
      !           Input
      !           -----
      !                 r : real*8,dimension(D,N)
      !                       Positions of the particles of the system.
      !           Output
      !           ------
      !                 f : real*8,dimension(D,N)
      !                       Force on the particles. F(:,i) is the force
      !                       vector of the forces acting on the i*th particle.
      !                       Each task only has the force for their assigned particles.                      
      !                 U : real*8
      !                       Potential energy of the system. Master has 
      !                       the sum of all contributions.
      !                 P : real*8
      !                       Potential contribution of the pressure.
      !                       Factor 1/(3V) is already included on output. 
      !                       Master has the sum of all contributions.
            implicit none
            include 'mpif.h'
            real*8,intent(in)  :: r(D,N)
            real*8,intent(out) :: f(D,N),U,P
            real*8             :: Ulocal,Plocal
            real*8             :: distv(D),dist,fij(D)
            integer            :: i,j,ierror

            !Initialize quantities.
            f = 0.d0
            U = 0.d0
            Ulocal = 0.d0
            P = 0.d0
            Plocal = 0.d0

            do i=imin,imax !Loop over assigned particles
                  do j=1,N
                        if(i/=j) then
                              !Compute distance and apply minimum image convention.
                              distv = r(:,i)-r(:,j)
                              call min_img(distv)
                              dist = sqrt(sum((distv)**2))

                              if(dist<rc) then !Cutoff
                                    !Compute forces and pressure.
                                    fij = (48.d0/dist**14 - 24.d0/dist**8)*distv
                                    f(:,i) = f(:,i) + fij

                                    Ulocal = Ulocal + 4.d0*((1.d0/dist)**12-(1.d0/dist)**6)-&
                                                4.d0*((1.d0/rc)**12-(1.d0/rc)**6)
                              end if
                        end if
                  end do
                  Plocal = Plocal + sum(r(:,i) * f(:,i))
            end do
            !Add 1/3V factor and collect contributions of U and P on Master.
            Plocal = 1.d0/(3.d0*L**3)*Plocal
            call MPI_REDUCE(Ulocal,U,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            if(taskid==master) U = U/2.d0 !Remove double counting
            call MPI_REDUCE(Plocal,P,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
      end subroutine compute_force_LJ

      subroutine andersen_therm(v,Temp)
         !Author: Alex Parraga
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
            include 'mpif.h'
            real*8,intent(inout) :: v(D,N)
            real*8,intent(in) :: Temp
            real*8 :: std,x1,x2,PI
            integer :: i,j

            std = sqrt(Temp) !Standard deviation of the gaussian.
            PI = 4d0*datan(1d0)

            do i=imin,imax
               if (rand()<nu) then ! Check if collision happens.
                  do j=1,D
                     x1 = rand()
                     x2 = rand()
                     v(j,i) = std*dsqrt(-2d0*dlog(1.d0-x1))*dcos(2d0*PI*x2)
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
        include 'mpif.h'
        real(8), intent(in)   :: v(D,N)
        real(8) :: ekinlocal
        real(8), intent (out) :: ekin, Tins
        integer :: i, ierror

        !Initialization
        ekin=0.0d0
        ekinlocal=0.d0

        do i=imin,imax
           ekinlocal=ekinlocal+sum(v(:,i)**2)/2.0d0  
        enddo

        call MPI_REDUCE(ekinlocal,ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
     
        !Store instant temperature information in master
        if (taskid==master) Tins=2.0d0*ekin/(3.0d0*N)

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
         include 'mpif.h'
         integer :: i, time_i
         real(8) :: r(D,N), rlocal(D,N), v(D,N), U, f(D,N), fold(D,N), t, dt, P
         integer :: ierror


         !Initialize variables
         rlocal = r

         !Change of positions
         do i=imin,imax
            rlocal(:,i)=rlocal(:,i)+v(:,i)*dt +fold(:,i)*dt**2*0.5d0
            call min_img(rlocal(:,i))  !Apply PBC
         enddo

         do i=1,D
            call MPI_Allgatherv(rlocal(i,imin:imax), local_size, MPI_DOUBLE_PRECISION, r(i,:), &
                                aux_size, aux_pos, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierror)
         enddo

         !forces at t+dt
         call compute_force_LJ(r,f,U,P)

         !Compute velocities
         do i=imin,imax
            v(:,i)=v(:,i)+(fold(:,i)+f(:,i))*dt*0.5d0
         enddo 
        
         fold=f
         t=time_i*dt !Update time
         return
      end subroutine verlet_v_step


      subroutine vvel_solver(Nt,dt,r,v,Temp,eunit,eunit_dim,eunit_g,eunit_g_dim,flag_g)
      !Author: Laia Barjuan
      !Co-Authors: David March (radial distribution), Arnau Jurado (interface
      !with forces)
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
         include 'mpif.h'
         integer, intent(in) :: Nt, eunit, eunit_dim, eunit_g, eunit_g_dim
         real(8) :: dt, Temp
         real(8) :: r(D,N), v(D,N), f(D,N), fold(D,N)
         real(8) :: ekin, U, t, Tins, Ppot, Ptot
         integer :: i,j
         ! Flags for writing g:
         integer, intent(in) :: flag_g
         integer :: Nshells, ierror
         real(8), dimension(:), allocatable :: g_avg, g_avg_final
         
         ! Initialization of the g(r) calculation:
         if(flag_g.ne.0) then
           Nshells = 100
           call prepare_shells(Nshells)
           allocate(g_avg(Nshells))
           g_avg = 0d0
           if(taskid == master) allocate(g_avg_final(Nshells))
         endif

         t = 0.d0
         call compute_force_LJ(r,f,U,Ppot) !Initial force, energy and pressure
         call energy_kin(v,ekin,Tins) !Compute initial kinetic energy
         if (taskid==master) Ptot = rho*Tins + Ppot !AJ: modifed to use the rho in parameters.

         !Write intial results.
         if (taskid==master) then
             write(eunit,*)"#time,   K,   U,  E,  T,  Ptot"
             write(eunit,*) t, ekin, U, ekin+U, Tins, Ptot

             write(eunit_dim,*)"#time,   K,   U,  E,  T,  Ptot"
             write(eunit_dim,*) t*unit_of_time,&
                        ekin*unit_of_energy, U*unit_of_energy, (ekin+U)*unit_of_energy,&
                        Tins*epsilon, Ptot*unit_of_pressure
         endif
         
         fold=f
         do i=1,Nt !Main time loop.
            call verlet_v_step(r,v,fold,t,i,dt,U,Ppot) !Perform Verlet step.
            call andersen_therm(v,Temp) !Apply thermostat
            call energy_kin(v,ekin,Tins)
            if (taskid==master) Ptot = rho*Tins + Ppot

            !Write to file.
            if (taskid==master) then
               write(eunit,*) t, ekin, U, ekin+U, Tins, Ptot

               write(eunit_dim,*) t*unit_of_time,&
                        ekin*unit_of_energy, U*unit_of_energy, (ekin+U)*unit_of_energy,&
                        Tins*epsilon, Ptot*unit_of_pressure
            endif
            
            ! Save snapshot of g(r) to average
            if(flag_g.ne.0) then
              call rad_dist_fun(r,Nshells)
              g_avg = g_avg + g
            endif

            if(mod(i,int(0.001*Nt))==0 .and. taskid==master) then
                  write (*,"(A,F5.1,A)",advance="no") "Progress: ",i/dble(Nt)*100.,"%"
                  if (i.le.Nt) call execute_command_line('echo -e "\033[A"')
            endif
            
         enddo
         
         if(flag_g.ne.0) then
           call MPI_REDUCE(g_avg,g_avg_final,Nshells,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
           if(taskid.eq.master) then
             g_avg_final = g_avg_final/dble(Nt)
             do j=1,Nshells
               write(eunit_g,*) (j-1)*grid_shells+grid_shells/2d0, g_avg_final(j)
               write(eunit_g_dim,*) ((j-1)*grid_shells+grid_shells/2d0)*unit_of_length, g_avg_final(j)
             enddo
             deallocate(g_avg_final)
           endif
           deallocate(g_avg)
           call deallocate_g_variables()
         endif

         return
      end subroutine vvel_solver


end module integraforces
