      program main
      use parameters
      use init
      use pbc
      use integraforces
      use statvis
      use rad_dist

      implicit none
      character(len=50) :: input_name
      real*8, allocatable :: pos(:,:), vel(:,:), fold(:,:)
      real*8, allocatable :: epotVEC(:), PVEC(:), ekinVEC(:), etotVEC(:), TinsVEC(:)
      real*8, allocatable :: epotVECins(:), g_avg(:), g_squared_avg(:)
      real*8, allocatable :: Xpos(:), Ypos(:), Zpos(:), posAUX(:,:), noPBC(:,:)
      
      real*8  :: time,ekin,epot,Tins,P,etot
      real*8  :: epotAUX,epotMEAN,PMEAN,epotVAR,PVAR
      real*8  :: ekinMEAN,ekinVAR,etotMEAN,etotVAR,TinsMEAN,TinsVAR
      real*8  :: Xmean,Ymean,Zmean,Xvar,Yvar,Zvar
      integer :: i,j,flag_g,k,cnt

      integer :: Nshells
      real*8  :: ti,tf,elapsed_time !AJ: collective timing of program.
      
      
      call cpu_time(ti)

      ! Per executar el programa cal fer >> main.x input_file. Si no, donara error.
      if (command_argument_count() == 0) stop "ERROR: Cridar fent >> ./main.x input_path"
      call get_command_argument(1, input_name)
      
      ! Obrim input i el llegim
      open(unit=10, file=input_name)
      call get_param(10)
      close(10)

      print*,"------------------------Parameters-------------------------------"
      print"(A,X,I5,2X,A,X,I1)", "N=",N,"D=",D
      print"(A,X,E14.7)","dt_sim=",dt_sim
      print"(A,X,I8)","seed=",seed
      print"(A,X,F4.2,2X,A,X,F6.2)","rho=",rho,"T=",T_ref
      print"(A,X,F10.5,2X,A,X,F10.5,2X,A,X,F10.5)","eps=",epsilon,"sigma=",sigma,"rc=",rc
      print"(A,X,I10)","n_equil=",n_equil
      print"(A,X,I3,2X,I10,2X,I10)","n_meas,n_conf,n_total=",n_meas,n_conf,n_total
      print*,"-----------------------------------------------------------------"
      
      ! Allocates
      allocate(pos(D,N))
      allocate(vel(D,N))  
      allocate(fold(D,N)) 
      allocate(epotVECins(n_meas))
      allocate(PVEC(n_conf))
      allocate(ekinVEC(n_conf))
      allocate(epotVEC(n_conf))
      allocate(etotVEC(n_conf))
      allocate(TinsVEC(n_conf))
      allocate(Xpos(N))
      allocate(Ypos(N))
      allocate(Zpos(N))
      allocate(posAUX(D,N))
      allocate(noPBC(D,N))

      ! Initialize positions and velocities
      call init_sc(pos)
      call init_vel(vel, 10.d0)  ! Cridem amb temperatura reduida (T'=10) molt alta per fer el melting


      ! Start melting of the system
      open(unit=10,file="results/thermodynamics_initialization.dat") ! AJ: open result for initialitzation.
      open(unit=12,file="results/dimensionalized/thermodynamics_initialization_dim.dat") 
      open(unit=11,file="results/init_conf.xyz")

      call writeXyz(D,N,pos,11)

      flag_g = 0 ! DM: don't write g(r)
      print*,"------Melting Start------"
      call vvel_solver(5000,1.d-4,pos,vel,1000.d0,10,12,0,0,flag_g) ! AJ: Initialization of system.
      print*,"------Melting Completed------"

      call writeXyz(D,N,pos,11) ! AJ: write initial configuration, check that it is random.

      close(10)
      close(11)
      close(12)

      ! Start dynamics
      ! Perform equilibration of the system
      call init_vel(vel, T_ref) ! AJ: reescale to target temperature

      open(unit=10,file="results/thermodynamics_equilibration.dat")
      open(unit=11,file="results/dimensionalized/thermodynamics_equilibration_dim.dat")

      print*,"------Equilibration Start------"
      call vvel_solver(n_equil,dt_sim,pos,vel,T_ref,10,11,0,0,flag_g)
      print*,"------Equilibration Completed------"

      close(10)
      close(11)

      ! Once the system is equilibrated, start dynamics of the system
      open(unit=10,file="results/thermodynamics.dat")
      open(unit=11,file="results/trajectory.xyz")
      open(unit=12,file="results/radial_distribution.dat")
      open(unit=13,file="results/mean_epot.dat")
      open(unit=14,file="results/diffcoeff.dat")
      open(unit=15,file="results/averages.dat")
      open(unit=21,file="results/ekinBIN.dat")
      open(unit=22,file="results/epotBIN.dat")
      open(unit=23,file="results/correlation_energy.dat")

      open(unit=16,file="results/dimensionalized/thermodynamics_dim.dat")
      open(unit=17,file="results/dimensionalized/trajectory_dim.xyz")
      open(unit=18,file="results/dimensionalized/mean_epot_dim.dat")
      open(unit=19,file="results/dimensionalized/diffcoeff_dim.dat")
      open(unit=24,file="results/dimensionalized/averages_dim.dat")
      open(unit=25,file="results/dimensionalized/radial_distribution_dim.dat")

      write(10,*)"#t,   K,   U,  E,  T,  v_tot,  Ptot"
      write(16,*)"#t,   K,   U,  E,  T,  Ptot"
      
      ! Set the variables for computing g(r) (defined in rad_dist module)
      Nshells = 100
      call prepare_shells(Nshells)
      ! allocate averaging vectors for g
      allocate(g_avg(Nshells))
      allocate(g_squared_avg(Nshells))
      g_avg = 0d0
      g_squared_avg = 0d0
      

      k = 0
      cnt = 0
      epotAUX = 0.d0
      noPBC = 0.d0

      print*,"------Simulation Start------"

      call compute_force_LJ(pos,fold,epot,P)
      do i = 1,n_total
     
          posAUX = pos

            call verlet_v_step(pos,vel,fold,time,i,dt_sim,epot,P)
            call andersen_therm(vel,T_ref)

          do j=1,N
              ! X
                if ((pos(1,j)-posAUX(1,j)).gt.(0.9d0*L)) then
                    noPBC(1,j) = noPBC(1,j) - L
              elseif ((pos(1,j)-posAUX(1,j)).lt.(0.9d0*L)) then
                    noPBC(1,j) = noPBC(1,j) + L
              endif
              ! Y
                if ((pos(2,j)-posAUX(2,j)).gt.(0.9d0*L)) then
                    noPBC(2,j) = noPBC(2,j) - L
              elseif ((pos(2,j)-posAUX(2,j)).lt.(0.9d0*L)) then
                    noPBC(2,j) = noPBC(2,j) + L
              endif
              ! Z
                if ((pos(3,j)-posAUX(3,j)).gt.(0.9d0*L)) then
                    noPBC(3,j) = noPBC(3,j) - L
              elseif ((pos(3,j)-posAUX(3,j)).lt.(0.9d0*L)) then
                    noPBC(3,j) = noPBC(3,j) + L
              endif
            enddo
            ! Càlcul del coeficient de difusió per cada dimensió (s'han evitat les PBC)
            Xpos(:) = pos(1,:) + noPBC(1,:)
            Ypos(:) = pos(2,:) + noPBC(2,:)
            Zpos(:) = pos(3,:) + noPBC(3,:)
            call estad(N,Xpos,Xmean,Xvar)
            call estad(N,Ypos,Ymean,Yvar)
            call estad(N,Zpos,Zmean,Zvar)
            write(14,*) 2.d0*dble(i), Xvar*dble(N), Yvar*dble(N), Zvar*dble(N)
            write(19,*) 2.d0*time, Xvar*dble(N)*unit_of_length**2,&
                  Yvar*dble(N)*unit_of_length**2,&
                  Zvar*dble(N)*unit_of_length**2

            k = k+1
            epotVECins(k) = epot

            if(mod(i,n_meas) == 0) then ! AJ : measure every n_meas steps
                  ! Average de epot cada n_meas. Ho escribim en un fitxer
                  k = 0
                  cnt = cnt+1
                  call estad(n_meas,epotVECins,epotMEAN,epotVAR)
                  epotAUX = (epotMEAN+epotAUX*dble(cnt-1))/dble(cnt)
                  write(13,*) i, epotAUX
                  write(18,*) i, epotAUX*unit_of_energy

                  call energy_kin(vel,ekin,Tins)
                  write(10,*) time, ekin, epot, ekin+epot, Tins, dsqrt(sum(sum(vel,2)**2)), P+rho*Tins
                  write(16,*) time*unit_of_time,&
                        ekin*unit_of_energy, epot*unit_of_energy, (ekin+epot)*unit_of_energy,&
                        Tins*epsilon, (P+rho*Tins)*unit_of_pressure
                  ekinVEC(cnt) = ekin
                  epotVEC(cnt) = epot
                  etotVEC(cnt) = ekin+epot
                  TinsVEC(cnt) = Tins
                  PVEC(cnt) = P+rho*Tins
                  call writeXyz(D,N,pos,11)
                  call writeXyz(D,N,pos*unit_of_length,17)

                  ! Compute g(r) and write to file
                  call rad_distr_fun(pos,Nshells)
                  g_avg = g_avg + g
                  g_squared_avg = g_squared_avg + g**2
            endif

            if(mod(i,int(0.001*n_total))==0) then
                   write (*,"(A,F5.1,A)",advance="no") "Progress: ",i/dble(n_total)*100.,"%"
                   if (i<n_total) call execute_command_line('echo "\033[A"')
            endif
      enddo
      write (*,*)
      write (*,*) "----Simulation Completed----"

      close(10)
      close(11)
      close(13)
      close(14)

      close(16)
      close(17)
      close(18)
      close(19)
      
      ! Average de la g(r)
      g_avg = g_avg/dble(n_conf)
      g_squared_avg = g_squared_avg/dble(n_conf)
      write(12,*) " # r (reduced units),   g(r),   std dev "
      write(25,*) " # r (Angstroms),   g(r),   std dev "
      do i=1,Nshells
        write(12,*) grid_shells*(i-1)+grid_shells/2d0, dsqrt(g_squared_avg(i) - g_avg(i)**2)
        write(25,*) (grid_shells*(i-1)+grid_shells/2d0)*sigma, g_avg(i), dsqrt(g_squared_avg(i) - g_avg(i)**2)
      enddo
      close(12)
      close(25)
        

      ! Averages finals
      deallocate(epotVECins)

      call estad(n_conf,ekinVEC,ekinMEAN,ekinVAR)
      call estad(n_conf,epotVEC,epotMEAN,epotVAR)
      call estad(n_conf,etotVEC,etotMEAN,etotVAR)
      call estad(n_conf,TinsVEC,TinsMEAN,TinsVAR)
      call estad(n_conf,PVEC,PMEAN,PVAR)
      write(15,*) "Sample mean and Statistical error"
      write(15,*) "Kinetic Energy", ekinMEAN, dsqrt(ekinVAR)
      write(15,*) "Potential Energy", epotMEAN, dsqrt(epotVAR)
      write(15,*) "Total Energy", etotMEAN, dsqrt(etotVAR)
      write(15,*) "Instant Temperature", TinsMEAN, dsqrt(TinsVAR)
      write(15,*) "Pressure", PMEAN, dsqrt(PVAR)
      close(15)

      write(24,*) "Sample mean and Statistical error"
      write(24,*) "Kinetic Energy", ekinMEAN*unit_of_energy, dsqrt(ekinVAR)*unit_of_energy
      write(24,*) "Potential Energy", epotMEAN*unit_of_energy, dsqrt(epotVAR)*unit_of_energy
      write(24,*) "Total Energy", etotMEAN*unit_of_energy, dsqrt(etotVAR)*unit_of_energy
      write(24,*) "Instant Temperature", TinsMEAN*epsilon, dsqrt(TinsVAR)*epsilon
      write(24,*) "Pressure", PMEAN*unit_of_pressure, dsqrt(PVAR)*unit_of_pressure
      close(24)

      ! Binning de les energies cinètica i potencial
      call binning(n_conf,ekinVEC,50,21)
      call binning(n_conf,epotVEC,50,22)
      
      ! Funció d'autocorrelació per l'energia total
      call corrtime(n_conf,etotVEC,23)

      close(21)
      close(22)
      close(23)

      ! Deallocates
      deallocate(pos) 
      deallocate(vel)
      deallocate(fold)
      deallocate(epotVEC)
      deallocate(PVEC)
      deallocate(ekinVEC)
      deallocate(etotVEC)
      deallocate(TinsVEC)
      deallocate(g_avg)
      deallocate(g_squared_avg)
      deallocate(Xpos)
      deallocate(Ypos)
      deallocate(Zpos)
      deallocate(posAUX)
      deallocate(noPBC)
      call deallocate_g_variables()

      call cpu_time(tf)
      elapsed_time = tf-ti
      print"(A,X,F14.7,X,A)","End program, time elapsed:",elapsed_time,"seconds"

      end program main
