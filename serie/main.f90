      program main
      use parameters
      use init
      use pbc
      use integraforces
      use statvis
      use rad_dist

      implicit none
      character(len=50) :: input_name
      real*8, allocatable :: pos(:,:), vel(:,:)
      real*8, allocatable :: epotVEC(:), PVEC(:), ekinVEC(:), etotVEC(:), TinsVEC(:)
      real*8, allocatable :: epotVECins(:), g_avg(:), g_squared_avg(:)
      real*8, allocatable :: Xpos(:), Ypos(:), Zpos(:)
      
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
      print"(A,X,F4.2,2X,A,X,F5.2)","rho=",rho,"T=",T_ref
      print"(A,X,F7.4,2X,A,X,F4.2)","eps=",epsilon,"sigma=",sigma,"rc=",rc
      print"(A,X,I3,2X,I5,2X,I5)","n_meas,n_conf,n_total=",n_meas,n_conf,n_total
      print*,"-----------------------------------------------------------------"
      
      ! Allocates
      allocate(pos(D,N))
      allocate(vel(D,N))   
      allocate(epotVECins(n_meas))
      allocate(PVEC(n_conf))
      allocate(ekinVEC(n_conf))
      allocate(epotVEC(n_conf))
      allocate(etotVEC(n_conf))
      allocate(TinsVEC(n_conf))
      allocate(Xpos(N))
      allocate(Ypos(N))
      allocate(Zpos(N))

      ! Initialize positions and velocities
      call init_sc(pos)
      call init_vel(vel, 10.d0)  ! Cridem amb temperatura reduida (T'=10) molt alta per fer el melting


      ! Start melting of the system
      open(unit=10,file="results/thermodynamics_initialization.dat") ! AJ: open result for initialitzation.
      open(unit=11,file="results/init_conf.xyz")

      call writeXyz(D,N,pos,11)

      flag_g = 0 ! DM: don't write g(r)
      call vvel_solver(5000,1.d-4,pos,vel,1000.d0,10,0,flag_g) ! AJ: Initialization of system.

      call writeXyz(D,N,pos,11) ! AJ: write initial configuration, check that it is random.

      close(10)
      close(11)

      ! Start dynamics
      ! Perform equilibration of the system
      call init_vel(vel, T_ref) ! AJ: reescale to target temperature

      open(unit=10,file="results/thermodynamics_equilibration.dat")

      call vvel_solver(n_equil,dt_sim,pos,vel,T_ref,10,0,flag_g)

      close(10)

      ! Once the system is equilibrated, start dynamics of the system
      open(unit=10,file="results/thermodynamics.dat")
      open(unit=11,file="results/trajectory.xyz")
      open(unit=12,file="results/radial_distribution.dat")
      open(unit=13,file="results/mean_epot.dat")
      open(unit=14,file="results/diffcoeff.dat")
      open(unit=15,file="results/averages.dat")
      
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
      print*,"------Simulation Start------"
      do i = 1,n_total
      
            call verlet_v_step(pos,vel,time,i,dt_sim,epot,P)
            call andersen_therm(vel,dt_sim,T_ref)

            ! Càlcul del coeficient de difusió per cada dimensió
            Xpos(:) = pos(1,:)
            Ypos(:) = pos(2,:)
            Zpos(:) = pos(3,:)
            call estad(N,Xpos,Xmean,Xvar)
            call estad(N,Ypos,Ymean,Yvar)
            call estad(N,Zpos,Zmean,Zvar)
            write(14,*) 2.d0*dble(i), Xvar*dble(N), Yvar*dble(N), Zvar*dble(N)

            k = k+1
            epotVECins(k) = epot

            if(mod(i,n_meas) == 0) then ! AJ : measure every n_meas steps
                  ! Average de epot cada n_meas. Ho escribim en un fitxer
                  k = 0
                  cnt = cnt+1
                  call estad(n_meas,epotVECins,epotMEAN,epotVAR)
                  write(13,*) i, (epotMEAN+epotAUX*dble(cnt-1))/dble(cnt)
                  epotAUX = (epotMEAN+epotAUX*dble(cnt-1))/dble(cnt)

                  call energy_kin(vel,ekin,Tins)
                  write(10,*) time, ekin, epot, ekin+epot, Tins, dsqrt(sum(sum(vel,2)**2)), P+rho*Tins
                  ekinVEC(cnt) = ekin
                  epotVEC(cnt) = epot
                  etotVEC(cnt) = ekin+epot
                  TinsVEC(cnt) = Tins
                  PVEC(cnt) = P+rho*Tins
                  call writeXyz(D,N,pos,11)

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
      
      ! Average de la g(r)
      g_avg = g_avg/dble(n_conf)
      g_squared_avg = g_squared_avg/dble(n_conf)
      write(12,*) " # r (reduced units), r (Angstroms),   g(r),   std dev "
      do i=1,Nshells
        write(12,*) grid_shells*(i-1), grid_shells*(i-1)*sigma, g_avg(i), dsqrt(g_squared_avg(i) - g_avg(i)**2)
      enddo
      close(12)
        

      ! Averages finals
      deallocate(epotVECins)

      call estad(n_conf,ekinVEC,ekinMEAN,ekinVAR)
      call estad(n_conf,epotVEC,epotMEAN,epotVAR)
      call estad(n_conf,etotVEC,etotMEAN,etotVAR)
      call estad(n_conf,TinsVEC,TinsMEAN,TinsVAR)
      call estad(n_conf,PVEC,PMEAN,PVAR)
      write(15,*) "Sample mean and Variance"
      write(15,*) "Kinetic Energy", ekinMEAN, ekinVAR
      write(15,*) "Potential Energy", epotMEAN, epotVAR
      write(15,*) "Total Energy", etotMEAN, etotVAR
      write(15,*) "Instant Temperature", TinsMEAN, TinsVAR
      write(15,*) "Pressure", PMEAN, PVAR

      ! Binning de les energies cinètica i potencial
      call binning(n_conf,ekinVEC,50,"results/ekinBIN.dat")
      call binning(n_conf,epotVEC,50,"results/epotBIN.dat")
      
      ! Funció d'autocorrelació per l'energia total
      call corrtime(n_conf,etotVEC,"results/correlation_energy.dat")

      close(15)

      ! Deallocates
      deallocate(pos) 
      deallocate(vel)
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
      call deallocate_g_variables()

      call cpu_time(tf)
      elapsed_time = tf-ti
      print"(A,X,F14.7,X,A)","End program, time elapsed:",elapsed_time,"seconds"

      end program main
