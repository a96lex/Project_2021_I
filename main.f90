      program main
      use parameters
      use init
      use pbc
      use integraforces
      use statvis
      use radial_distribution

      implicit none
      character(len=50) :: input_name
      real*8, allocatable :: pos(:,:), vel(:,:)
      real*8, allocatable :: epotVEC(:), PVEC(:), ekinVEC(:), etotVEC(:), TinsVEC(:)
      real*8, allocatable :: g_avg(:), g_squared_avg(:)
      
      real*8 :: time,ekin,epot,Tins,P,etot,mes
      real*8 :: epotAUX,PAUX,epotMEAN,PMEAN,epotVAR,PVAR
      real*8 :: ekinMEAN,ekinVAR,etotMEAN,etotVAR,TinsMEAN,TinsVAR
      integer :: i,j,flag_g,k,cnt

      integer :: Nshells
      

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
      print"(A,X,F4.2,2X,A,X,F4.2)","rho=",rho,"T=",T_ref
      print"(A,X,F7.4,2X,A,X,F4.2)","eps=",epsilon,"sigma=",sigma,"rc=",rc
      print"(A,X,I3,2X,I5,2X,I5)","n_meas,n_conf,n_total=",n_meas,n_conf,n_total
      print*,"-----------------------------------------------------------------"
      
      ! Allocates
      allocate(pos(D,N))
      allocate(vel(D,N))   
      allocate(epotVEC(n_meas))
      allocate(PVEC(n_meas))
      allocate(ekinVEC(n_total/n_meas))
      allocate(etotVEC(n_total/n_meas))
      allocate(TinsVEC(n_total/n_meas))

      ! Initialize positions and velocities
      call init_sc(pos)
      call init_vel(vel, 10.d0)  ! Cridem amb temperatura reduida (T'=10) molt alta per fer el melting


      ! Start melting of the system
      open(unit=10,file="results/thermodynamics_initialization.dat") ! AJ: open result for initialitzation.
      open(unit=11,file="results/init_conf.xyz")

      call writeXyz(D,N,pos,11)

      flag_g = 0 ! DM: don't write g(r)
      call vvel_solver(5000,1.d-4,pos,vel,10.d0,10,0,flag_g) ! AJ: Initialization of system.

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
      open(unit=14,file="results/mean_press.dat")
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
      PAUX = 0.d0
      print*,"------Simulation Start------"
      do i = 1,n_total
      
            call verlet_v_step(pos,vel,time,dt_sim,epot,P)
            call andersen_therm(vel,dt_sim,T_ref)

            k = k+1
            epotVEC(k) = epot
            PVEC(k) = P

            if(mod(i,n_meas) == 0) then ! AJ : measure every n_meas steps
                  ! Average de epot i P cada n_meas. Ho escribim en un fitxer cada un
                  k = 0
                  cnt = cnt+1
                  call estad(n_meas,epotVEC,epotMEAN,epotVAR)
                  call estad(n_meas,PVEC,PMEAN,PVAR)
                  write(13,*) i, (epotMEAN+epotAUX*dble(cnt-1))/dble(cnt)
                  write(14,*) i, (PMEAN+PAUX*dble(cnt-1))/dble(cnt)
                  epotAUX = (epotMEAN+epotAUX*dble(cnt-1))/dble(cnt)
                  PAUX = (PMEAN+PAUX*dble(cnt-1))/dble(cnt)

                  call energy_kin(vel,ekin,Tins)

                  write(10,*) time, ekin, epot, ekin+epot, Tins, sum(vel,2), P
                  !We sould also need pressure here.
                  call writeXyz(D,N,pos,11)
                  ! Compute g(r) and write to file
                  call rad_distr_fun(pos,Nshells)
                  g_avg = g_avg + g
                  g_squared_avg = g_squared_avg + g**2
            endif

            if(mod(i,int(0.1*n_total))==0) then
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
      g_avg = g_avg/(n_total/n_meas)
      g_squared_avg = g_squared_avg/(n_total/n_meas)
      write(12,*) " # r (reduced units), r (Angstroms),   g(r),   std dev "
      do i=1,Nshells
        write(12,*) grid_shells*(i-1), grid_shells*(i-1)*sigma, g_avg(i), dsqrt(g_squared_avg(i) - g_avg(i)**2)
      enddo
      close(12)
        

      ! Averages finals (faltarà la pressió)
      deallocate(epotVEC)
      deallocate(PVEC)
      allocate(epotVEC(n_total/n_meas))
      allocate(PVEC(n_total/n_meas))
      open(10,file="results/thermodynamics.dat",status="old")
      do i=1,n_total/n_meas
        read(10,*) time, ekin, epot, etot, Tins, mes, P
        ekinVEC(i) = ekin
        epotVEC(i) = epot
        etotVEC(i) = etot
        TinsVEC(i) = Tins
        PVEC(i) = P
      enddo

      call estad(n_total/n_meas,ekinVEC,ekinMEAN,ekinVAR)
      call estad(n_total/n_meas,epotVEC,epotMEAN,epotVAR)
      call estad(n_total/n_meas,etotVEC,etotMEAN,etotVAR)
      call estad(n_total/n_meas,TinsVEC,TinsMEAN,TinsVAR)
      call estad(n_total/n_meas,PVEC,PMEAN,PVAR)
      write(15,*) "Sample mean and Variance"
      write(15,*) "Kinetic Energy", ekinMEAN, ekinVAR
      write(15,*) "Potential Energy", epotMEAN, epotVAR
      write(15,*) "Total Energy", etotMEAN, etotVAR
      write(15,*) "Instant Temperature", TinsMEAN, TinsVAR
      write(15,*) "Pressure", PMEAN, PVAR

      ! Binning de les energies cinètica i potencial
      call binning(n_total/n_meas,ekinVEC,50,"results/ekinBIN.dat")
      call binning(n_total/n_meas,epotVEC,50,"results/epotBIN.dat")

      close(10)
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
      call deallocate_g_variables()

      end program main
