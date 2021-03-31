      program main
      use parameters
      use init
      use pbc
      use integraforces
      use statvis
      use rad_dist
  
      implicit none
      include 'mpif.h'
      character(len=50)   :: input_name
      real*8, allocatable :: pos(:,:), vel(:,:), fold(:,:), g_avg(:), g_squared_avg(:)
      real*8, allocatable :: epotVECins(:), epotVEC(:), PVEC(:), ekinVEC(:), etotVEC(:), TinsVEC(:)
      real*8, allocatable :: Xpos(:), Ypos(:), Zpos(:)
      
      real*8              :: time,ekin,epot,Tins,P,etot
      real*8              :: epotAUX,epotMEAN,PMEAN,epotVAR,PVAR
      real*8              :: ekinMEAN,ekinVAR,etotMEAN,etotVAR,TinsMEAN,TinsVAR
      real*8              :: Xmean,Ymean,Zmean,Xvar,Yvar,Zvar
      integer             :: i,ierror,Nshells,flag_g,k,cnt
      real*8              :: ti_global,tf_global,elapsed_time !AJ: collective timing of program.
      
      ! Init MPI
      call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
      ti_global = MPI_WTIME()

      if (taskid == master) then
          allocate(aux_pos(numproc))
          allocate(aux_size(numproc))
      end if
  
      ! To execute the program >> main.x input_file. Otherwise, an error will occur.
      if (command_argument_count() == 0) stop "ERROR: call using >> ./main.x input_path"
      call get_command_argument(1, input_name)
      
      ! Open and read input 
      open(unit=10, file=input_name)
      call get_param(10)
      close(10)
  
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
      
      if(taskid==master) then
            print*,"------------------------Parameters-------------------------------"
            print"(A,X,I5,2X,A,X,I1)", "N=",N,"D=",D
            print"(A,X,E14.7)","dt_sim=",dt_sim
            print"(A,X,I8)","seed=",seed
            print"(A,X,F4.2,2X,A,X,F5.2)","rho=",rho,"T=",T_ref
            print"(A,X,F7.4,2X,A,X,F4.2)","eps=",epsilon,"sigma=",sigma,"rc=",rc
            print"(A,X,I3,2X,I5,2X,I10)","n_meas,n_conf,n_total=",n_meas,n_conf,n_total
            print*,"-----------------------------------------------------------------"
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
  
      ! Initialize positions and velocities
      call init_sc_gather(pos)
      call init_vel_gather(vel, 1000.d0)
  
      !Start melting
      if(taskid==master) then
            open(unit=10,file="results/thermodynamics_initialization.dat")
            open(unit=11,file="results/init_conf.xyz")
            call writeXyz(D,N,pos,11)
      end if
     
      flag_g = 0
      if(taskid==master)print*,"------Melting Start------"
      call vvel_solver(5000,1.d-4,pos,vel,1000.d0,10,0,flag_g)
      if(taskid==master)print*,"------Melting Completed------"
      !End melting

      ! Start dynamics
      ! Perform equilibration of the system
      call init_vel_gather(vel, T_ref) ! Reescale to target temperature

      if(taskid==master) then
            close(10)
            open(unit=10,file="results/thermodynamics_equilibration.dat")
      end if

      if(taskid==master)print*,"------Equilibration Start------"
      call vvel_solver(n_equil,dt_sim,pos,vel,T_ref,10,0,flag_g)
      if(taskid==master)print*,"------Equilibration Completed------"

   
      !Prepare files for main simulation   
      if(taskid==master) then
            close(10)
            close(11)
            open(unit=10,file="results/thermodynamics.dat")
            open(unit=11,file="results/trajectory.xyz")
            open(unit=12,file="results/radial_distribution.dat")
            open(unit=13,file="results/mean_epot.dat")
            open(unit=14,file="results/diffcoeff.dat")
            open(unit=15,file="results/averages.dat")
      
            write(10,*)"#t,   K,   U,  E,  T,  v_tot,  Ptot"
      end if
  
      !Prepare g(r) variables
      Nshells = 100
      call prepare_shells(Nshells)
      if(taskid==master) then
            allocate(g_avg(Nshells))
            allocate(g_squared_avg(Nshells))
            g_avg = 0d0
            g_squared_avg = 0d0
            
            print*,"------Simulation Start------"
      endif
      
      
      k = 0
      cnt = 0
      epotAUX = 0.d0
      
      call compute_force_LJ(pos,fold,epot,P)
      do i = 1,n_total
      
            call verlet_v_step(pos,vel,fold,time,i,dt_sim,epot,P)
            call andersen_therm(vel,T_ref)
            
            ! Càlcul del coeficient de difusió per cada dimensió
            Xpos(:) = pos(1,:)
            Ypos(:) = pos(2,:)
            Zpos(:) = pos(3,:)
            call estad(N,Xpos,Xmean,Xvar)
            call estad(N,Ypos,Ymean,Yvar)
            call estad(N,Zpos,Zmean,Zvar)
            if (taskid.eq.master) then
                write(14,*) 2.d0*time, Xvar*dble(N), Yvar*dble(N), Zvar*dble(N)
            endif

            k = k+1
            epotVECins(k) = epot
      
            if(mod(i,n_meas) == 0) then ! AJ : measure every n_meas steps
                  ! Average de epot cada n_meas. Ho escribim en un fitxer
                  k = 0
                  cnt = cnt+1
                  call estad(n_meas,epotVECins,epotMEAN,epotVAR)
                  if (taskid.eq.master) then
                        epotAUX = (epotMEAN+epotAUX*dble(cnt-1))/dble(cnt)
                        write(13,*) i, epotAUX
                  endif
                  
                  call energy_kin(vel,ekin,Tins)
                  if(taskid==master) then
                        write(10,*) time, ekin, epot, ekin+epot, Tins, dsqrt(sum(sum(vel,2)**2)), P+rho*Tins
                        ekinVEC(cnt) = ekin
                        epotVEC(cnt) = epot
                        etotVEC(cnt) = ekin+epot
                        TinsVEC(cnt) = Tins
                        PVEC(cnt) = P+rho*Tins
                        call writeXyz(D,N,pos,11)
                  end if
                  call rad_dist_fun_pairs(pos,Nshells)
                  if(taskid==master) then
                        g_avg = g_avg + g
                        g_squared_avg = g_squared_avg + g**2
                  endif
            endif
      
            if(mod(i,int(0.001*n_total))==0 .and. taskid==master) then
                  write (*,"(A,F5.1,A)",advance="no") "Progress: ",i/dble(n_total)*100.,"%"
                  if (i<n_total) call execute_command_line('echo "\033[A"')
            endif
      enddo
      
      if(taskid==master) then
            write (*,*)
            write (*,*) "----Simulation Completed----"
            close(10)
            close(11)
            close(13)
            close(14)
      end if
      
      ! Average de la g(r)
      if(taskid==master) then
            g_avg = g_avg/dble(n_conf)
            g_squared_avg = g_squared_avg/dble(n_conf)
            write(12,*) " # r (reduced units), r (Angstroms),   g(r),   std dev "
            do i=1,Nshells
                  write(12,*) grid_shells*(i-1), grid_shells*(i-1)*sigma, g_avg(i), dsqrt(g_squared_avg(i) - g_avg(i)**2)
            enddo
            close(12)
      endif

      ! Averages finals
      if (allocated(epotVECins)) deallocate(epotVECins)

      call estad(n_conf,ekinVEC,ekinMEAN,ekinVAR)
      call estad(n_conf,epotVEC,epotMEAN,epotVAR)
      call estad(n_conf,etotVEC,etotMEAN,etotVAR)
      call estad(n_conf,TinsVEC,TinsMEAN,TinsVAR)
      call estad(n_conf,PVEC,PMEAN,PVAR)
      
      if (taskid.eq.master) then
          write(15,*) "Sample mean and Variance"
          write(15,*) "Kinetic Energy", ekinMEAN, ekinVAR
          write(15,*) "Potential Energy", epotMEAN, epotVAR
          write(15,*) "Total Energy", etotMEAN, etotVAR
          write(15,*) "Instant Temperature", TinsMEAN, TinsVAR
          write(15,*) "Pressure", PMEAN, PVAR
          close(15)
      endif

      ! Binning de les energies cinètica i potencial
      call binning(n_conf,ekinVEC,50,"results/ekinBIN.dat")
      call binning(n_conf,epotVEC,50,"results/epotBIN.dat")
      
      ! Funció d'autocorrelació per l'energia total
      call corrtime(n_conf,etotVEC,"results/correlation_energy.dat")
  
  
      if (allocated(pos)) deallocate(pos)
      if (allocated(vel)) deallocate(vel)
      if (allocated(fold)) deallocate (fold)
      if (allocated(aux_pos)) deallocate(aux_pos)
      if (allocated(aux_size)) deallocate(aux_size)
      if (allocated(epotVEC)) deallocate(epotVEC)
      if (allocated(PVEC)) deallocate(PVEC)
      if (allocated(ekinVEC)) deallocate(ekinVEC)
      if (allocated(etotVEC)) deallocate(etotVEC)
      if (allocated(TinsVEC)) deallocate(TinsVEC)
      if (allocated(Xpos)) deallocate(Xpos)
      if (allocated(Ypos)) deallocate(Ypos)
      if (allocated(Zpos)) deallocate(Zpos)
      if (allocated(g_avg)) deallocate(g_avg)
      if (allocated(g_squared_avg)) deallocate(g_squared_avg)
      call deallocate_g_variables()
  
      tf_global = MPI_WTIME()
      call MPI_REDUCE(tf_global-ti_global,elapsed_time,1,MPI_DOUBLE_PRECISION,MPI_MAX,master,MPI_COMM_WORLD,ierror)
      if(taskid==master) then
            print"(A,X,F14.7,X,A)","End program, time elapsed:",elapsed_time,"seconds"
      end if
  
      call MPI_FINALIZE(ierror)
      end program main
  
