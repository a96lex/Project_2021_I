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
    real*8, allocatable :: pos(:,:), vel(:,:)
    real*8              :: time,epot,P
    integer             :: i,ierror,Nshells
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
    
    if(taskid==master) then
        print*,"------------------------Parameters-------------------------------"
        print"(A,X,I5,2X,A,X,I1)", "N=",N,"D=",D
        print"(A,X,E14.7)","dt_sim=",dt_sim
        print"(A,X,I8)","seed=",seed
        print"(A,X,F4.2,2X,A,X,F5.2)","rho=",rho,"T=",T_ref
        print"(A,X,F7.4,2X,A,X,F4.2)","eps=",epsilon,"sigma=",sigma,"rc=",rc
        print"(A,X,I3,2X,I5,2X,I5)","n_meas,n_conf,n_total=",n_meas,n_conf,n_total
        print*,"-----------------------------------------------------------------"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierror) 

    ! Initialize positions and velocities
    call init_sc_gather(pos)
    call init_vel_gather(vel, 1000.d0)

    if(taskid==master) then
        open(10,file="results/init_conf.xyz")
        call writeXyz(D,N,pos,10)
        close(10)
    end if
   
   !Test v_verlet
    do i=1,1000
      call verlet_v_step(pos,vel,time,i,dt_sim,epot,P)
    enddo
 
    ! Start g(r) test: david: torno a cridar init_sc perque amb 2 processadors el resultat de verlet em dona problemes
    call init_sc_gather(pos)
    Nshells = 100
    call prepare_shells(Nshells)
    call rad_dist_fun(pos,Nshells)
    
    ! BARRIER FOR TESTING
    !call sleep(floor(2d0))
    !call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
    call divide_particles_pairs()
    call rad_dist_fun_pairs(pos,Nshells)
    if(taskid == master) then
        open(11, file="results/radial_distribution.dat")
        open(111, file="results/radial_distribution_pairs.dat")
        do i=1,Nshells
            write(11,*) (i-1)*grid_shells+grid_shells/2d0,g(i)
            write(111,*) (i-1)*grid_shells+grid_shells/2d0,g_p(i)
        enddo
    endif
    call deallocate_g_variables()
    !End g(r) test

    if (allocated(pos)) deallocate(pos)
    if (allocated(vel)) deallocate(vel)
    if (allocated(aux_pos)) deallocate(aux_pos)
    if (allocated(aux_size)) deallocate(aux_size)

    tf_global = MPI_WTIME()
    call MPI_REDUCE(tf_global-ti_global,elapsed_time,1,MPI_DOUBLE_PRECISION,MPI_MAX,master,MPI_COMM_WORLD,ierror)
    if(taskid==master) then
        print"(A,X,F14.7,X,A)","End program, time elapsed:",elapsed_time,"seconds"
    end if

    call MPI_FINALIZE(ierror)
end program main
