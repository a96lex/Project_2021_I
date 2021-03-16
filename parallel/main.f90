program main
    use parameters
    use init
    ! use pbc
    use integraforces
    use statvis
    ! use rad_dist

    implicit none
    include 'mpif.h'
    character(len=50)   :: input_name
    real*8, allocatable :: pos(:,:), vel(:,:)

    integer :: ierror
    
    ! Init MPI
    call MPI_INIT(ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


    ! Per executar el programa cal fer >> main.x input_file. Si no, donara error.
    if (command_argument_count() == 0) stop "ERROR: Cridar fent >> ./main.x input_path"
    call get_command_argument(1, input_name)
    
    ! Obrim input i el llegim
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
        print"(A,X,F4.2,2X,A,X,F5.2)","rho=",rho,"T=",T_ref
        print"(A,X,F7.4,2X,A,X,F4.2)","eps=",epsilon,"sigma=",sigma,"rc=",rc
        print"(A,X,I3,2X,I5,2X,I5)","n_meas,n_conf,n_total=",n_meas,n_conf,n_total
        print*,"-----------------------------------------------------------------"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Initialize positions and velocities
    ! call init_sc_paralel(pos)
    if (taskid == master) call init_sc(pos)

    if(taskid==master) open(1,file="results/init_conf.xyz")

    if(taskid==master) then
        call writeXyz(D,N,pos,1)
    end if

    if (allocated(pos)) deallocate(pos)
    if (allocated(vel)) deallocate(vel)

    call MPI_FINALIZE(ierror)
end program main
