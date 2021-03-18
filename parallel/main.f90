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

    !For force tests, remove for non-testing version
    real*8,allocatable  :: f(:,:)
    real*8              :: U,P

    integer :: i,ierror,Nshells
    
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
    call init_sc_outer(pos)

    if(taskid==master) then
        open(10,file="results/init_conf.xyz")
        call writeXyz(D,N,pos,10)
        close(10)
    end if

    !Start force test
    allocate(f(D,N))
    call compute_force_LJ(pos,f,U,P)
    if(taskid==master) then
        open(50,file="results/force_test.dat")
        write(50,"(2(E14.7,2X))")U,P
        do i=1,N
            write(50,"(I3,3(E14.7,2X))")i,f(:,i)
        end do
        close(50)
    end if
    deallocate(f)
    !End force test
    
    !Start g(r) test
    Nshells = 100
    call prepare_shells_and_procs(Nshells,numproc)
    call rad_distr_fun(pos,Nshells)
    if(taskid == master) then
        open(11, file="results/radial_distribution.dat")
        do i=1,Nshells
            write(11,*) (i-1)*grid_shells+grid_shells/2d0,g(i)
        enddo
    endif
    call deallocate_g_variables()
    !End g(r) test

    if (allocated(pos)) deallocate(pos)
    if (allocated(vel)) deallocate(vel)

    call MPI_FINALIZE(ierror)
end program main
