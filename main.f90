      program main
      use parameters
      use init
      use pbc
      use integraforces
      !AJ: afegir statvis quan sigui un modul

      implicit none
      character(len=50) :: input_name
      real*8, allocatable :: pos(:,:), vel(:,:)

      ! Per executar el programa cal fer >> main.x input_file. Si no, donara error.
      if (command_argument_count() == 0) stop "ERROR: Cridar fent >> ./main.x input_path"
      call get_command_argument(1, input_name)
      
      ! ES: Obrim els arxius necessaris tots a la vegada (?)
      open(unit=10, file=input_name)

      ! Read input parameters
      call get_param(10)
      close(10)

      print*,"------Parameters------"
      print*, "N=",N,"D=",D
      print*,"dt_sim=",dt_sim
      print*,"rho=",rho,"L=",L
      print*,"eps=",epsilon,"sigma=",sigma,"rc=",rc
      print*,"n_meas,n_conf,n_total=",n_meas,n_conf,n_total
      
      ! Allocates
      allocate(pos(D, N))
      allocate(vel(D,N))   

      ! Initialize positions and velocities
      call init_sc(pos)
      call init_vel(vel, 10.d0)  ! Cridem amb temperatura reduida (T'=10) molt alta per fer el melting


      





      ! Deallocates
      deallocate(pos) 
      deallocate(vel)

      end program main
