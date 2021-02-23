      program main
      use parameters
      use init
      use pbc
      use integraforces
      !AJ: afegir statvis quan sigui un modul

      implicit none
      character(len=50) :: input_name
      real*8, allocatable :: pos(:,:) ! ES: Les declarem aqui o a un modul?

      ! Per executar el programa cal fer >> main.x input_file. Si no, es provara de llegir l'arxiu estandar input.txt.
      if (command_argument_count() == 0) then
            input_name = "input.txt" ! Nom de input predeterminat
      else
            call get_command_argument(1, input_name)
      end if

      ! ES: Obrim els arxius necessaris tots a la vegada (?)
      open(unit=10, file=input_name)

      ! Initialize system
      call get_param(10)
      close(10)
      allocate(pos(D, N))  ! ES: Fem el allocate aqui o la poso al inclosa al get_param? En cas del segon, caldria passar tambe
                           ! pos i vel a get_param() o fer un modul on hi siguin.
      call init_sc(pos)
      !call init_vel(vel, 10.)  ! Cridem amb temperatura reduida (T'=10) molt alta per fer el melting

      print*,"------Parameters------"
      print*, "N=",N,"D=",D
      print*,"dt=",dt
      print*,"rho=",rho,"L=",L
      print*,"eps=",epsilon,"sigma=",sigma,"rc=",rc
      print*,"n_meas,n_conf,n_total=",n_meas,n_conf,n_total


      end program main
