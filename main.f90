      program main
      use parameters
      use init
      use pbc
      use integraforces
      !AJ: afegir statvis quan sigui un modul

      call get_param("input_template.txt") !AJ: canviar a unit

      print*,"------Parameters------"
      print*, "N=",N,"D=",D
      print*,"dt=",dt
      print*,"rho=",rho,"L=",L
      print*,"n_meas,n_conf,n_total=",n_meas,n_conf,n_total


      end program main