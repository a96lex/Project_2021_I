      module integraforces
            use parameters
            contains

            subroutine compute_force_LJ(r,f,U,P)
            ! Computes the force, potential energy and pressure of a system
            ! of N particles. Needs the external parameters M,D,L,rc to work
                  implicit none
                  real*8,intent(in) :: r(D,N)
                  real*8,intent(out) :: f(D,N),U,P
                  real*8 :: distv(D),dist
                  integer :: i,j
                  f = 0.d0
                  U = 0.d0
                  P = 0.d0
                  do i=1,N
                        do j=i+1,N
                              distv = r(:,i)-r(:,j)
                              ! call pbc(distv) !PBC goes here
                              dist = sqrt(sum((distv)**2))
                              if(dist<rc) then
                                    f(:,i) = f(:,j)& 
                                    + (48.d0/dist**14 - 24.d0/dist**8)*distv
                                    f(:,j) = f(:,j)& 
                                    - (48.d0/dist**14 - 24.d0/dist**8)*distv
      
                                    U = U + 4.d0*((1.d0/dist)**12-(1.d0/dist)**6)-&
                                            4.d0*((1.d0/rc)**12-(1.d0/rc)**6)
                                    P = P + sum(distv * f(:,i))
                              endif
                        enddo
                  enddo
                  ! P = P/(dble(N)*(dble(N)-1.d0)/2.d0)

                  P = rho*T_ref + 1.d0/(3.d0*L**3)*P
            end subroutine compute_force_LJ

      end module integraforces