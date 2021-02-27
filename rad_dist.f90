module radial_distribution
  real(8), parameter :: pi = 4d0*datan(1d0)
  ! grid_shells will save the width of the spherical shells:
  real(8) :: grid_shells
  ! rad_distr will save the number of particles on each shell
  ! shells_vect will save the volume of each spherical shell
  real(8), dimension(:), allocatable :: rad_distr, shells_vect
  contains
    
    subroutine prepare_shells(Nshells)
      ! Given the number of shells desired, this subrutine has to be called right before starting the simulation.
      ! Computes/allocates the varialbes declared in the module that will be used in the computation of g(r)
      use parameters, only : L, rho
      implicit none
      integer, intent(in) :: Nshells
      integer i
      allocate(rad_distr(Nshells))
      allocate(shells_vect(Nshells))
      ! Define the grid parameter:
      grid_shells = L/(2*dble(Nshells))
      ! Compute volume of each shell:
      do i=1,Nshells
        shells_vect(i) = 4d0*pi/3d0 * grid_shells**3 * rho * (i**3 - (i-1)**3)
      enddo
   end subroutine prepare_shells
   
   
   subroutine rad_distr_fun(pos)
    ! computes g(r) in a histogram-like way, saving it in the defined rad_distr array
    ! pass the shell volumes * density in a vector to avoid calculating them every time! (denominator, g = N/(dens*V))
    use parameters, only : N,L,D
    implicit none
    real(8), intent(in) :: pos(D,N)
    ! internal:
    integer i,j,k
    real(8) dist,distv,inner_radius,outer_radius
    
    rad_distr = 0d0
    ! Compute the radial distribution function by averaging the r.d.f. over all particles.
    do i=1,N
      do j=1,N
        if(i.ne.j) then
          ! Check in wich cell particle j falls repsective from particle i:
          distv = pos(:,i) - pos(:,j)
          call min_img_2(distv)
          dist = sqrt(sum((distv)**2))
          do k=1,Nshells
            outer_radius = k*grid_space
            inner_radius = (k-1)*grid_space
            if(dist.lt.outer_radius.and.dist.gt.inner_radius) then
              g(k) = g(k) + (1d0/shell_vols(k))/dble(N) ! normalize to the number of particles, divide by N
            endif
          enddo
        endif
      enddo
    enddo
    
    return
  end subroutine rad_distr_fun
end module radial_distribution
  
