! AUTHOR: David March

module rad_dist
  real(8), parameter :: pi = 4d0*datan(1d0)
  ! grid_shells will save the width of the spherical shells:
  real(8) :: grid_shells
  ! g will save the actual radial distr function, g = N/(dens*V))
  ! shells_vect will save the volume*density of each spherical shell, to avoid computing it at every call of rad_distr_fun
  real(8), dimension(:), allocatable :: g, shells_vect
  contains
    
    subroutine prepare_shells(Nshells)
      ! Given the number of shells desired, this subrutine has to be called right before starting the simulation.
      ! Computes/allocates the varialbes declared in the module that will be used in the computation of g(r).
      ! INPUT:
      !      Nshells:    number of bins where g(r) will be computed
      ! OUTPUT: (module variable)
      !      shells_vect:    holds the volume*density product of each sperical shell
      use parameters, only : L
      implicit none
      integer, intent(in) :: Nshells
      integer i
      allocate(g(Nshells))
      allocate(shells_vect(Nshells))
      ! Define the grid parameter:
      grid_shells = (L/2d0)/(dble(Nshells))
      ! Compute volume of each shell:
      do i=1,Nshells
        shells_vect(i) = 4d0*pi/3d0 * grid_shells**3 * (i**3 - (i-1)**3)
      enddo
   end subroutine prepare_shells
   
   subroutine rad_distr_fun(pos,Nshells)
    ! computes g(r) in a histogram-like way, saving it in the defined g array
    ! uses the module variable shells_vect
    ! INPUT:
    !      pos:    position matrix of the particles
    !      Nshells:    number of bins where to compute g
    ! OUTPUT: (module variable)
    !      g:    radial distribution function
    use parameters, only : N,D
    use pbc
    implicit none
    integer, intent(in) :: Nshells
    real(8), intent(in) :: pos(D,N)
    ! internal:
    integer i,j,k
    real(8) dist,inner_radius,outer_radius
    real(8), dimension(D) :: distv(D)
    
    g = 0d0
    ! Compute the radial distribution function by averaging the r.d.f. over all particles.
    do i=1,N-1
        do j=i+1,N
            ! Check in wich cell particle j falls repsective from particle i:
            distv = pos(:,i) - pos(:,j)
            call min_img(distv)
            dist = sqrt(sum((distv)**2))
            ! Given dist, add contribution to g(r):
            do k=1,Nshells
                outer_radius = k*grid_shells
                inner_radius = (k-1)*grid_shells
                if(dist.lt.outer_radius.and.dist.gt.inner_radius) then
                    g(k) = g(k) + 1d0/(rho * shells_vect(k))
                endif
            enddo
        enddo
    enddo
    g = 2d0*g/dble(N)
    return
  end subroutine rad_distr_fun
  
  subroutine deallocate_g_variables()
  ! Deallocate arrays used for computing g. Called before program end.
    deallocate(g)
    deallocate(shells_vect)
  end subroutine

end module rad_dist
  
