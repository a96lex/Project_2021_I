! AUTHOR: David March

module rad_dist
  real(8), parameter :: pi = 4d0*datan(1d0)
  real(8) :: grid_shells
  real(8), dimension(:), allocatable :: g, shells_vect
  
  ! module variables info:
  ! grid_shells will save the width of the spherical shells:
  ! g will save the actual radial distr function, g = N/(dens*V))
  ! shells_vect will save the volume*density of each spherical shell, to avoid computing it at every call of rad_distr_fun
  contains
  
  subroutine prepare_shells(Nshells)
     ! Given the number of shells desired, this subrutine has to be called right before starting the simulation.
     ! Computes/allocates the varialbes declared in the module that will be used in the computation of g(r).
     ! INPUT: --------------------------------------------------------------------------
     !      Nshells:    number of bins where g(r) will be computed
     ! OUTPUT: (module variable) -------------------------------------------------------
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

 subroutine rad_dist_fun(pos,Nshells)
 ! THIS ONE USES i=imin_p,imax_p j=jmin_p(i),jmax_p(i) in the nested loop. Equal number of paris per processor.
 ! computes g(r) in a histogram-like way, saving it in the defined g array
    ! uses the module variable shells_vect
    ! INPUT: ----------------------------------------------------------------------------
    !      pos:    position matrix of the particles
    !      Nshells:    number of bins where to compute g
    ! OUTPUT: (module variable) ---------------------------------------------------------
    !      g:    radial distribution function
    use parameters
    use pbc
    implicit none
    include 'mpif.h'
    integer, intent(in) :: Nshells
    real(8), intent(in) :: pos(D,N)
    ! internal:
    integer i,j,k,part_ini,part_end,ierror,request
    real(8) dist,inner_radius,outer_radius
    real(8), dimension(D) :: distv(D)
    real(8), dimension(N) :: coord
    
    g = 0d0        
    ! Nested loop for i,j pairs only
    do i=imin_p,imax_p
       do j=jmin_p(i),jmax_p(i)
          distv = pos(:,i) - pos(:,j)
          call min_img(distv)
          dist = dsqrt(sum((distv)**2))
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
 end subroutine rad_dist_fun
 

  subroutine deallocate_g_variables()
  ! Deallocate arrays used for computing g. Called before program end.
    deallocate(g)
    deallocate(shells_vect)
  end subroutine

end module rad_dist      
