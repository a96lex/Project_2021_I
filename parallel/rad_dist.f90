! AUTHOR: David March

module rad_dist
  real(8), parameter :: pi = 4d0*datan(1d0)
  real(8) :: grid_shells
  real(8), dimension(:), allocatable :: g, shells_vect
  ! this is for testing purposes only:
  real(8), dimension(:), allocatable :: g_p
  
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
      allocate(g_p(Nshells))
      allocate(shells_vect(Nshells))
      ! Define the grid parameter:
      grid_shells = (L/2d0)/(dble(Nshells))
      ! Compute volume of each shell:
      do i=1,Nshells
          shells_vect(i) = 4d0*pi/3d0 * grid_shells**3 * (i**3 - (i-1)**3)
      enddo
   end subroutine prepare_shells
   
 subroutine rad_dist_fun(pos,Nshells)
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
    real(8), dimension(D,N) :: pos_local
    real(8), dimension(Nshells) :: glocal
    
    g = 0d0
    glocal = 0d0
    pos_local = pos
    
    do i=1,3
       coord = pos_local(i,:)
       call MPI_BCAST(coord,N,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
       pos_local(i,:) = coord
    enddo
    
    ! test lines:
    !print*, " about to start g(r) loop, task ", taskid,", ranges: ", imin, imax
    !print*, " is ther something in positions? ", pos_local(1,1), pos_local(2,1), pos_local(3,1)
    
    do i=imin,imax
       do j=1,N
          if(j.ne.i) then
             distv = pos_local(:,i) - pos_local(:,j)
             call min_img(distv)
             dist = dsqrt(sum((distv)**2))
             do k=1,Nshells
                outer_radius = k*grid_shells
                inner_radius = (k-1)*grid_shells
                if(dist.lt.outer_radius.and.dist.gt.inner_radius) then
                   glocal(k) = glocal(k) + 1d0/(rho * shells_vect(k))
                endif
             enddo
          endif
       enddo
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !print*, "Task", taskid, " g(40) ", glocal(40)
    call MPI_REDUCE(glocal,g,Nshells,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
    ! and normalize for N particles
    if(taskid.eq.master) g = g/dble(N)
 end subroutine rad_dist_fun
 
 
 subroutine rad_dist_fun_pairs(pos,Nshells)
 ! THIS ONE USES imin_p,imax_p to use j=i+1,N in the nested loop !
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
    real(8), dimension(D,N) :: pos_local
    real(8), dimension(Nshells) :: glocal
    
    g_p = 0d0
    glocal = 0d0
    pos_local = pos
    
    do i=1,3
       if(taskid.eq.master) coord = pos_local(i,:)
       call MPI_BCAST(coord,N,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
       pos_local(i,:) = coord
    enddo
    

    ! test lines:
    !print*, " about to start g(r) loop (pairs), task ", taskid,", ranges: ", imin_p, imax_p
    !print*, " is ther something in positions? ", pos_local(1,1), pos_local(2,1), pos_local(3,1)
    
    ! Nested loop for i,j pairs only
    do i=imin_p,imax_p
       do j=i+1,N
          distv = pos_local(:,i) - pos_local(:,j)
          call min_img(distv)
          dist = dsqrt(sum((distv)**2))
          do k=1,Nshells
             outer_radius = k*grid_shells
             inner_radius = (k-1)*grid_shells
             if(dist.lt.outer_radius.and.dist.gt.inner_radius) then
                glocal(k) = glocal(k) + 1d0/(rho * shells_vect(k))
             endif
          enddo
       enddo
    enddo
    !print*, "Task", taskid, " g(40) ", glocal(40)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(glocal,g_p,Nshells,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
    ! Add the duplicated contribution left by not counting the j,i pairs in the nested loop, normalize for N particles
    if(taskid.eq.master) g_p = 2d0*g_p/dble(N)
 end subroutine rad_dist_fun_pairs
 

  subroutine deallocate_g_variables()
  ! Deallocate arrays used for computing g. Called before program end.
    deallocate(g)
    deallocate(shells_vect)
  end subroutine

end module rad_dist      
