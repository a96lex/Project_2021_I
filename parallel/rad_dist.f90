! AUTHOR: David March

module rad_dist
  real(8), parameter :: pi = 4d0*datan(1d0)
  ! grid_shells will save the width of the spherical shells:
  real(8) :: grid_shells
  ! g will save the actual radial distr function, g = N/(dens*V))
  ! shells_vect will save the volume*density of each spherical shell, to avoid computing it at every call of rad_distr_fun
  real(8), dimension(:), allocatable :: g, shells_vect
  integer, dimension(:,:), allocatable :: ranges_proc
  contains
    
    subroutine prepare_shells_and_procs(Nshells,avail_proc)
      ! Given the number of shells desired, this subrutine has to be called right before starting the simulation.
      ! Computes/allocates the varialbes declared in the module that will be used in the computation of g(r).
      ! Sets the particle range that each processor has to compute
      ! INPUT:
      !      Nshells:    number of bins where g(r) will be computed
      !      avail_proc:  number of available processors
      ! OUTPUT: (module variable)
      !      shells_vect:    holds the volume*density product of each sperical shell
      use parameters, only : N, L, numproc
      implicit none
      integer, intent(in) :: Nshells, avail_proc
      integer i,j
      real(8) total_pairs, pairs_per_proc, sum_pairs
      real(8), dimension(N) :: num_pairs
      allocate(g(Nshells))
      allocate(shells_vect(Nshells))
      ! Define the grid parameter:
      grid_shells = (L/2d0)/(dble(Nshells))
      ! Compute volume of each shell:
      do i=1,Nshells
          shells_vect(i) = 4d0*pi/3d0 * grid_shells**3 * (i**3 - (i-1)**3)
      enddo
      ! Distribute particles bewteen processors (aprox equal number per processor):
      do i=1,N
          num_pairs(i) = N-i
      enddo
      allocate(ranges_proc(avail_proc,2))
      total_pairs = sum(dble(num_pairs))
      pairs_per_proc = total_pairs/dble(numproc)
      ranges_proc(1,1) = 1
      ranges_proc(numproc,2) = N-1
      do i=1,numproc-1
          sum_pairs = 0d0
          limits: do j=ranges_proc(i,1),N
              sum_pairs = sum_pairs + dble(num_pairs(j))
              if(sum_pairs.gt.pairs_per_proc) then
                  ranges_proc(i,2) = j
                  ranges_proc(i+1,1) = j+1
                  exit limits
              endif
          enddo limits
      enddo
    ! Check lines below: could be erased  
    ! do i=1,numproc
    !     print*, "Limits for proc ",i,":", ranges_proc(i,:), "Num paris: ",sum(num_pairs(ranges_proc(i,1):ranges_proc(i,2)))
    ! enddo
   end subroutine prepare_shells_and_procs
   
   subroutine rad_distr_fun(pos,Nshells)
    ! computes g(r) in a histogram-like way, saving it in the defined g array
    ! uses the module variable shells_vect
    ! INPUT:
    !      pos:    position matrix of the particles
    !      Nshells:    number of bins where to compute g
    ! OUTPUT: (module variable)
    !      g:    radial distribution function
    use parameters, only : N,L,D,rho,numproc,taskid,master
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
    part_ini = ranges_proc(taskid+1,1)
    part_end = ranges_proc(taskid+1,2)
    pos_local = pos
    ! Cast the position vector to all processors:
    do i=1,3
        coord = pos_local(i,:)
        call MPI_BCAST(coord,N,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
        pos_local(i,:) = coord
    enddo
    ! Compute the radial distribution function by averaging the r.d.f. over all particles.
    do i=part_ini,part_end
        do j=i+1,N
            if(i.ne.j) then
                ! Check in wich cell particle j falls repsective from particle i:
                distv = pos(:,i) - pos(:,j)
                call min_img_2(distv)
                dist = sqrt(sum((distv)**2))
                ! Given dist, add contribution to g(r):
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
    g = 2d0*g/dble(N)    
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    ! Check lines below:
    !print*, "taskid ",taskid," g local 20 30 40", glocal(20), glocal(30), glocal(40)
    !print*, "taskid",taskid,"g",glocal(:)
    
    call MPI_REDUCE(glocal,g,Nshells,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
  end subroutine rad_distr_fun
  
  subroutine deallocate_g_variables()
  ! Deallocate arrays used for computing g. Called before program end.
    deallocate(g)
    deallocate(shells_vect)
  end subroutine

end module rad_dist
  
