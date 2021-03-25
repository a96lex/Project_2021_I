! AUTHOR: David March

module pbc
  use parameters
  contains  
  subroutine min_img(distv)
  ! finds the closest image of a particle further from L/2 applying PBC:
  ! INPUT: 
  !      distv:    distance vector (x,y,z) between particles
  ! OUTPUT: 
  !      distv:    distance vector (x,y,z) of closest image, if some component is > L/2, else leaves distv unchanged
        implicit none
        real(8), dimension(D), intent(inout) :: distv
        integer i
        do i=1,D
              if(distv(i).gt.(L/2d0)) distv(i) = distv(i) - L
              if(distv(i).lt.(-L/2d0)) distv(i) = distv(i) + L
        enddo
  end subroutine min_img
            
! (currently unused)  		    
  subroutine apply_pbc(pos)
    ! applies pbc to the position matrix once the positions are updated in the integration algorithm
    ! INPUT:
    !      pos:   position matrix
    ! OUTPUT:
    !       pos:   corrected position matrix, if some paticle was outside the box
    implicit none
    real(8), dimension(D,N), intent(inout) :: pos
    integer i, j, factor
    real(8) :: lower_bound, upper_bound
    
    ! Parametres de la caixa de simulacio:
    lower_bound = 0d0
    upper_bound = L
  
    ! Check if the particle has scaped the simulation box, and if so apply PBC:
    do i = 1,N
      do j = 1,D
        if(pos(j,i).gt.upper_bound) then
          ! how much greater? : factor
          factor = int(pos(j,i)/L)
          pos(j,i) = pos(j,i) - factor*L
        endif
        if(pos(j,i).lt.lower_bound) then
          ! how much lower? : factor
          factor = int(dabs(pos(j,i))/L)
          pos(j,i) = pos(j,i) + (factor+1)*L
        endif
      enddo
    enddo
  end subroutine apply_pbc
  
end module
