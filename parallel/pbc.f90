! AUTHOR: David March

module pbc
  use parameters
  contains  
  subroutine min_img_2(distv)
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
  end subroutine min_img_2  
end module
