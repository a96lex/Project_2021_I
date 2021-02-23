module pbc
  use parameters
  contains
 ! ---------------------------------------------------- Minimum image search -----------------------------------------------------------  
  subroutine min_img(x,y,z,a,b,c,dist,L,v)
    ! finds the closest image of a particle further from L/2 applying PBC
    ! x,y,z are the coordinates of the fixed particle
    ! a,b,c are the coordinates of the of the particle whose image to find
    ! dist is the distance between them
    ! L is the side of the simulation box
    ! v(3) returns a vector with the coordinates of the image
    implicit none
    real(8), intent(in) :: x,y,z,a,b,c,L
    real(8), dimension(3), intent(out) :: v
    real(8), intent(inout) :: dist
    real(8) :: dist_prova,signe !AJ, necessari per compilacio
    
    ! Aixo és una funció que es podria passar amb un modul o algo:
    real(8) dist_3d
    v(1) = a
    v(2) = b
    v(3) = c
    ! If the particle lies in the inferior half of the box for the coordinate x (y,z) try negative displacements.
    ! The contrarty if it is on the superior part of the box.
    ! ------------- x ---------------
    if(x.lt.L/2d0) signe=-1d0
    if(x.gt.L/2d0) signe=1d0 ! realment no es donara aquest cas, pero per mantenir la generalitat
    ! Try a displacement on x and accept it if it reduces de distance

    ! dist_prova = dist_3d(x,y,z,a+signe*L,b,c)
    if(dist_prova.lt.dist) then ! accept new x coordinate
      dist=dist_prova
      v(1) = v(1) + signe*L
    endif
    ! ------------- y ---------------
    if(y.lt.L/2d0) signe=-1d0
    if(y.gt.L/2d0) signe=1d0
    ! Try a displacement on y and accept it if it reduces de distance
    ! dist_prova = dist_3d(x,y,z,v(1),b+signe*L,c)

    if(dist_prova.lt.dist) then ! accept new y coordinate
      dist=dist_prova
      v(2) = v(2) + signe*L
    endif
    ! ------------- z ---------------
    if(z.lt.L/2d0) signe=-1d0
    if(z.gt.L/2d0) signe=1d0
    ! Try a displacement on z and accept it if it reduces de distance
    
    ! dist_prova = dist_3d(x,y,z,v(1),v(2),c+signe*L)
    if(dist_prova.lt.dist) then ! accept new z coordinate
      dist=dist_prova
      v(3) = v(3) + signe*L
    endif
    return
  end subroutine min_img
  
  
  subroutine min_img_2(distv)
  ! Subrutina adaptada al calcul de forces (subroutine compute_LJ_force)
        implicit none
        real(8), dimension(D), intent(inout) :: distv
        integer i
        do i=1,D
              if(distv(i).gt.(L/2d0)) distv(i) = distv(i) - L
              if(distv(i).lt.(-L/2d0)) distv(i) = distv(i) + L
        enddo
  end subroutine min_img_2
            
  		    
  
  ! ---------------------------------------------------- PBC -----------------------------------------------------------------------
  subroutine apply_pbc(pos)
    ! applies pbc to the position matrix once the positions are updated in the integration algorithm
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
