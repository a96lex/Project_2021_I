program main
    implicit none
    call plotEnergy(10,10,"dades/test.dat")
end program main


! subrutina per a escriure l'estat espaial del sistema 
! en format xyz. L'arxiu s'ha d'obrir (amb action="write") 
! i tancar al programa des d'on es cridi la subrutina 
subroutine writeXyz(D,N,r,unit)
    implicit none
    integer, intent(in)  :: D,N,unit
    double precision, intent(in) :: r(D,N)
    integer   :: i
    write(unit=unit,fmt=*) N
    write(unit=unit,fmt=*)
    do i=1,N
        write(unit=unit,fmt=*) "He",r(i,1),r(i,2),r(i,3)
    enddo    
end subroutine writeXyz

! encara no funciona. No aconsegueixo fer plots amb més d'una
! linea temporal a fortran. Potser es millor fer gnuplot normal
subroutine plotEnergy(timesteps,unit,filename)
    implicit none
    integer, intent(in) :: timesteps,unit
    character(len=100), intent(in) :: filename
    real*8 :: Ecin, Epot, Etot
    integer :: i

    open(unit=unit, file=trim(filename), action="read")
   
    print '(a)', 'unset key'             ! Disable legend.
    print '(a)', 'set title "Plot test"' ! Set title.
    print '(a)', 'plot "-" using 0:3 with lines, "" using 0:2 with lines'
    
    do i = 1, timesteps
        read(unit,*) Ecin, Epot, Etot
        write(*, *) Ecin, Epot, Etot
    end do

    close(unit=unit)

end subroutine plotEnergy