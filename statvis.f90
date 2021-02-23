    module statvis
        contains
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

	    ! Subrutina que retorna el valor esperat ('mean') i la variància ('var')
            ! d'un vector ('vec') de dimensió 'd'
	    subroutine estad(d,vec,mean,var)
	    implicit none
	    integer d,i
	    double precision vec(d)
	    double precision mean,var,tot
  	      mean=sum(vec)/dble(d)
  	      tot=0.d0
  	      do i=1,d
    	        tot=tot+(vec(i)-mean)**(2.d0)
  	      enddo
  	      var=tot/(dble(d)*(dble(d)-1.d0))
	    return
	    end subroutine estad

	    ! Subrutina que escriu en un fitxer output (de nom 'filename', que és un input)
            ! el valor esperat i la desviació quadràtica a partir de binnejar el vector
            ! input 'vec' de dimensió 'd' (com en la subrutina estad)
            ! Recomano numBINmin=50 (a MoMo em funcionava realment bé)
	    subroutine binning(d,vec,numBINmin,filename)
	    implicit none
	    integer numBINmin,m,i,j,c,numBIN,test,d
	    double precision,allocatable :: binned(:),aux(:)
	    character(len=*) :: filename
	    double precision total,mean,var,vec(d),var1

	    ! Obrim el fitxer en el que escriurem els resultats
	      open(10,file=filename)
	    ! Definim per primer cop el vector auxiliar
	      allocate(aux(size(vec)))
	      aux=vec

	    ! Comencem pel cas unibinari
	      m=1
	      call estad(size(vec),vec,mean,var)
	    ! I escribim els resultats
	      write(10,*) m,mean,dsqrt(var)
	      var1=var
	    ! Fem per la resta de casos. Iniciem pel cas bibinari
	      m=2
	      numBIN=size(vec)/m
	    ! Els càlculs es duran a terme sempre i quan el número de bins sigui major al número mínim d'aquests
	      do while (numBIN > numBINmin)
	    !   La dimensió del vector binnejat dependrà de si estem a la meitat del vector total o no
	        test=size(aux)-numBIN*2
	        if (test.eq.(0)) then
	          allocate(binned(numBIN))
	        else
	          allocate(binned(numBIN+1))
	        endif

	    !   Construïm ja el vector binnejat
	        c=0
	        binned=0.d0
	        do i=1,numBIN
	          total=0.d0
	          do j=1,2
	            c=c+1
	            total=total+aux(c)
	          enddo
	          binned(i)=total/2.d0
	        enddo
	        do i=1,test
	          c=c+1
	          binned(numBIN+1)=binned(numBIN+1)+aux(c)/dble(test)
	        enddo

	    !   Escrivim els resultats per cada cas de número de bins
	        call estad(size(binned),binned,mean,var)
	        write(10,*) m,mean,dsqrt(var)
	    !   I redefinim el vector auxiliar com el vector binnejat per tal de poder fer els càlculs per la següent iteració
	        deallocate(aux)
	        allocate(aux(size(binned)))
	        aux=binned
	    !   Redefinim també les noves variables per la següent iteració (duplicació dels bins)
	        deallocate(binned)
	        m=m*2
	        numBIN=size(aux)/2
	      enddo
	      write(10,*)
	      write(10,*)
	      write(10,*) "Autocorrelation time =", var/var1
	      close(10)
	    return
	    end subroutine binning

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
    end module statvis
