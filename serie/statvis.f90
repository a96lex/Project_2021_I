    module statvis
        contains
 
        ! Author: Alex
        ! subrutina per a escriure l'estat espaial del sistema 
        ! en format xyz. L'arxiu s'ha d'obrir (amb action="write") 
        ! i tancar al programa des d'on es cridi la subrutina
        subroutine writeXyz(dimensions,n_particles,system_positions,file_unit)
            implicit none
            integer, intent(in)  :: dimensions,n_particles,file_unit
            double precision, intent(in) :: system_positions(dimensions,n_particles)
            integer   :: i
            double precision :: zeros(3-dimensions)
            write(unit=file_unit,fmt=*) n_particles
            write(unit=file_unit,fmt=*)

            if (dimensions < 3) zeros = 0.d0
            do i=1,n_particles
                if (dimensions < 3) then
                    write(unit=file_unit,fmt=*) "He", system_positions(:,i), zeros(:)
                else
                    write(unit=file_unit,fmt=*) "He",system_positions(1,i),system_positions(2,i),system_positions(3,i)
                end if
            enddo    
        end subroutine writeXyz

        
        subroutine estad(d,vec,mean,var)
        !Author: Jaume Ojer
        ! Subrutina que retorna el valor esperat ('mean') i la variància ('var')
        ! d'un vector ('vec') de dimensió 'd'
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

        
        subroutine binning(d,vec,numBINmin,filename)
        !Author: Jaume Ojer
        ! Subrutina que escriu en un fitxer output (de nom 'filename', que és un input)
        ! el valor esperat i la desviació quadràtica a partir de binnejar el vector
        ! input 'vec' de dimensió 'd' (com en la subrutina estad)
        ! Recomano numBINmin=50 (a MoMo em funcionava realment bé)
        implicit none
        integer numBINmin,m,i,j,c,numBIN,test,d
        double precision,allocatable :: binned(:),aux(:)
        character(len=*) :: filename
        double precision total,mean,var,vec(d),var1

            open(10,file=filename)
            ! Definim per primer cop el vector auxiliar
            allocate(aux(size(vec)))
            aux=vec

            ! Comencem pel cas unibinari
            m=1
            call estad(size(vec),vec,mean,var)
            write(10,*) m,mean,dsqrt(var)
            var1=var
            ! Fem per la resta de casos. Iniciem pel cas bibinari
            m=2
            numBIN=size(vec)/m
            ! Els càlculs es duran a terme sempre i quan el número de bins sigui major al número mínim d'aquests
            do while (numBIN > numBINmin)
                ! La dimensió del vector binnejat dependrà de si estem a la meitat del vector total o no
                test=size(aux)-numBIN*2
                if (test.eq.(0)) then
                    allocate(binned(numBIN))
                else
                    allocate(binned(numBIN+1))
                endif

                ! Construïm ja el vector binnejat
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

                ! Escrivim els resultats per cada cas de número de bins
                call estad(size(binned),binned,mean,var)
                write(10,*) m,mean,dsqrt(var)
                ! I redefinim el vector auxiliar com el vector binnejat per tal de poder fer els càlculs per la següent iteració
                deallocate(aux)
                allocate(aux(size(binned)))
                aux=binned
                ! Redefinim també les noves variables per la següent iteració (duplicació dels bins)
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
        
        
	    subroutine corrtime(d,vec,filename)
	    !Author: Jaume Ojer
	    ! Construcció de la funció i temps d'autocorrelació (específic per tau = 500)
        ! El vector 'vec' a analitzar de dimensió 'd' és l'input.
        ! La subrutina escriu en un fitxer de nom 'filename' la funció i temps d'autocorrelació
	    implicit none
	    integer d,tau,n
	    character(len=*) :: filename
	    double precision vec(d),mean,var,corr(500),corsum,time
            open(50,file=filename)
            call estad(size(vec),vec,mean,var)
            ! Fem un bucle per cada un dels lags
            do tau=1,500
                corsum=0.d0
                ! I apliquem la definició de variància amb el lag
                do n=1,d-tau
                    corsum=corsum+(vec(n)-mean)*(vec(n+tau)-mean)
                enddo
              !  if (mod(tau,200).eq.(0)) then
              !    write(*,*) tau
              !  endif
                corr(tau)=corsum/(dble(d-tau)*var)
                write(50,*) tau,corr(tau)
            enddo
            ! Calculem també el temps d'autocorrelació integrat
            time=1.d0+2.d0*sum(corr)
            write(50,*)
            write(50,*)
            write(50,*) "Integrated Autocorrelation Time =", time
            close(50)
	    return
	    end subroutine corrtime
    end module statvis
