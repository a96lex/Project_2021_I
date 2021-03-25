    module statvis
        use parameters
        contains

        subroutine writeXyz(dimensions,n_particles,system_positions,file_unit)
        !Author: Alex 
        ! subrutina per a escriure l'estat espaial del sistema 
        ! en format xyz. L'arxiu s'ha d'obrir (amb action="write") 
        ! i tancar al programa des d'on es cridi la subrutina 
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
        include 'mpif.h'
        integer d,i,request,ierror,imin,imax
        double precision vec(d)
        double precision mean,var,tot,totlocal
            ! Passem el vector input a tots els processadors perquè es reparteixin el bucle de la variància
            call MPI_BCAST(vec,d,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
            mean=sum(vec)/dble(d)
            imin = taskid * d / numproc + 1
            imax = (taskid + 1) * d / numproc
            totlocal=0.d0
            tot=0.d0
            do i=imin,imax
                totlocal=totlocal+(vec(i)-mean)**(2.d0)
            enddo
            ! Juntem cada contribució del bucle en una i li passem al master, el qual calcularà la variància total
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(totlocal,tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            if (taskid.eq.master) then
                var=tot/(dble(d)*(dble(d)-1.d0))
            endif
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
        include 'mpif.h'
        integer d,tau,n,request,ierror,tmin,tmax
        character(len=*) :: filename
        double precision vec(d),mean,var,corrlocal(500),corr(500),corsum,time
            ! Passem el vector input a tots els processadors perquè es reparteixin el bucle dels lags
            call MPI_BCAST(vec,d,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
            call estad(size(vec),vec,mean,var)
            call MPI_BCAST(var,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
            tmin = taskid * 500 / numproc + 1
            tmax = (taskid + 1) * 500 / numproc
            corr=0.d0
            corrlocal=0.d0
            ! Bucle dels lags
            do tau=tmin,tmax
                corsum=0.d0
                ! I definició de variància
                do n=1,d-tau
                    corsum=corsum+(vec(n)-mean)*(vec(n+tau)-mean)
                enddo
                corrlocal(tau)=corsum/(dble(d-tau)*var)
            enddo
            ! Juntem cada contribució del bucle en una i li passem al master, el qual escriurà el fitxer
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(corrlocal,corr,500,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            if (taskid.eq.master) then
                open(50,file=filename)
                do tau=1,500
                    write(50,*) tau,corr(tau)
                enddo
                time=1.d0+2.d0*sum(corr)
                write(50,*)
                write(50,*)
                write(50,*) "Integrated Autocorrelation Time =", time
                close(50)
            endif
        return
        end subroutine corrtime
    end module statvis
