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
            call MPI_REDUCE(totlocal,tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            if (taskid.eq.master) then
                var=tot/(dble(d)*(dble(d)-1.d0))
            endif
        return
        end subroutine estad
        
        
        subroutine estadcoeff(d,vec,mean,var)
        !Author: Jaume Ojer
        ! Subrutina específica pels coeficients de difusió que retorna 
        ! el valor esperat ('mean') i la variància ('var') d'un vector ('vec') de dimensió 'd'
        implicit none
        include 'mpif.h'
        integer d,i,request,ierror,imin,imax
        double precision vec(d)
        double precision mean,var,tot,totlocal
            mean=sum(vec)/dble(d)
            imin = taskid * d / numproc + 1
            imax = (taskid + 1) * d / numproc
            totlocal=0.d0
            tot=0.d0
            do i=imin,imax
                totlocal=totlocal+(vec(i)-mean)**(2.d0)
            enddo
            ! Juntem cada contribució del bucle en una i li passem al master, el qual calcularà la variància total
            call MPI_REDUCE(totlocal,tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            if (taskid.eq.master) then
                var=tot/(dble(d)*(dble(d)-1.d0))
            endif
        return
        end subroutine estadcoeff
        
        
        subroutine binning(d,vec,numBINmin,file_unit)
        !Author: Jaume Ojer
        ! Subrutina que escriu en un fitxer output (de nom 'filename', que és un input)
        ! el valor esperat i la desviació quadràtica a partir de binnejar el vector
        ! input 'vec' de dimensió 'd' (com en la subrutina estad)
        ! Recomano numBINmin=50 (a MoMo em funcionava realment bé)
        implicit none
        include 'mpif.h'
        integer numBINmin,file_unit,m,i,j,c,numBIN,test,d,dimbin,request,ierror
        double precision,allocatable :: binned(:),aux(:)
        double precision total,mean,var,vec(d),var1
        
            ! Definim per primer cop el vector auxiliar
            allocate(aux(d))
            aux=vec

            ! Comencem pel cas unibinari
            m=1
            call estad(d,vec,mean,var)
            if (taskid.eq.master) then
                write(file_unit,*) m,mean,dsqrt(var)
                var1=var
            endif
            ! Fem per la resta de casos. Iniciem pel cas bibinari
            m=2
            numBIN=d/m
            ! Els càlculs es duran a terme sempre i quan el número de bins sigui major al número mínim d'aquests
            do while (numBIN > numBINmin)
                ! La dimensió del vector binnejat dependrà de si estem a la meitat del vector total o no
                if (taskid.eq.master) then
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
                    dimbin = size(binned)
                endif

                ! Escrivim els resultats per cada cas de número de bins
                call MPI_BCAST(dimbin,1,MPI_INTEGER,master,MPI_COMM_WORLD,request,ierror)
                if (taskid.ne.master) then
                    allocate(binned(dimbin))
                endif
                call MPI_BCAST(binned,dimbin,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
                call estad(size(binned),binned,mean,var)
                
                if (taskid.eq.master) then
                    write(file_unit,*) m,mean,dsqrt(var)
                endif
                ! I redefinim el vector auxiliar com el vector binnejat per tal de poder fer els càlculs per la següent iteració
                deallocate(aux)
                allocate(aux(size(binned)))
                aux=binned
                ! Redefinim també les noves variables per la següent iteració (duplicació dels bins)
                deallocate(binned)
                m=m*2
                numBIN=size(aux)/2
            enddo
            
            if (taskid.eq.master) then
                write(file_unit,*)
                write(file_unit,*)
                write(file_unit,*) "Autocorrelation time =", var/var1
            endif
        return
        end subroutine binning
        
        
        subroutine corrtime(d,vec,file_unit)
        !Author: Jaume Ojer
        ! Construcció de la funció i temps d'autocorrelació (tau = d/10)
        ! El vector 'vec' a analitzar de dimensió 'd' és l'input.
        ! La subrutina escriu en un fitxer de nom 'filename' la funció i temps d'autocorrelació
        implicit none
        include 'mpif.h'
        integer d,file_unit,tau,n,lag,request,ierror,tmin,tmax
        double precision vec(d),mean,var,corsum,time
        double precision,allocatable :: corr(:),corrlocal(:)
            ! Passem el vector input a tots els processadors perquè es reparteixin el bucle dels lags
            call MPI_BCAST(vec,d,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
            call estad(size(vec),vec,mean,var)
            call MPI_BCAST(mean,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,request,ierror)
            lag = d/10
            allocate(corrlocal(lag))
            allocate(corr(lag))
            tmin = taskid * lag / numproc + 1
            tmax = (taskid + 1) * lag / numproc
            corr=0.d0
            corrlocal=0.d0
            ! Bucle dels lags
            do tau=tmin,tmax
                corsum=0.d0
                ! I definició de variància
                do n=1,d-tau
                    corsum=corsum+(vec(n)-mean)*(vec(n+tau)-mean)
                enddo
                corrlocal(tau)=corsum/dble(d-tau)
            enddo
            ! Juntem cada contribució del bucle en una i li passem al master, el qual escriurà el fitxer
            call MPI_REDUCE(corrlocal,corr,lag,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierror)
            if (taskid.eq.master) then
                write(file_unit,*) 0, 1.d0
                do tau=1,lag
                    write(file_unit,*) tau, (corr(tau))/(var*dble(d))
                enddo
                time=1.d0+2.d0*(sum(corr))/(var*dble(d))
                write(file_unit,*)
                write(file_unit,*)
                write(file_unit,*) "Integrated Autocorrelation Time =", time
            endif
            deallocate(corrlocal)
            deallocate(corr)
        return
        end subroutine corrtime
    end module statvis
