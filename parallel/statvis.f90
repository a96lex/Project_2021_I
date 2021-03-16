    module statvis
        contains
            ! subrutina per a escriure l'estat espaial del sistema 
            ! en format xyz. L'arxiu s'ha d'obrir (amb action="write") 
            ! i tancar al programa des d'on es cridi la subrutina 
            ! Author: Alex 
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
    end module statvis
