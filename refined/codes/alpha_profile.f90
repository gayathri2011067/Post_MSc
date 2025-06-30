module alpha_profile

use parameters
use time_grid
use physical_grid
use make_a_grid

implicit none


!****************************************************************************************************************       
! this is for the kinetic alpha effect only
! once R_alpha=!0, the kinetic alpha effect is turned on
! this module need not be used unless kinetic alpha form needs to be altered
! Currently using the model as in lc+14, i.e., alpha = R_alpha * sin(pi*x)
!****************************************************************************************************************       




    double precision, dimension(nx) :: alpha_k,alpha_cap,alpha, alpha_m

contains

    subroutine construct_alpha_profile

        alpha_cap = sin(pi*x)
        alpha_k = R_alpha*alpha_cap
 

    end subroutine construct_alpha_profile

end module alpha_profile




