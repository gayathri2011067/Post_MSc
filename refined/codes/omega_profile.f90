module omega_profile

use parameters
use time_grid
use physical_grid
use make_a_grid
use eta_profile

implicit none

!*****************************************************************************************************************
!This module is not used at all. Omega is already defined in parameters file, and in this project there is no vertical stratification for omega.
!*****************************************************************************************************************

  double precision, dimension(nx) :: omega_fzr

contains
  subroutine construct_omega_profile

    omega_fzr = 1

  end subroutine construct_omega_profile

end module omega_profile