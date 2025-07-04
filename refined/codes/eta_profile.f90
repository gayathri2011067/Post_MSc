module eta_profile
  use parameters
  use time_grid
  use physical_grid
  use make_a_grid

  implicit none

  double precision, dimension(nx) :: eta_fz 
  double precision, parameter :: eta__0 = 1 ! dimensionless
  double precision, parameter :: eta_1 = 0.95

!****************************************************************************************************************  
!This module is not used at all. Eta is already defined in field initial.
!****************************************************************************************************************

contains

    subroutine construct_eta_profile

        eta_fz=eta__0*(1-eta_1*exp(-(x)**2))

    end subroutine construct_eta_profile

end module eta_profile