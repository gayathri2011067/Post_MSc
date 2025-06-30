module seed
  
!****************************************************************************************************************
!This module initializes the random seed for the simulation.
!It is used to generate random numbers for various purposes in the simulation.
!****************************************************************************************************************
use parameters
use time_grid
use physical_grid
use make_a_grid
implicit none

    double precision, dimension(nx) :: xseed
    contains
      subroutine init_random_seed
        call random_number(xseed)
      end subroutine init_random_seed

end module seed