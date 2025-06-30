
!******************************************************************************************************************************************************
!This module constructs the one dimensional grid for simulation.
!DO NOT EDIT THE PHYSICAL GRID PARAMETERS UNLESS YOU KNOW WHAT YOU ARE DOING.
!Time grid can be editted as required for different time resolutions.
!******************************************************************************************************************************************************

module time_grid

use parameters
implicit none

  integer, parameter :: Nt= 50000.                                          !points per diffusion time
  integer, parameter :: n1= 5000.                                           !Number of snapshots
  double precision, parameter :: total_t= 30.                               !unit diffusion time
  integer :: growth_index_guess = n1 / (2. * total_t)                       !idx for calculating the slope at 0.5 unit time
  integer :: n2 = nint(total_t*Nt/n1)                                       !Number of timesteps between snapshots
  double precision, parameter :: dt= 1./Nt                                  !time step
  double precision :: t=0.                                                  !time variable initialization
  double precision :: first=0.                                              !for Runge-Kutta routine

end module time_grid



!**************************************************************************************************************************************************

module physical_grid

use parameters
use time_grid
implicit none

  integer, parameter :: nxphys= 101                                        !Resolution in z
  integer, parameter :: nxghost= 3                                         !Number of ghost cells at each end in z
  integer, parameter :: mid_x= (nxphys+1)/2                                !midpoint in z
  integer, parameter :: nx= nxphys +2*nxghost                              !Resolution in z
  double precision, dimension(nx) :: x
  double precision :: dx, alp, beta

end module physical_grid

!******************************************************************************************************************************************************


module make_a_grid

implicit none
integer :: i

contains

  subroutine construct_grid

    use parameters
    use time_grid
    use physical_grid


    double precision, parameter :: len= 2.*h
    double precision, dimension(nx) :: space

    dx=len/(nxphys-1)                                                       !x corresponds to z
    alp=dt/dx**2                                                            !alpha for finite difference scheme
    beta= dt/dx                                                             !beta for fin diff 

    do i=1,nx
      x(i)= -(h+nxghost*dx) +(i-1)*dx                                       !+0.01 !NOTE:+0.01 to avoid 0
      x(i)= x(i)                                                            !dimensionless now
    enddo

  endsubroutine construct_grid

end module 

!******************************************************************************************************************************************************