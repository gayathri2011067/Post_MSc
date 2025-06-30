module initial_field
use parameters
use time_grid
use physical_grid
use make_a_grid
use eta_profile
use alpha_profile
use velocity_profile
use seed
use spatial_derivatives
implicit none

!****************************************************************************************************************
! This module initializes the magnetic field and other related variables like turbulent velocity profile and eta.
! Note that the case of velocity stratification is controlled by "A" parameter in parameters.f90.
!****************************************************************************************************************


double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq
double precision, dimension(nx) :: Fr, Fphi, Er, Ephi, Tee, Phii, dphii, d2phii,small_u!,small_u_0
double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br, B_0,eta
double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi
double precision, dimension(nx) :: term1, term2, term3,analytic_B
double precision, dimension(n1) :: log_B_mid, d_log_B_mid_dt, d2_log_B_mid_dt2
double precision :: Bseed
double precision, dimension(nx) :: l
double precision, dimension(nx) :: rho,rho_0


contains

    subroutine field_initialization
        call construct_alpha_profile
        call construct_velocity_profile
        call init_random_seed
        call spatial_derivative(phii, 6, dphii, d2phii)

            Bseed     = 100.  !CONCERN: Bseed is 100 because we are getting 10e-3 data for this number

            !constant rho  !CONCERN: I am not at all sure about this part of the code. I am just setting it to 1 in dimensionless units. 
            !in dimension units rho is not defined at all. !CONCERN: I SHOULD COME BACK TO THIS.
            rho_0     = 1.
            rho       = 1.

            small_u   = small_u_0*exp(A*x**2/2.)
            l = tau*small_u
            eta       = (1.*tau/(3.))*small_u**2
            
            B_eq      = 4.*pi*rho*small_u**2
            B_0       = 4.*pi*rho_0*small_u_0**2
            
            !*****************INITIAL FIELDS********************
            B_r       = xseed*exp(-x**2)*(1.-x**2)*Bseed
            B_phi     = xseed*exp(-x**2)*(1.-x**2)*Bseed
            !***************************************************


    end subroutine field_initialization
    
end module initial_field