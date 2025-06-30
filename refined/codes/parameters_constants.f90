module parameters
implicit none

!***************************************************************************************************************

!---------> this module contains the all the parameter values used for simulation.
!---------> it also contains the conversion factors for the dimensional values to dimensionless values
!---------> NOTE: Only change the parameter values given in "CONTROL PARAMETERS" section, rest will always be kept the same

!***************************************************************************************************************


! CONSTANT VALUES ♥ for conversions, all the dimensional values are in kpc and gyr terms

      double precision, parameter :: pi= 3.14156295358
      double precision, parameter :: g_mp=   1.67262158d-24
      double precision, parameter :: cm_kpc= 3.08567758d21
      double precision, parameter :: km_kpc= 3.08567758d16
      double precision, parameter :: cm_km=  1.d5
      double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600
      double precision, parameter :: G_muG= 1.d6

!***************************************************************************************************************

! WITH DIMENSION ♥ 

      double precision, parameter :: small_u_0_dim = 20.*km_kpc/s_Gyr ! km/s.kpc --> 1/Gyr
      double precision, parameter :: tau_dim = 10.*0.001 !Gyr ! TODO_LATER: Find correct value
      double precision, parameter :: eta_0_dim = (tau_dim*small_u_0_dim**2)/3. ! km/s.kpc --> 1/Gyr
      ! double precision, parameter :: small_u_0_dim = 45.*km_kpc/s_Gyr ! km/s.kpc --> 1/Gyr
      !above expression for new u from ss21
      double precision, parameter :: radius_dim = 4.d0 ! kpc
      double precision, parameter :: r_d_dim = 10.d0 ! kpc
      double precision, parameter :: h_d_dim = 0.5d0 ! kpc
      double precision, parameter :: h_dim = 0.5!h_d_dim*(sqrt(1+(radius_dim/r_d_dim)**2)) !kpc !disc flaring
      double precision, parameter :: omega_0_dim = 127.*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
      double precision, parameter :: r__omega_dim = 2. ! kpc
      double precision, parameter :: l_dim = 0.1 ! kpc
      double precision, parameter :: omega_dim = omega_0_dim*(1.+(radius_dim/r__omega_dim)**2.)**(-0.5) ! 1/Gyr
      double precision, parameter :: G_dim = -45.6*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
      double precision, parameter :: alpha_0_dim = 1.50 ! kpc/Gyr
      double precision, parameter :: k_dim = 0.3*s_Gyr/km_kpc !km.kpc/s --> kpc**2/Gyr
      double precision, parameter :: R_dim = 20.!kpc
      double precision, parameter :: z_i_dim = -h_dim!kpc
      double precision, parameter :: z_f_dim = +h_dim !kpc
      double precision, parameter :: B_0_dim = 0.82/G_muG !muG to G 
      double precision, parameter :: t_d_dim = 0.78445507463348940 !h_dim**2/eta_0_dim! Gyr

!***************************************************************************************************************

! DIMENSIONLESS ♥

      double precision, parameter :: radius = radius_dim/h_dim
      double precision, parameter :: r_d = r_d_dim/h_dim
      double precision, parameter :: h_d = h_d_dim/h_dim
      double precision, parameter :: h = h_dim/h_dim
      double precision, parameter :: eta_0 = eta_0_dim/eta_0_dim
      double precision, parameter :: t_d = t_d_dim/(h_dim**2./eta_0_dim)
      double precision, parameter :: omega_0 = omega_0_dim*h_dim/(eta_0_dim) 
      double precision, parameter :: r__omega = r__omega_dim/h_dim
      ! double precision, parameter :: l = l_dim/h_dim
      double precision, parameter :: omega = omega_dim*h_dim/(eta_0_dim)
      double precision, parameter :: G = G_dim*(h_dim**2/(eta_0_dim))
      double precision, parameter :: alpha_0 = alpha_0_dim*h_dim/eta_0_dim
      double precision, parameter :: k = k_dim/(eta_0_dim) 
      double precision, parameter :: R = R_dim/h_dim
      double precision, parameter :: z_i = z_i_dim/h_dim
      double precision, parameter :: z_f = z_f_dim/h_dim      
      double precision, parameter :: R_m_inv = 0. 
      double precision, parameter :: tau = tau_dim/(h_dim**2./eta_0_dim)
      double precision, parameter :: small_u_0 = small_u_0_dim*h_dim/eta_0_dim

!****************************************************************************************************************


! CONTROL PARAMETERS ♥ 
!---------> these are the parameters that can be changed to control the simulation

      double precision, parameter :: R_alpha = 1.                      !to set the kinetic alpha effect, set as 0 to turn off the effect
      double precision, parameter :: R_omega = -20.                    !to set the differential rotation, fiducial is set to be -20, refer to lc+14(toolbox)
      double precision, parameter :: Dynamo_number = R_alpha * R_omega !dynamo number
      
      double precision, parameter :: C1 = 7./45.  
      double precision, parameter :: C2 = -203./5400.    
      double precision, parameter :: C3 = 403./8100.  
      double precision, parameter :: C4 = -1./6.    

      double precision, parameter :: R_k = 0.3                         !fickian diffusion as givden in lc+14, this matters only if old_diffusion_control is set to 1
      double precision, parameter :: R_U = 0.                          !to set the advection term, set as 0 to turn off the effect

      logical :: use_general = .true.
      double precision, parameter :: A = -3.                            !To set the cases, 1 for case 1, -1 for case 2 and so on

      double precision, parameter :: new_diffusive_flux_control = 1.   !for new diffusion form as given in gs+23
      double precision, parameter ::old_diffusion_control =0.          !for fickian diffusion
      double precision, parameter :: random_advective_flux_control= 0. !for random advection term as given in gs+23. NOTE: this is not verified in the current version of the code, but can be used in future
      double precision, parameter :: vishniac_term_control = 1.        !for vishniac term as given in gs+23. Set to zero for conventional dynamo action

      double precision :: xi_0 = 6.0 !control of the small scale magnetic field strength ON VISHNIAC FLUX EXPRESSION
      double precision :: xi_diff    =0.0        !control of the small scale magnetic field strength ON NEW DIFFUSION EXPRESSION

      logical :: alp_square = .false.


!****************************************************************************************************************


end module parameters


