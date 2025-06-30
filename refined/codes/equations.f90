module equations


      use parameters
      use eta_profile           !not used in simulation
      use velocity_profile      !not used since advection is set to zero
      use alpha_profile          
      use omega_profile         !not used since omega profile is assumed to be constant vertically
      use initial_field
      use time_grid
      use physical_grid
      use make_a_grid
      use spatial_derivatives
      use temporal_derivatives



!*****************************************************************************************************************
! This module contains the dynamical quenching model equations for dynamo, FOSA approximation used.
! The induction equations are as per lc+14, and helicity evolution equation is lc+14 modified by including the new helicity fluxes
!*****************************************************************************************************************

      implicit none
      double precision, dimension(nx) :: dBrdt, dBphidt, dalpdt, dFrdt, dFphidt, dErdt, dEphidt 

      contains



      subroutine FOSA(B_r_dummy, B_phi_dummy, alpha_m_dummy)

            double precision, intent(inout), dimension(nx) :: B_r_dummy, B_phi_dummy, alpha_m_dummy
            double precision,  dimension(nx) :: d_alpha_Bphi, d_alpha_Br
            double precision,  dimension(nx) :: d2Br , dBr, dBphi, d2Bphi
            double precision,  dimension(nx) :: d_Uz_Bphi, d_Uz_Br
            double precision,  dimension(nx) :: Uz_Bphi, Uz_Br
            double precision,  dimension(nx) :: alpha_Bphi, alpha_Br
            double precision, dimension(nx) :: d_alpha_m, d2_alpha_m
            double precision, dimension(nx) :: alpha_total,vishniac_general,new_diff_general
            double precision, dimension(nx) :: vishniac_term,new_diffusive_flux, random_advective_flux
            double precision, dimension(nx) :: coeff1_fd, coeff2_fd,alpm_1_fd, alpm_2_fd
            double precision, dimension(nx) :: fra_1, fra_2, fra_3, hb
            double precision, dimension(nx) :: der_u, d_sq_u
            double precision, dimension(nx) :: d2_alpha_Br, d2_alpha_Bphi
            double precision, dimension(nx) :: d2_Uz_Br, d2_Uz_Bphi
            double precision :: R_alpha_inv = 0.

            character(len=30) :: ghost_zone_type = 'anti-symmetric'
            character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
            character(len=30) :: ghost_zone_type3 = 'symmetric'
            character(len=30) :: ghost_zone_type4 = 'makezero'

            call impose_boundary_conditions(B_r_dummy, ghost_zone_type)
            call impose_boundary_conditions(B_phi_dummy, ghost_zone_type)
            call impose_boundary_conditions(alpha_m_dummy, ghost_zone_type)

            alpha_total = alpha_m_dummy + R_alpha*alpha_cap
            alpha_Bphi = alpha_total * B_phi_dummy
            alpha_Br = alpha_total * B_r_dummy

            Uz_Br = B_r_dummy *U_z
            Uz_Bphi = B_phi_dummy *U_z

            call spatial_derivative(B_r_dummy, 6, dBr, d2Br)
            call spatial_derivative(B_phi_dummy, 6, dBphi, d2Bphi)
            call spatial_derivative(alpha_Br, 6, d_alpha_Br, d2_alpha_Br)
            call spatial_derivative(alpha_Bphi, 6, d_alpha_Bphi, d2_alpha_Bphi)
            call spatial_derivative(Uz_Br, 6, d_Uz_Br, d2_Uz_Br)
            call spatial_derivative(Uz_Bphi, 6, d_Uz_Bphi, d2_Uz_Bphi)
            call spatial_derivative(alpha_m_dummy, 6, d_alpha_m, d2_alpha_m)
            call spatial_derivative(U_z, 6, der_u, d_sq_u)





                  !-----------------------------VISHNIAC TERM-----------------------------------------------------------------------------!
                  !                                                                                                                       !
                  ! -----> This is the vishniac flux term, derived from gs+23.                                                            ! 
                  ! -----> automatic switching between cases                                                                              !
                  !                                                                                                                       !
                  if (A == 1) then                                                                                                        !
                  vishniac_term = (0.079 - (7.0/45.0)*xi_0*exp(-2.0*x**2)) * (small_u_0**4 * tau**2 * 2.0*x*xi_0) * (R_omega / 9.0)       !
                  else if (A == -1) then                                                                                                  !
                  vishniac_term = (-0.156*xi_0 + 0.154) * 2.0*x*(eta**2)*xi_0*R_omega*exp(-2.0*x**2)         
                  else if (A == -3) then
                  vishniac_term = -(0.311*xi_0*exp(2.0*x**2) - 0.693) * x*(eta**2)*xi_0*R_omega*exp(2.0*x**2)         
                  else                                                                                                                    !
                  vishniac_term = 0.0                                                                                                     !
                  end if                                                                                                                  !
                  !-----------------------------------------------------------------------------------------------------------------------!




                  !-----------------------------NEW DIFFUSION TERM------------------------------------------------------------------------!
                  !                                                                                                                       !
                  ! -----> This is the vishniac flux term, derived from gs+23.                                                            ! 
                  ! -----> automatic switching between cases                                                          !  
                  !                                                                                                                       !
                  if (A == 1) then  
                  ! xi_diff = xi_0                                                                                                          !                 
                  coeff1_fd = ((exp(x**2)-exp(-x**2)*xi_diff)*(2.*x)*small_u_0**2*7.*tau)/(27.)                                           !
                  coeff2_fd = ((7.*small_u_0**2*tau)/(27.))*(exp(x**2)+exp(-x**2)*xi_diff)                                                ! 
                  alpm_1_fd = 4.*pi*3.*rho*tau*small_u**2*(d_alpha_m + 2.*x*alpha_m_dummy)/(36.*pi*eta*rho)                               !
                  alpm_2_fd = 4.*pi*3.*rho*tau*small_u**2*(d2_alpha_m + 4.*x*d_alpha_m + 2.*alpha_m_dummy +&                              !
                  4.*x**2*alpha_m_dummy)/(36.*pi*eta*rho)                                                                                 !
                  new_diffusive_flux = (coeff1_fd*alpm_1_fd + coeff2_fd*alpm_2_fd)                                                        !
                  else if (A == -1) then 
                  ! xi_diff = xi_0                                                                                                     !
                  coeff1_fd = -(exp(-x**2)*(1+xi_diff)*(2.*x)*small_u_0**2*7.*tau)/(27.)                                                  !
                  coeff2_fd = ((7.*small_u_0**2*tau)/(27.))*(exp(-x**2))*(1+xi_diff)                                                      !  
                  alpm_1_fd = 4.*pi*3.*rho*tau*small_u**2*(d_alpha_m - 2.*x*alpha_m_dummy)/(36.*pi*eta*rho)                               !
                  alpm_2_fd = 4.*pi*3.*rho*tau*small_u**2*(d2_alpha_m - 4.*x*d_alpha_m - 2.*alpha_m_dummy +&                              !
                  4.*x**2*alpha_m_dummy)/(36.*pi*eta*rho)                                                                                !
                  new_diffusive_flux = (coeff1_fd*alpm_1_fd + coeff2_fd*alpm_2_fd)                                                                                                                         !
                  else if (A == -3) then  
                  ! xi_diff = xi_0
                  coeff1_fd = -14.*x*eta*(3. + xi_diff*exp(2.*x**2))  /(9.)                                                                                                   !
                  coeff2_fd = 7.*eta*(1. + xi_diff*exp(2.*x**2))/(9.)  
                  alpm_1_fd = -6.*x*alpha_m_dummy + d_alpha_m
                  alpm_2_fd = -12.*x*d_alpha_m - 6.*alpha_m_dummy + 36.*x**2*alpha_m_dummy + d2_alpha_m       
                  new_diffusive_flux = (coeff1_fd*alpm_1_fd + coeff2_fd*alpm_2_fd)                                                                                                                         !
                  else                                                                                                                    !
                  new_diffusive_flux = 0.0                                                                                                !
                  end if                                                                                                                  !
                  !                                                                                                                       !
                  ! ----------------------------------------------------------------------------------------------------------------------!    

vishniac_general = 2.*R_omega*eta**2*x*xi_0*exp(-1.*(1+A)*x**2)*(A*C2-C1*xi_0*exp(-1.*(1+A)*x**2)-C3-C4) ! This is the vishniac term in general, not used in the code but can be used for future reference
new_diff_general = (7.*eta/9.)*(1+xi_0*exp(-1.*(1+A)*x**2))*(d2_alpha_m - 4.*x**2*A**2*alpha_m_dummy + 4.*A*x*d_alpha_m)&
+(7.*eta/9.)*(A-xi_0*exp(-1.*(1+A)*x**2))*(2.*x)*(d_alpha_m + 2.*A*x*alpha_m_dummy)! This is the new diffusive flux in general, not used in the code but can be used for future reference



                  !CONCERN: Random advective flux is not verified (ie, the derivation is not finalized yet), so we always set it to zero
                  ! !-----------------------------RANDOM ADVECTION---------------------------------
                  !       hb = 3.*rho*tau*small_u**2*alpha_m_dummy*4.*pi/(36.*pi*eta*rho)
                  !       fra_1 = (tau*small_u**2/18.)*(2.+ 4.*x**2) + (7.*tau/27.)*exp(-x**2)*xi_diff*small_u_0**2*(4.*x**2-2.)
                  !       fra_2 = (2*x*small_u_0**2)*((tau/18.)*exp(x**2)-(7.*tau/27.)*exp(-x**2)*xi_diff)
                  !       fra_3 = (tau/18.)*(41./5.)*R_omega*eta
                  !       random_advective_flux = (-fra_1*hb - fra_2*alpm_2_fd + fra_3*(d_alpha_m*hb + alpm_2_fd*alpha_m_dummy ))
                  ! !----------------------------------------------------------------------------------------------------    





                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<   THE DYNAMICAL QUENCHING MODEL   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
                  !                                                                                                                         !
                  dBrdt = -d_alpha_Bphi +eta*d2Br - d_Uz_Br                                                                                 !
                  !                                                                                                                         !
                  if (alp_square) then                                                                                                      !
                  dBphidt = R_omega*B_r_dummy + eta*d2Bphi - d_Uz_Bphi + d_alpha_Br                                                         !
                  else                                                                                                                      !
                  dBphidt = R_omega*B_r_dummy + eta*d2Bphi - d_Uz_Bphi                                                                      !
                  end if                                                                                                                    !
                  !                                                                                                                         !
                  !   
                  if (use_general) then 
                  dalpdt =    (-2./(3.*tau))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)/B_eq**2&                                          !
                  - eta*(B_phi_dummy*dbr-B_r_dummy*dBphi)/B_eq**2+ R_m_inv*alpha_m_dummy) &                                                       !
                  -U_z*d_alpha_m- alpha_m_dummy*der_u &                                                                                     !
                  + vishniac_general*vishniac_term_control&                                                                                    !
                  + new_diff_general*new_diffusive_flux_control &                                                                         !
                  +R_k*d2_alpha_m*old_diffusion_control 
                  else                                                                                                                            !
                  dalpdt =    (-2./(3.*tau))*((alpha_total)*(B_r_dummy**2+B_phi_dummy**2)/B_eq**2&                                          !
                  - eta*(B_phi_dummy*dbr-B_r_dummy*dBphi)/B_eq**2+ R_m_inv*alpha_m_dummy) &                                                       !
                  -U_z*d_alpha_m- alpha_m_dummy*der_u &                                                                                     !
                  + vishniac_term*vishniac_term_control&                                                                                    !
                  + new_diffusive_flux*new_diffusive_flux_control &                                                                         !
                  +R_k*d2_alpha_m*old_diffusion_control  
                  end if                                                                                   !
                  ! ------------------------------------------------------------------------------------------------------------------------!    



      end subroutine FOSA

end module equations

