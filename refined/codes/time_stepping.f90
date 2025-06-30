module timestepping
    use parameters
    use eta_profile
    use alpha_profile
    use omega_profile
    use initial_field
    use time_grid
    use physical_grid
    use make_a_grid
    use equations

    implicit none
    
!****************************************************************************************************************
!This module implements the RK3 implicit timestepping scheme for the MHD equations.
!It works well, do not change.
!Also, before changing read the paper: Brandenburg, A. (2001). Computational aspects of astrophysical MHD and turbulence. ArXiv. https://doi.org/10.48550/arXiv.astro-ph/0109497
!****************************************************************************************************************


    double precision, dimension(nx) :: k1r,k1phi,k2r,k2phi,k3r,k3phi,k4r,k4phi
    double precision, dimension(nx) :: k1Fr,k2FR,k3Fr,k4Fr
    double precision, dimension(nx) :: k1Fphi,k2Fphi,k3Fphi,k4Fphi
    double precision, dimension(nx) :: k1Er,k2Er,k3Er,k4Er
    double precision, dimension(nx) :: k1Ephi,k2Ephi,k3Ephi,k4Ephi
    double precision, dimension(nx) :: k1alpha,k2alpha,k3alpha,k4alpha
    double precision,dimension(nx) :: g_1, g_2, g_3, z_1, z_2
    contains
        



        subroutine RK3_implicit

            implicit none
            

            double precision, dimension(nx) :: Br_g, Bphi_g, alpha_m_g
            double precision, dimension(nx) :: Fr_g, Fphi_g, Er_g, Ephi_g

            !-----------------------!
            g_1 = 8.0 / 15.0        !
            g_2 = 5.0 / 12.0        !
            g_3 = 3.0 / 4.0         !----> REFER: Brandenburg, A. (2001).ArXiv. 
            z_1 = -17.0 / 60.0      !             Computational aspects of astrophysical MHD and turbulence.  
            z_2 = -5.0 / 12.0       !             https://doi.org/10.48550/arXiv.astro-ph/0109497
            !-----------------------!


            call FOSA(B_r, B_phi, alpha_m)
            k1r = dt * dBrdt
            k1phi = dt * dBphidt
            k1alpha = dt * dalpdt
            B_r = B_r + g_1 * k1r
            B_phi = B_phi + g_1 * k1phi
            alpha_m = alpha_m + g_1 * k1alpha
            Br_g = B_r + z_1 * k1r
            Bphi_g = B_phi + z_1 * k1phi
            alpha_m_g = alpha_m + z_1 * k1alpha
            
            call FOSA(B_r, B_phi, alpha_m)
            k2r = dt * dBrdt
            k2phi = dt * dBphidt
            k2alpha = dt * dalpdt
            B_r = Br_g + g_2 * k2r
            B_phi = Bphi_g + g_2 * k2phi
            alpha_m = alpha_m_g + g_2 * k2alpha
            Br_g = B_r + z_2 * k2r
            Bphi_g = B_phi + z_2 * k2phi
            alpha_m_g = alpha_m + z_2 * k2alpha
            
            call FOSA(B_r, B_phi, alpha_m)
            k3r = dt * dBrdt
            k3phi = dt * dBphidt
            k3alpha = dt * dalpdt
            B_r = Br_g + g_3 * k3r
            B_phi = Bphi_g + g_3 * k3phi
            alpha_m = alpha_m_g + g_3 * k3alpha

            t = t +  dt
            alpha = alpha_m + alpha_k

        end subroutine RK3_implicit










    


end module timestepping