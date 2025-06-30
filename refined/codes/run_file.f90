program run_all

!************************************************************************
! This is where everything runs. The labour force(?).
! Add all print statements and everything here. 
!************************************************************************

use alpha_profile
use velocity_profile
use initial_field
use seed
use parameters
use time_grid
use physical_grid
use make_a_grid
use eta_profile
use omega_profile
use equations
use spatial_derivatives
use temporal_derivatives
use timestepping

implicit none



  character(len=32) :: arg1
  integer :: kk, j
  character(len=30) :: data_path, filename, xfile, omegafile, alphafile, Br_ini_file,analytic,&
  B_r_final_file, B_phi_ini_file, B_phi_final_file, time_file, alpham_final_file,turb_vel_file,&
  vishniac_term_file,alp_tot_file, alp_k_file, log_B_mid_file, d_log_B_mid_dt_file, d2_log_B_mid_dt2_file


    call construct_grid
    call construct_eta_profile
    call construct_omega_profile
    call construct_alpha_profile
    call field_initialization


    ! Define the output file names
    data_path = '../run_files/'
    filename =  'eta_fz_values.txt'
    xfile=  'z_values.txt'
    omegafile=  'omega_values.txt'
    alphafile=  'alpha_values.txt'
    Br_ini_file=  'Br_ini.txt'
    B_phi_ini_file=  'B_phi_ini.txt'
    B_r_final_file=  'Br_final.txt'
    B_phi_final_file=  'B_phi_final.txt'
    time_file=  'time.txt'
    alpham_final_file=  'alpham_final.txt'
    turb_vel_file=  'turb_vel.txt'
    vishniac_term_file=  'vishniac_term.txt'
    alp_tot_file = 'alp_tot.txt'
    alp_k_file = 'alp_k.txt'
    log_B_mid_file = 'log_B_mid.txt'
    d_log_B_mid_dt_file = 'd_log_B_mid_dt.txt'
    d2_log_B_mid_dt2_file = 'd2_log_B_mid_dt2.txt'
    analytic=  'analytic.txt'

    ! Open the file for writing
    open(unit=10, file=trim(data_path) // filename)
    open(unit=17, file=trim(data_path) // xfile)
    open(unit=19, file=trim(data_path) // alphafile)
    open(unit=20, file=trim(data_path) // Br_ini_file)
    open(unit=21, file=trim(data_path) // B_phi_ini_file)
    open(unit=22, file=trim(data_path) // B_r_final_file)
    open(unit=23, file=trim(data_path) // B_phi_final_file)
    open(unit=24, file=trim(data_path) // time_file)
    open(unit=25, file=trim(data_path) // alpham_final_file)
    open(unit=26, file=trim(data_path) // turb_vel_file)
    open(unit=27, file=trim(data_path) // vishniac_term_file)
    open(unit=28, file=trim(data_path) // alp_tot_file)
    open(unit=29, file=trim(data_path) // alp_k_file)
    open(unit=30, file=trim(data_path) // log_B_mid_file)
    open(unit=31, file=trim(data_path) // d_log_B_mid_dt_file)
    open(unit=32, file=trim(data_path) // d2_log_B_mid_dt2_file)
    open(unit=33, file=trim(data_path) // analytic)


    do i = 1, nx
        write(10, '(F12.8)') eta_fz(i)
        write(17, '(F12.8)') x(i) 
        write(19, '(F12.8)') alpha_k(i)
        write(20, '(F12.8)') B_r(i)
        write(21, '(F12.8)') B_phi(i)
        write(26, '(F12.8)') small_u(i)
        write(33, '(F12.8)') analytic_B(i)


    end do

    ! Close these files
    close(10)
    close(17)
    close(19)
    close(20)
    close(21)
    close(26)
    close(33)




print *, 't_d ', t_d_dim
print*, 'const', R_k* pi**7  /(R_omega**2* 32.*sqrt(2.)*(eta_0_dim)**2)
print*, 'alp', (R_omega* sqrt(2.)*xi_0*eta_0_dim*(-0.155*xi_0 +0.079))/(R_k* pi**2* h_dim)


!**************************************************************************************!
!                           𝕺𝖓𝖊 𝖑𝖔𝖔𝖕 𝖙𝖔 𝖗𝖚𝖑𝖊 𝖙𝖍𝖊𝖓 𝖆𝖑𝖑                                        !
!**************************************************************************************!
do kk = 1, n1 ! for n1 iterations                                                      !
    ! print*, 'iteration ', kk -1, 'completed'                                         !
    do j = 1, n2 ! for n2 time steps                                                   !
      call RK3_implicit                                                                !
    end do                                                                             !
    call temporal_derivative(log_B_mid, 6, d_log_B_mid_dt, d2_log_B_mid_dt2)           !
      write (22, *) B_r/B_0                                                            !
      write (23, *) B_phi/B_0                                                          !
      write (24, *) t*t_d_dim                                                          !
      write (25, *) alpha_m                                                            !
      write (28, *) alpha                                                              !
      write (29, *) alpha_k                                                            !
      write (30, *) log_B_mid                                                          !
      write (31, *) d_log_B_mid_dt                                                     !
      write (32, *) d2_log_B_mid_dt2                                                   !
end do                                                                                 !
!**************************************************************************************!


close(22)
close(23)
close(24)
close(25)
close(28)
close(29)
close(30)
close(31)
close(32)

print*, 't=', t
print *, 'File sucessfully run'








end program run_all
