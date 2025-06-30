module temporal_derivatives

    use time_grid
    use physical_grid
    use make_a_grid

    implicit none
 
contains

!*****************************************************************************************************************
! CONCERN: NOT IMPLEMENTED. DO NOT EVEN LOOK AT THIS. GO TO SPATIAL_DERIVATIVES.F90
!*****************************************************************************************************************


  subroutine temporal_derivative(fun, fd_order, first_derivative, second_derivative)
    
    integer, intent(in) :: fd_order
    double precision, dimension(nx), intent(in) :: fun
    double precision, dimension(nx), intent(inout) :: first_derivative, second_derivative
    integer :: i

  
    ! Compute finite difference derivatives
    select case (fd_order)
    case (2)
      do i = 1+nxghost, nxphys+nxghost
        first_derivative(i) = (fun(i+1) - fun(i-1)) / (2.0 * dt)
        second_derivative(i) = (fun(i-1) - 2.0 * fun(i) + fun(i+1)) / (dt ** 2)
      end do
    case (4)
      do i = 1+nxghost, nxphys+nxghost
        first_derivative(i) = (fun(i-2) - 8.0 * fun(i-1) + 8.0 * fun(i+1) - fun(i+2)) / (12.0 * dt)
        second_derivative(i) = (-fun(i-2) + 16.0 * fun(i-1) - 30.0 * fun(i) &
        + 16.0 * fun(i+1) - fun(i+2)) / (12.0 * (dt ** 2))
      end do
    case (6)
      do i = 1+nxghost, nxphys+nxghost
        first_derivative(i) = (-fun(i-3) + 9.0 * fun(i-2) - 45.0 * fun(i-1) + 45.0 * fun(i+1) &
        - 9.0 * fun(i+2) + fun(i+3)) / (60.0 * dt)
        second_derivative(i) = (2.0 * fun(i-3) - 27.0 * fun(i-2) + 270.0 * fun(i-1) &
        - 490.0 * fun(i) + 270.0 * fun(i+1) - 27.0 * fun(i+2) + 2.0 * fun(i+3)) / (180.0 * (dt ** 2))
      end do
      !       boundary points
      first_derivative(1)=(1./(60.*dt))*(-147.*fun(1)+360.*fun(2)-450.*fun(3)+400.*fun(4)-225.*fun(5)+72.*fun(6)-10.*fun(7))
      first_derivative(2)=(1./(60.*dt))*( -10.*fun(1)- 77.*fun(2)+150.*fun(3)-100.*fun(4)+ 50.*fun(5)-15.*fun(6)+ 2.*fun(7))
      first_derivative(3)=(1./(60.*dt))*(   2.*fun(1)- 24.*fun(2)- 35.*fun(3)+ 80.*fun(4)- 30.*fun(5)+ 8.*fun(6)-    fun(7))
!
!      outer points
!
      first_derivative(nx  )=(1./(60.*dt))*(147.*fun(nx)-360.*fun(nx-1)+450.*fun(nx-2)&
      -400.*fun(nx-3)+225.*fun(nx-4)-72.*fun(nx-5)+10.*fun(nx-6))
      first_derivative(nx-1)=(1./(60.*dt))*( 10.*fun(nx)+ 77.*fun(nx-1)-150.*fun(nx-2)&
      +100.*fun(nx-3)- 50.*fun(nx-4)+15.*fun(nx-5)- 2.*fun(nx-6))
      first_derivative(nx-2)=(1./(60.*dt))*( -2.*fun(nx)+ 24.*fun(nx-1)+ 35.*fun(nx-2)&
      - 80.*fun(nx-3)+ 30.*fun(nx-4)- 8.*fun(nx-5)+    fun(nx-6))



!       do boundary points (also 6th order)
!
      second_derivative(1)=(1./(60.*dt))*(812.*fun(1)-3132.*fun(2)+5265.*fun(3)&
      -5080.*fun(4)+2970.*fun(5)-972.*fun(6)+137.*fun(7))
      second_derivative(2)=(1./(60.*dt))*(137.*fun(1)- 147.*fun(2)- 255.*fun(3)&
      + 470.*fun(4)- 285.*fun(5)+93. *fun(6)- 13.*fun(7))
      second_derivative(3)=(1./(60.*dt))*(-13.*fun(1)+ 228.*fun(2)- 420.*fun(3)&
      + 200.*fun(4)+  15.*fun(5)-12. *fun(6)+  2.*fun(7))
  !
  !  outer points
  !
      second_derivative(nx  )=(1./(60.*dt))*(812.*fun(nx)-3132.*fun(nx-1)+5265.*fun(nx-2)&
      -5080.*fun(nx-3)+2970.*fun(nx-4)-972.*fun(nx-5)+137.*fun(nx-6))
      second_derivative(nx-1)=(1./(60.*dt))*(137.*fun(nx)- 147.*fun(nx-1)- 255.*fun(nx-2)&
      + 470.*fun(nx-3)- 285.*fun(nx-4)+ 93.*fun(nx-5)- 13.*fun(nx-6))
      second_derivative(nx-2)=(1./(60.*dt))*(-13.*fun(nx)+ 228.*fun(nx-1)- 420.*fun(nx-2)&
      + 200.*fun(nx-3)+  15.*fun(nx-4)- 12.*fun(nx-5)+  2.*fun(nx-6))
  ! 
    ! case (8)
    !   do i = 1+nxghost, nxphys+nxghost
    !     first_derivative(i) = (3.0 * fun(i-4) - 32.0 * fun(i-3) + 168.0 * fun(i-2) &
    !     - 672.0 * fun(i-1) + 672.0 * fun(i+1) - 168.0 * fun(i+2) &
    !     + 32.0 * fun(i+3) - 3.0 * fun(i+4)) / (840.0 * dt)
    !     second_derivative(i) = (-9.0 * fun(i-4) + 128.0 * fun(i-3) - 1008.0 * fun(i-2) &
    !     + 8064.0 * fun(i-1) - 14350.0 * fun(i) + 8064.0 * fun(i+1) &
    !     - 1008.0 * fun(i+2) + 128.0 * fun(i+3) - 9.0 * fun(i+4)) / (5040.0 * (dt ** 2))
    !   end do
    ! case (10)
    !   do i = 1+nxghost, nxphys+nxghost
    !     first_derivative(i) = (-2.0 * fun(i-5) + 25.0 * fun(i-4) - 150.0 * fun(i-3) &
    !     + 600.0 * fun(i-2) - 2100.0 * fun(i-1) + 2100.0 * fun(i+1) - 600.0 * fun(i+2) &
    !     + 150.0 * fun(i+3) - 25.0 * fun(i+4) + 2.0 * fun(i+5)) / (2520.0 * dt)
    !     second_derivative(i) = (8.0 * fun(i-5) - 125.0 * fun(i-4) + 1000.0 * fun(i-3) &
    !     - 6000.0 * fun(i-2) + 42000.0 * fun(i-1) - 73766.0 * fun(i) + 42000.0 * fun(i+1) &
    !     - 6000.0 * fun(i+2) + 1000.0 * fun(i+3) - 125.0 * fun(i+4) + 8.0 * fun(i+5)) / (25200.0 * (dt ** 2))
    !   end do !NOTE: The 10th  and 8th order FD scheme is not implemented yet, but can be added after changing the number of ghost zones.


    case default
      print *, 'Invalid order of finite difference'
      return
    end select

  end subroutine temporal_derivative
      
end module temporal_derivatives
