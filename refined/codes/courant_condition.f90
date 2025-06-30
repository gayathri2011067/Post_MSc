module cfl_checker
  use physical_grid
  use time_grid
  use parameters
  use initial_field
  implicit none


  ! double precision,parameter :: vmax=1.
  ! double precision, dimension(nx) :: v
  ! double precision :: cfl_max
  ! double precision :: v0

  ! !====== CONFIGURE THIS SHIT ======!
  ! ! v0 = 1.0D0  ! <- adjust this value if needed

  ! !====== BUILD VELOCITY PROFILE ======!




  ! !====== CFL CALCULATION ======!
  ! cfl_max = vmax * dt / dx
  ! print *, "===================="
  ! print *, "  CFL CHECKER"
  ! print *, "===================="
  ! print *, "dt      = ", dt
  ! print *, "dx      = ", dx
  ! print *, "v_max   = ", vmax
  ! print *, "CFL_max = ", cfl_max
  ! print *, "===================="

  ! if (cfl_max > 1.0D0) then
  !    print *, "WARNING: CFL > 1.0 — you mad bastard."
  !    print *, "  Your scheme might survive, but accuracy will cry."
  ! else if (cfl_max > 0.5D0) then
  !    print *, "Moderate CFL — you're in the spicy-but-okay zone."
  ! else
  !    print *, "CFL is safe. You're a boring, responsible adult."
  ! endif

end module cfl_checker
