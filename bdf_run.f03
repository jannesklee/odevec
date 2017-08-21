program bdf_program
  use bdf_method_mod
  implicit none
  integer, parameter :: neq=3
  integer, parameter :: nvector=1
  double precision, dimension(neq,neq) :: P
  double precision, dimension(neq) :: r, y
  double precision :: t_start, t_stop, dt

  t_start = 0.0d0
  t_stop = 40.0d0
!  t_stop = 1d4
!  t_stop = 1d11
  dt = 1d-4
  y = (/ double precision :: 1, 0, 0 /)

!  call SolveLinearSystem(nvector,neq,r,y, dt)
  call SolveODE_BDF(neq, nvector, y, t_start, t_stop, dt)

end program bdf_program
