program bdf_program
  use bdf_method_mod
  implicit none
  integer :: i
  integer, parameter :: neq=3
  integer, parameter :: nvector=256
  double precision, dimension(nvector,neq) :: y
  double precision :: t_start, t_stop, dt

  t_start = 0.0d0
  t_stop = 40.0d0
  dt = 1d-6

  do i=1,nvector
    y(i,:) = (/ 1d0, 0d0, 0d0 /)
  end do

  call SolveODE_BDF(nvector, neq, y, t_start, t_stop, dt)

end program bdf_program
