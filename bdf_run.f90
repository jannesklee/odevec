program bdf_program
  use bdf_method_mod
  implicit none
  integer :: i, m
  integer, parameter :: neq=3
  integer, parameter :: nvector=256
  double precision, dimension(nvector,neq) :: y
  double precision :: t_start, t_stop, dt, rtol, atol, start, finish

  do m=28,28
    rtol = 10d0**(-2d0-m*0.25d0)   ! relative tolerance
    atol = 1d-6*rtol               ! absolute tolerance

    t_start = 0.0d0
    t_stop = 40.0d0
    dt = 1d-6

    do i=1,nvector
      y(i,:) = (/ 1d0, 0d0, 0d0 /)
    end do

    call cpu_time(start)

    call SolveODE_BDF(nvector, neq, y, t_start, t_stop, dt, rtol, atol)

    call cpu_time(finish)

!    print *, rtol, atol, finish-start

  end do

end program bdf_program
