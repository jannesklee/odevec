program bdf_program
  use bdf_method_mod
  implicit none
  type(bdf_type) :: BDF
  integer :: i, m
  integer, parameter :: neq=3
  integer, parameter :: nvector=1
  double precision, dimension(nvector,neq) :: y
  double precision :: t_start, t_stop, dt, rtol, atol, start, finish


  call InitBDF(BDF,nvector,neq)


  do m=28,28
    rtol = 10d0**(-2d0-m*0.25d0)   ! relative tolerance
    atol = 1d-6*rtol               ! absolute tolerance

    t_start = 0.0d0
    t_stop  = 40d0
    dt = 1d-6

    do i=1,nvector
      y(i,:) = (/ 1d0, 0d0, 0d0 /)
    end do

    call cpu_time(start)

    call SolveODE_BDF(BDF, t_start, t_stop, dt, rtol, atol, y)

    call cpu_time(finish)

!    print *, rtol, atol, finish-start

  end do

  call CloseBDF(BDF)

end program bdf_program
