program rober
  use bdf_method_mod
#include "tap.h"
  implicit none
  type(bdf_type) :: BDF
  integer :: i, m
  integer, parameter :: neq=3
  integer, parameter :: nvector=256
  integer, parameter :: rtol_start = 27, rtol_stop=28
  double precision, dimension(nvector,neq) :: y
  double precision :: t_start, t_stop, dt, start, finish, rtol, atol
  real :: eps1 = 1e-4, eps2=1e-8, eps3=1e-4
  character (len=10) :: rtol_char


  call InitBDF(BDF,nvector,neq,rtol,atol)

  ! TAP output
  TAP_PLAN(4*(rtol_stop-rtol_start+1))

  do m=rtol_start,rtol_stop
    BDF%rtol = 10d0**(-2d0-m*0.25d0)      ! relative tolerance
    BDF%atol = 1d-6*BDF%rtol              ! absolute tolerance

    t_start = 0.0d0
    t_stop  = 40d0
    dt = 1d-6

    do i=1,nvector
      y(i,:) = (/ 1d0, 0d0, 0d0 /)
    end do

    call cpu_time(start)

    call SolveODE_BDF(BDF, t_start, dt, t_stop, y)

    call cpu_time(finish)

    print *, BDF%rtol, BDF%atol, finish-start

    ! check the solution
    write(rtol_char,'(F5.2)') log10(BDF%rtol)
    TAP_CHECK_CLOSE(40.0,real(t_start),0.1, "ROBER - time, log(rtol):"//rtol_char)
    TAP_CHECK_CLOSE(0.71582657204273303,real(y(1,1)),eps1,"ROBER - y(1), log(rtol):"//rtol_char)
    TAP_CHECK_CLOSE(9.1855154650076860E-006,real(y(1,2)),eps2,"ROBER - y(2), log(rtol):"//rtol_char)
    TAP_CHECK_CLOSE(0.28416424223264120,real(y(1,3)),eps3,"ROBER - y(3), log(rtol):"//rtol_char)
  end do

  TAP_DONE

  call CloseBDF(BDF)

end program rober
