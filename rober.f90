program rober
  use odevec_main
  use odevec_commons
  implicit none
  integer :: i, m, k
  type(odevec) :: ode
  integer :: rtol_start=0, rtol_stop=0
  double precision :: t_start, t_stop, start, finish,dt,t_step
  real :: eps1 = 1e-4, eps2=1e-8, eps3=1e-4
  character (len=10) :: rtol_char


  call InitOdevec(ode)

  do m=rtol_start,rtol_stop
    ode%rtol = 10d0**(-2d0-m*0.25d0)          ! relative tolerance
    ode%atol = 1d-6*ode%rtol                  ! absolute tolerance

    call cpu_time(start)
    do k=1,1
      t_start = 0.0d0
      t_stop  = 1d10
      dt = 1d-5

      do i=1,ode%nvector
        ode%y(i,:) = (/ 1d0, 0d0, 0d0 /)
      end do

      t_step = 1d-4
      print *, t_start, ode%y(1,:)
      do while (t_stop  < 1e12)
        t_stop = t_start + t_step
        dt = 1d-6
        call SolveODE(ode,t_start, dt, t_stop, ode%y,GetRHS,GetJac,GetLU)
        t_start = t_start
        t_step = t_step*1.3
        print *, t_stop, ode%y(1,:)
      end do

    end do
    call cpu_time(finish)

!    print *, ode%rtol, ode%atol, (finish-start)/1d0
  end do

  call CloseOdevec(ode)

end program rober
