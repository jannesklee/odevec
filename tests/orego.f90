program orego
  use odevec_main
  use odevec_commons
  implicit none
  integer            :: i, m, k
  integer            :: rtol_start=32, rtol_stop=32
  type(odevec)       :: ode
  double precision   :: t_start, t_stop, start, finish,dt,t_step
  real               :: eps1 = 1e-4, eps2=1e-8, eps3=1e-4
  character (len=10) :: rtol_char


  call InitOdevec(ode)

  do m=rtol_start,rtol_stop
    call cpu_time(start)

    ode%rtol = 10d0**(-2d0-m*0.25d0)          ! relative tolerance
    ode%atol = 1d-6*ode%rtol                  ! absolute tolerance

    t_start = 0.0d0

    do i=1,ode%nvector
      ode%y(i,:) = (/ 1.0, 2.0, 3.0/)
    end do

    t_stop = 1.0
    print *, t_start, ode%y(1,:)
    do while (t_stop  <= 360.0)
      dt = 1d-10
      call SolveODE(ode,t_start, dt, t_stop, ode%y,GetRHS,GetJac,GetLU)
      print *, t_stop, ode%y(1,:)
      t_stop = t_stop + 0.01
    end do

    call cpu_time(finish)

    print *, ""
    print *, "#-------------------------------------------------------#"
    print *, "simulation time:", (finish-start)/1d0
    print *, "used tolerances (rel.,abs.):", ode%rtol, ode%atol
    print *, "#-------------------------------------------------------#"
    print *, ""
  end do

  call CloseOdevec(ode)

end program orego
