program primordial
  use odevec_main
  use odevec_commons
  implicit none
  type(odevec) :: ode
  double precision :: t0, time, dt, t_step, t_stop
  integer :: i, j


  call InitOdevec(ode)

  ode%rtol = 1e-4
  ode%atol = 1e-20

  ! initialization
  ode%y(:,:) = 1e-20
  do i=1,ode%nvector
    ode%y(i,3) = 0.2d0
    ode%y(i,5) = 0.4d0
  end do

  t0 = 0.0d0
  time = t0
  t_step = 1d7
  dt = 1d7

  write (*,*) time, ode%y(1,:)
  do while (time < 1e10)
    t_stop = time + t_step

    call SolveODE(ode,time,dt,t_stop,ode%y,GetRHS,GetJac,GetLU)

    write (*,*) time, ode%y(1,:)
  end do


!  call CloseBDF()

end program primordial
