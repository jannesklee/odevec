PRogram primordial
  use odevec_main
  use odevec_commons
  implicit none
  type(odevec) :: ode
  double precision :: t0, time, dt, t_step, t_stop
  double precision :: start,finish
  integer :: i, j

  call cpu_time(start)

  call InitOdevec(ode)

  ode%rtol = 1d-4
  ode%atol = 1d-20
  ode%dt_min = 1d-20

  ! initialization
  ode%y(:,:) = 1d-20
  do i=1,ode%nvector
    ode%y(i,3) = 0.2d0
    ode%y(i,5) = 0.8d0
  end do

  t0 = 0.0d0
  time = t0
  t_step = 1d8
  dt = 1d8*1e-6
  t_stop = time + t_step

  write (*,*) time, ode%y(1,:)
  do while (time < 1e10)

    call SolveODE(ode,time,dt,t_stop,ode%y,GetRHS,GetJac,GetLU)

    t_stop = t_stop + t_step

    write (*,*) time, ode%y(1,:)
  end do
  call cpu_time(finish)

  print *,
  print *, "#-------------------------------------------------------#"
  print *, "simulation time:", (finish-start)/1d0
  print *, "used tolerances (rel.,abs.):", ode%rtol, ode%atol
  print *, "#-------------------------------------------------------#"
  print *,


  call CloseOdevec(ode)

end program primordial
