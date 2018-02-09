program primordial
  use bdf_method_mod
  implicit none
  type(odevec) :: ode
  double precision :: t0, time, dt, t_step, t_stop
  integer :: i, j


  call InitBDF(ode)

  ode%rtol = 1e-4
  ode%atol = 1e-20

  ! initialization
  ode%y(:,:) = 1e-20
  do i=1,ode%nvector
    ode%y(i,3) = 0.2d0
    ode%y(i,5) = 0.4d0
!    y(i,:) = (/ 0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, &
!                0d0, 0d0, 0d0, 0d0, 0d0, &
!                0d0, 0d0, 0d0, 0d0, 0d0, 0d0   /)
  end do

  t0 = 0.0d0
  time = t0
  t_step = 1d7
  dt = 1d7

  write (*,*) time, ode%y(:,3), ode%y(:,5)
  do while (time < 1e10)
    t_stop = time + t_step

    call SolveODE_BDF(ode, time, dt, t_stop, ode%y)

    write (*,*) time, ode%y(:,3), ode%y(:,5)
  end do


! implicit none
!  5   integer,parameter :: nsp=krome_nmols !number of species (common)
!  4   real*8            :: Tgas,dt,x(nsp),spy
!  3   real*8 :: t_start, t_stop, time, t_step
!  2 
!  1 
!  1   call krome_init() !init krome (mandatory)
!  2 
!  3   Tgas = 1d3 !gas temperature (K)
!  4 
!  7   ! initialization
!  8   x(:) = 1d-20 !default abundances
!  9   x(krome_idx_H) = 0.2d0
! 10   x(krome_idx_H2) = 0.4d0
! 11   
! 12   time = 0.0d0
! 13   dt = 1d7
! 14 
! 15   write (*,*) time, x(:)
! 16   do while (time < 1e10)
! 17     Tgas = 3d3
! 18     
! 19     call krome(x(:), Tgas, dt) !call KROME
! 20 
! 21     time = time + dt
! 22     
! 23     write (*,*) time, x(:)
! 24 
! 25   end do



!  call CloseBDF()

end program primordial
