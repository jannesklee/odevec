Program primordial
  use odevec_main
  use odevec_commons
  implicit none
  type(odevec) :: ode
  integer, parameter :: nmols=12,nspec=16
  integer, parameter :: nvector=8192
  double precision :: t0, time, dt, t_step, t_stop
  double precision :: start,finish
  double precision, pointer :: rho_c(:),x(:,:)
  integer :: i, j

  call cpu_time(start)

  call InitOdevec(ode)

  allocate(rho_c(ode%nvector), x(ode%nvector,nmols))

  ode%rtol = 1d-4
  ode%atol = 1d-20
  ode%dt_min = 1d-20

  x(:,:) = 1e-20

  do i=1,ode%nvector
    x(i,3)  = 0.75615    !H
    x(i,1)  = 4.46406e-8 !E
    x(i,8)  = 8.1967e-5  !H+
    x(i,6)  = 1.0e-40    !D
    x(i,11) = 1.0e-40    !D+
    x(i,4)  = 0.24375    !He
    x(i,9)  = 1.0e-40    !He+
    x(i,10) = 1.0e-40    !H2+
    x(i,5)  = 1.5123e-6  !H2
    x(i,7)  = 1.0e-40    !HD
    x(i,2)  = 1.0e-40    !H-
    x(i,12) = 1.0e-40    !He++
  end do

  do i=1,ode%nvector
    x(i,:) = x(i,:) / sum(x(i,:))
  end do


  t0 = 0.0d0
  time = t0
  t_step = 1d8
!  dt = t_step*1e-6
  dt = 2.2954626201190183E-026
  t_stop = time + t_step

  do while (time < 1e10)
    rho_c(:) = 1e-10

    do i=1,ode%nvector
      ode%y(i,:) = krome_x2n(x(i,:),rho_c(i))
    end do
    write (*,*) time, ode%y(1,:)

    call SolveODE(ode,time,dt,t_stop,ode%y,GetRHS,GetJac,GetLU)

    do i=1,ode%nvector
      x(i,:) = krome_n2x(ode%y(i,:),rho_c(i))
    end do


    t_stop = t_stop + t_step

  end do
  write (*,*) time, ode%y(1,:)
  call cpu_time(finish)

  print *, ""
  print *, "#-------------------------------------------------------#"
  print *, "simulation time:", (finish-start)/1d0
  print *, "used tolerances (rel.,abs.):", ode%rtol, ode%atol
  print *, "#-------------------------------------------------------#"
  print *, ""


  deallocate(rho_c, x)


  call CloseOdevec(ode)

  contains

  !return an array sized krome_nmols containing
  ! the mass fractions (#), computed from the number
  ! densities (1/cm3) and the total density in g/cm3
  function krome_n2x(n,rhogas)
    implicit none
    real*8 :: n(nmols),krome_n2x(nmols)
    real*8,value :: rhogas

    krome_n2x(:) = n(:) * krome_get_mass() / rhogas

  end function krome_n2x

  !********************
  !return an array sized krome_nmols containing
  ! the number densities (1/cm3), computed from the mass
  ! fractions and the total density in g/cm3
  function krome_x2n(x,rhogas)
    implicit none
    real*8 :: x(nmols)
    real*8 :: krome_x2n(nmols)
    real*8,value :: rhogas

    !compute densities from fractions
    krome_x2n(:) = rhogas * x(:) * krome_get_imass()

  end function krome_x2n

  !get an array of double containing the inverse
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !get inverse of the species masses (1/g)
  function get_imass()
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 1.09776932527d+27  !E
    get_imass(2) = 5.9721335838d+23 !H-
    get_imass(3) = 5.97538433901d+23  !H
    get_imass(4) = 1.49430705554d+23  !HE
    get_imass(5) = 2.9876921695d+23 !H2
    get_imass(6) = 2.98861411108d+23  !D
    get_imass(7) = 1.99220448934d+23  !HD
    get_imass(8) = 5.97863863505d+23  !H+
    get_imass(9) = 1.4945104915d+23 !HE+
    get_imass(10) = 2.98850552203d+23 !H2+
    get_imass(11) = 2.98942796572d+23 !D+
    get_imass(12) = 1.49471398286d+23 !HE++
    get_imass(13) = 0.d0  !CR
    get_imass(14) = 0.d0  !g
    get_imass(15) = 0.d0  !Tgas
    get_imass(16) = 0.d0  !dummy

  end function get_imass

  function krome_get_mass()
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !get species masses (g)
  function get_mass()
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 9.10938188d-28  !E
    get_mass(2) = 1.67444345638d-24 !H-
    get_mass(3) = 1.67353251819d-24 !H
    get_mass(4) = 6.69206503638d-24 !HE
    get_mass(5) = 3.34706503638d-24 !H2
    get_mass(6) = 3.34603251819d-24 !D
    get_mass(7) = 5.01956503638d-24 !HD
    get_mass(8) = 1.67262158d-24  !H+
    get_mass(9) = 6.69115409819d-24 !HE+
    get_mass(10) = 3.34615409819d-24  !H2+
    get_mass(11) = 3.34512158d-24 !D+
    get_mass(12) = 6.69024316d-24 !HE++
    get_mass(13) = 0.d0 !CR
    get_mass(14) = 0.d0 !g
    get_mass(15) = 0.d0 !Tgas
    get_mass(16) = 0.d0 !dummy

  end function get_mass


end program primordial
