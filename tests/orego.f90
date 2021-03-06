!##############################################################################
!#                                                                            #
!#    Copyright (C) 2019                                                      #
!#    Jannes Klee <jklee@astrophysik.uni-kiel.de>                             #
!#                                                                            #
!#    This file is part of OdeVec.                                            #
!#                                                                            #
!#    OdeVec is free software: you can redistribute it and/or modify          #
!#    it under the terms of the GNU General Public License as published by    #
!#    the Free Software Foundation, either version 3 of the License, or       #
!#    (at your option) any later version.                                     #
!#                                                                            #
!#    OdeVec is distributed in the hope that it will be useful,               #
!#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
!#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
!#    GNU General Public License for more details.                            #
!#                                                                            #
!#    You should have received a copy of the GNU General Public License       #
!#    along with OdeVec.  If not, see <https://www.gnu.org/licenses/>.        #
!#                                                                            #
!##############################################################################
!> \brief OREGO test - a repeating network
!!
!! \author Jannes Klee
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

    ode%RelativeTolerance = 10d0**(-2d0-m*0.25d0)          ! relative tolerance
    ode%AbsoluteTolerance = 1d-6*ode%RelativeTolerance                  ! absolute tolerance

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
    print *, "used tolerances (rel.,abs.):", ode%RelativeTolerance, ode%AbsoluteTolerance
    print *, "#-------------------------------------------------------#"
    print *, ""
  end do

  call CloseOdevec(ode)

end program orego
