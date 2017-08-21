module bdf_method_mod
  implicit none
contains

  subroutine SolveODE_BDF(neq, nvector, y, t, t_stop, dt)
    implicit none
    integer :: neq,nvector, order, iterator
    double precision, dimension(neq) :: y, rhs, el
    double precision, dimension(neq,6) :: y_old
    double precision :: t, t_old, t_stop, dt
    double precision, dimension(neq) :: dy, res
    intent(inout) :: y

    ! load initial condition into Nordsieck array
    order = 1
    y_old(:,1) = y

    iterator = 0
    do while (t <= t_stop)
      print *, t, y


      call SolveLinearSystem(nvector,neq,order,y,y_old,dt)
      t = t + dt

!      call CalcOrder(order,iterator)

      iterator = iterator + 1
    end do
  end subroutine


  subroutine SolveLinearSystem(nvector,neq,order,y,y_old,h)
    implicit none
    integer :: neq, nvector, order
    double precision, dimension(neq) :: dy, res, y, rhs, dy_old, y_null, el, rhs_init, del, y_safe
    double precision, dimension(neq,6) :: y_old
    double precision, dimension(neq,neq) :: L, U
    double precision, dimension(7,6) :: coeff
    double precision :: theta, sigma, tol, error, absdy, h, t, beta
    integer :: i, j, iterator
    intent(inout) :: y_old, y
    intent(in) :: nvector, neq, h

! NORDSIECK VERSION
    coeff = reshape((/ double precision :: &
          1.       , 1. , 0.        , 0.       , 0.        , 0.      , 0.       , &
          2./3.    , 1. , 1./3.     , 0.       , 0.        , 0.      , 0.       , &
          6./11.   , 1. , 6./11.    , 1./11.   , 0.        , 0.      , 0.       , &
          12./25.  , 1. , 7./10.    , 1./5.    , 1./50.    , 0.      , 0.       , &
          60./137. , 1. , 225./274. , 85./274. , -15./274. , 1./274. , 0.       , &
          20./49.  , 1. , 58./63.   , 5./12.   , -25./252. , 1./84.  , 1./1764. /), &
          shape(coeff))

    ! general variables
    sigma = 0.01d0
    beta = 1.0d0
    tol = 1d-12
    error = HUGE(tol)

    ! initialize solution arrays
    y_null = y

    ! use the LU matrices from the predictor value
    call GetLU(neq,y,t,beta,h,L,U)

    ! Calculate initial right hand side
    call CalcRHS(neq,y_old(:,1),rhs_init)
    y_old(:,2) = h*rhs_init

    ! some initializations
    iterator = 0
    dy_old = 1.0
    del = 0.0


    ! predictor step
    call PredictSolution(neq,order,y_old,el)

!    do while (error > sigma*tol)
    do while (iterator < 10)

      ! calculates residuum
      call CalcResiduum(nvector,neq,y,y_null,y_old,el,rhs_init,order,h,beta,res)

      ! calculates the solution for dy with given residuum
      call SolveLU(neq,L,U,res,del)

      ! add to error
      el = el + del

      dy = coeff(1,order)*el

      y = y_null + dy

      ! additional stuff
      theta = norm2(dy)/norm2(dy_old)
      error = theta/(1.-theta)*norm2(dy)
      dy_old = dy
      iterator = iterator + 1
    end do

    ! rewrite history array
    call UpdateNordsieck(neq, order, coeff, el, y_old)
  end subroutine

  subroutine UpdateNordsieck(neq, order, coeff, el, y_old)
    integer :: neq,i,j, order
    double precision, dimension(neq) :: el
    double precision, dimension(neq,6) :: y_old
    double precision, dimension(7,6) :: coeff
    double precision :: h, one, et

    do j = 1,order+1
      do i = 1,neq
        y_old(i,j) = y_old(i,j) + el(i)*coeff(order,j)
      end do
    end do
  end subroutine

  subroutine PredictSolution(neq,order,y_old,el)
    implicit none
    integer :: i,j,k,order,neq
    double precision, dimension(neq,6) :: y_old
    double precision, dimension(neq) :: el
    intent (in) :: neq,order
    intent (inout) :: y_old
    intent (out) :: el

    ! this loop effectively solves the Pascal Triangle without multiplications
    do k = 0,order-1
      do j = order,k-1
        y_old(:,j-1) = y_old(:,j) + y_old(:,j-1)
      end do
    end do
    el = 0.0
  end subroutine


  !> Solve the System \f$ LUx=r \f$
  subroutine SolveLU(neq,L,U,res,dy)
    implicit none
    integer :: neq
    double precision, dimension(neq,neq) :: L, U
    double precision, dimension(neq) :: dy, res
    double precision, dimension(neq) :: y
    intent(in) :: res, L, U, neq
    intent(out) :: dy
    integer :: j,k

    y = 0.0
    ! lower system
    do k = 1,neq,1
      do j = 1,k-1
        y(k) = y(k) + L(k,j)*y(j)
      end do
      y(k) = res(k) - y(k)
      y(k) = y(k)/L(k,k)
    end do

    ! upper system
    dy = 0.0
    do k = neq,1,-1
      do j = k+1,neq
        dy(k) = dy(k) + U(k,j)*dy(j)
      end do
      dy(k) = y(k) - dy(k)
      dy(k) = dy(k)/U(k,k)
    end do
  end subroutine


  !> Get the Matrices L and U either from A or given
  !!
  !! \todo The Python preprocessor needs to include here the code
  subroutine GetLU(neq,y,t,beta,h,L,U)
    implicit none
    integer :: neq
    double precision :: t, beta, h
    double precision, dimension(neq) :: y
    double precision, dimension(neq,neq) :: L,U
    intent (in) :: neq, y, t, beta, h
    intent (out) :: L, U

    L = transpose(reshape( &
      (/ double precision :: 1.0, 0.0, 0.0, &
       0.4*beta*h/(0.04*beta*h+1.), 1.0, 0.0 , &
       0.0, -60000000.0*beta*h*y(2)/(-400.0*beta**2*h**2*y(3)/(0.04*beta*h + 1) &
       - beta*h*(-60000000.0*y(2) - 10000.0*y(3)) + 1), 1.0 /), &
       (/3,3/)))

    U = transpose(reshape( &
      (/ double precision :: 0.04*beta*h + 1, -10000.0*beta*h*y(3), -10000.0*beta*h*y(2) , & ! end row 1
       0, -400.0*beta**2*h**2*y(3)/(0.04*beta*h + 1) - beta*h*(-60000000.0*y(2) - 10000.0*y(3)) + 1, &
       -400.0*beta**2*h**2*y(2)/(0.04*beta*h + 1) + 10000.0*beta*h*y(2), &! end row 2
       0, 0, 60000000.0*beta*h*y(2)*(-400.0*beta**2*h**2*y(2)/(0.04*beta*h + 1) &
        + 10000.0*beta*h*y(2))/(-400.0*beta**2*h**2*y(3)/(0.04*beta*h + 1) -  &
        beta*h*(-60000000.0*y(2) - 10000.0*y(3)) + 1) + 1 /), & ! end row 3
       (/3,3/)))

  end subroutine

  !> Calculates the residuum
  !!
  !! Calculates
  !! \f[
  !!    r = y_n - \sum_{i=1}^{order} \alpha_i * y_{n-i} - h \beta_0 f(t_n, y_n).
  !! \f]
  !! For a converging system \f$ r \approx 0 \f$ should become true.
  subroutine CalcResiduum(nvector,neq,y,y_null,y_old,el,rhs_init,order,h,beta,res)
    implicit none
    integer :: neq,nvector,order
    double precision, dimension(neq) :: y, y_null, rhs, a, res, el, rhs_init
    double precision, dimension(neq,6) :: y_old
    double precision :: h, beta, t
    integer :: i, j
    double precision, dimension(7,6) :: coeff           !< coefficients of bdf-formula
    intent (in) :: y, y_old, nvector, neq, order, h, beta, el
    intent (out) :: res

    call CalcRHS(neq,y,rhs)

    res = 0
    do i = 0, order
      res = h*rhs - h*y_old(:,2) - el
    end do
  end subroutine

  !> Example case for the right-hand-side
  !!
  !! \todo needs to be deleted later!
  !! This is Robertson's example
  subroutine CalcRHS(neq,y,rhs)
    implicit none
    integer :: neq
    double precision, dimension(neq) :: y
    double precision, dimension(neq) :: rhs
    intent (in) :: y, neq
    intent (out) :: rhs

    rhs(1) = -0.04*y(1) + 1d4*y(2)*y(3)
    rhs(2) = 0.04*y(1) - 3d7*y(2)*y(2) - 1d4*y(2)*y(3)
    rhs(3) = 3d7*y(2)*y(2)
  end subroutine

  subroutine CalcError()
    implicit none

    !page 31 hindmarsh
    !ewt(i,n) = rtol(i)*abs(y(i,n-1) + atol(i))

    !sqrt(sum(d(:,n)/ewt(:,n)))

  end subroutine

  subroutine CalcOrder(order,iterator)
    implicit none
    integer :: order, maxorder, iterator
    intent(inout) :: order

    ! \todo{chose order by timestep estimation}
    maxorder = 5
!    if (iterator < 3) then
!      order = 1
!    else
      if (order < maxorder) then
        order = order + 1
      else
        order = maxorder
      end if
!    end if
  end subroutine

  subroutine CalcTimestep()
    implicit none

    ! empty routine

  end subroutine

end module bdf_method_mod
