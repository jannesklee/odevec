module bdf_method_mod
  implicit none
contains

  subroutine SolveODE_BDF(neq, nvector, y, t, t_stop, dt)
    implicit none
    integer :: neq,nvector, order, iterator
    double precision, dimension(neq) :: y, rhs, en
    double precision, dimension(neq,6) :: y_NS
    double precision, dimension(7,6) :: coeff
    double precision :: t, t_old, t_stop, dt, r, rtol, atol
    double precision, dimension(neq) :: dy, res
    intent(inout) :: y

    ! bdf-matrix with Nordsieck representation
    coeff = reshape((/ double precision :: &
          1.       , 1. , 0.        , 0.       , 0.        , 0.      , 0.       , &
          2./3.    , 1. , 1./3.     , 0.       , 0.        , 0.      , 0.       , &
          6./11.   , 1. , 6./11.    , 1./11.   , 0.        , 0.      , 0.       , &
          12./25.  , 1. , 7./10.    , 1./5.    , 1./50.    , 0.      , 0.       , &
          60./137. , 1. , 225./274. , 85./274. , 15./274. , 1./274. , 0.       , &
          20./49.  , 1. , 58./63.   , 5./12.   , 25./252. , 1./84.  , 1./1764. /), &
          shape(coeff))

    ! load initial condition into Nordsieck array
    order = 1
    y_NS(:,1) = y


    ! initial conditions
    iterator = 0
    print *, t, y
    t = t + dt

    ! main solve - solve the linear system
    do while (t <= t_stop)

      ! solve the system
      call SolveLinearSystem(nvector,neq,order,coeff,y,y_NS,dt)

      print *, t, y

      t = t + dt
      iterator = iterator + 1
    end do
  end subroutine


  !> Solves a linear system for a given timestep
  !!
  !! Does the following:
  !! 1. initialization
  !!    - especially getting LU (needs to be prepared by a script or inserted by hand)
  !! 2. predict a solution
  !! 3. correct it (while loop)
  !!    - calculate change from residuum
  !!    - loop until convergence is reached
  !! 4. aftermath
  !!    - update Nordsieck array
  !!    - checking for sanity
  !!    - checking step-size & order
  !!      - change appropriately
  subroutine SolveLinearSystem(nvector,neq,order,coeff,y,y_NS,dt)
    implicit none
    integer :: neq, nvector, order
    double precision, dimension(neq) :: dy, res, y, rhs, dy_NS, y_null, en, rhs_init, den, y_safe, err_weight
    double precision, dimension(neq,6) :: y_NS
    double precision, dimension(neq,neq) :: L, U
    double precision, dimension(7,6) :: coeff
    double precision :: theta, sigma, tol, error, absdy, dt, t, beta, ss_ratio, rtol, atol, error_old
    integer :: i, j, iterator
    intent(inout) :: y_NS, y, dt
    intent(in) :: nvector, neq


    ! 1. initialization ---------------------------------------------------------------------------!
    ! general variables
    sigma = 0.01d0
    beta = 1.0d0
    theta = HUGE(1d0)
    rtol = 1d-4
    atol = 1d-20

    ! initialize solution arrays
    y_null = y

    ! use the LU matrices from the predictor value
    call GetLU(neq,y,t,beta,dt,L,U)

    ! Calculate initial right hand side
    call CalcRHS(neq,y_NS(:,1),rhs_init)
    y_NS(:,2) = dt*rhs_init

    ! some initializations
    iterator = 0
    dy_NS = 1.0
    den = 0.0

    call CalcErrorWeight(neq,rtol,atol,y,err_weight)


    ! 2. predictor --------------------------------------------------------------------------------!
    ! predictor step
    call PredictSolution(neq,order,y_NS,en)

    ! 3. corrector --------------------------------------------------------------------------------!
    iterator  = 0
    error = 0.7d0

    do while (error > tau(order,coeff)/(2d0*(order + 2d0)))
      ! calculates residuum
      call CalcResiduum(nvector,neq,y,y_null,y_NS,en,rhs_init,order,dt,beta,res)

      ! calculates the solution for dy with given residuum
      ! \todo this function needs to be enhanced for sparse matrices
      call SolveLU(neq,L,U,res,den)

      ! add correction to solution vector
      en = en + den
      dy = coeff(1,order)*en
      y = y_null + dy

      ! convergence tests
      call CheckConvergence(neq,res,err_weight,error)

      iterator = iterator + 1
      error_old = error
    end do

    ! 4. sanity checks/step size/order ------------------------------------------------------------!
    ! check error
    theta = WeightedNorm(neq,en,err_weight)/tau(order,coeff)
    ss_ratio = (1d0/theta)**(1d0/(order+1d0))
    print *, ss_ratio

    ! rewrite history array
    call UpdateNordsieck(neq, order, coeff, en, y_NS)
    ! if necessary calculate new order and ad column to history array
    call CalcOrder(neq,coeff,order,y_NS)
    ! order change must be done before, because step size changes
    call SetStepSize(neq,order,ss_ratio,y_NS,dt)
  end subroutine SolveLinearSystem


  !> Calculates the error for convergence
  subroutine CheckConvergence(neq,res,weight,error)
    implicit none
    integer :: neq
    double precision :: conv_rate, border, error_tmp, error
    double precision, dimension(neq) :: res, weight, y
    intent(inout) :: error
    intent(in) :: neq, res, weight

    error_tmp = WeightedNorm(neq,res,weight)

    conv_rate = max(0.2d0*conv_rate,error_tmp/error)

    error = error_tmp*min(1d0,1.5d0*conv_rate)

  end subroutine CheckConvergence


  !>
  function WeightedNorm(neq,en,weight)
    implicit none
    integer :: neq
    double precision, dimension(neq) :: en, weight
    double precision :: WeightedNorm

    WeightedNorm = sqrt(sum(en*weight*en*weight)/neq)
  end function WeightedNorm

  !>
  subroutine CalcErrorWeight(neq,rtol,atol,y,err_weight)
    implicit none
    integer :: neq
    double precision, dimension(neq) :: y,err_weight
    double precision :: rtol, atol
    intent(in) :: neq, rtol, atol, y
    intent(out) :: err_weight

    err_weight = 1./(rtol*abs(y) + atol)
  end subroutine CalcErrorWeight


  !> Sets the step size.
  !!
  !! When the step-size is changed also the corresponding history array needs to be adjusted.
  !! Since we are using the Nordsieck representation this is easily done by multiplying the
  !! track \f$r,rr,rrr,...\f$ to the Nordsieck array.
  subroutine SetStepSize(neq,order,r,y_NS,dt)
    implicit none
    integer :: neq,order,j
    double precision :: r,rr,dt
    double precision, dimension(neq,6) :: y_NS
    intent(inout) :: y_NS,dt
    intent(in) :: order,r,neq

    rr = 1.0
    do j = 1,order
      rr = r*rr
      y_NS(:,j) = y_NS(:,j)*rr
    end do

    dt = dt*r
  end subroutine SetStepSize



  subroutine UpdateNordsieck(neq, order, coeff, en, y_NS)
    integer :: neq,i,j, order
    double precision, dimension(neq) :: en
    double precision, dimension(neq,6) :: y_NS
    double precision, dimension(7,6) :: coeff
    double precision :: h, one, et

    do j = 1,order+1
      do i = 1,neq
        y_NS(i,j) = y_NS(i,j) + en(i)*coeff(order,j)
      end do
    end do
  end subroutine UpdateNordsieck


  !> predicts the solution
  !!
  !! This subroutine predicts the solution in Nordsieck representation by
  !! \f[
  !!    \mathbf{z}_n^{[0]} &= \mathbf{z}_{n-1} \mathbf{A} \\
  !!    e_n^{0} &= 0.
  !! \f]
  subroutine PredictSolution(neq,order,y_NS,en)
    implicit none
    integer :: i,j,k,order,neq
    double precision, dimension(neq,6) :: y_NS
    double precision, dimension(neq) :: en
    intent (in) :: neq,order
    intent (inout) :: y_NS
    intent (out) :: en

    ! this loop effectiveny solves the Pascal Triangle without multiplications
    do k = 0,order-1
      do j = order,k-1
        y_NS(:,j-1) = y_NS(:,j) + y_NS(:,j-1)
      end do
    end do
    en = 0.0
  end subroutine PredictSolution


  !> Solve the System \f$ LUx=r \f$
  !!
  !! \todo this is until now for dense systems. it should be changed to sparse matrices
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
  end subroutine SolveLU


  !> Get the Matrices L and U either from A or given
  !!
  !! \todo The Python preprocessor needs to include here the code
  subroutine GetLU(neq,y,t,beta,dt,L,U)
    implicit none
    integer :: neq
    double precision :: t, beta, dt
    double precision, dimension(neq) :: y
    double precision, dimension(neq,neq) :: L,U
    intent (in) :: neq, y, t, beta, dt
    intent (out) :: L, U

    L = transpose(reshape( &
      (/ double precision :: 1.0, 0.0, 0.0, &
       0.4*beta*dt/(0.04*beta*dt+1.), 1.0, 0.0 , &
       0.0, -60000000.0*beta*dt*y(2)/(-400.0*beta**2*dt**2*y(3)/(0.04*beta*dt + 1) &
       - beta*dt*(-60000000.0*y(2) - 10000.0*y(3)) + 1), 1.0 /), &
       (/3,3/)))

    U = transpose(reshape( &
      (/ double precision :: 0.04*beta*dt + 1, -10000.0*beta*dt*y(3), -10000.0*beta*dt*y(2) , & ! end row 1
       0, -400.0*beta**2*dt**2*y(3)/(0.04*beta*dt + 1) - beta*dt*(-60000000.0*y(2) - 10000.0*y(3)) + 1, &
       -400.0*beta**2*dt**2*y(2)/(0.04*beta*dt + 1) + 10000.0*beta*dt*y(2), &! end row 2
       0, 0, 60000000.0*beta*dt*y(2)*(-400.0*beta**2*dt**2*y(2)/(0.04*beta*dt + 1) &
        + 10000.0*beta*dt*y(2))/(-400.0*beta**2*dt**2*y(3)/(0.04*beta*dt + 1) -  &
        beta*dt*(-60000000.0*y(2) - 10000.0*y(3)) + 1) + 1 /), & ! end row 3
       (/3,3/)))

  end subroutine GetLU


  !> Calculates the residuum
  !!
  !! Calculates
  !! \f[
  !!    r = y_n - \sum_{i=1}^{order} \alpha_i * y_{n-i} - h \beta_0 f(t_n, y_n).
  !! \f]
  !! For a converging system \f$ r \approx 0 \f$ should become true.
  subroutine CalcResiduum(nvector,neq,y,y_null,y_NS,en,rhs_init,order,dt,beta,res)
    implicit none
    integer :: neq,nvector,order
    double precision, dimension(neq) :: y, y_null, rhs, a, res, en, rhs_init
    double precision, dimension(neq,6) :: y_NS
    double precision :: dt, beta, t
    integer :: i, j
    double precision, dimension(7,6) :: coeff           !< coefficients of bdf-formula
    intent (in) :: y, y_NS, nvector, neq, order, dt, beta, en
    intent (out) :: res

    call CalcRHS(neq,y,rhs)

    do i = 0, order
      res = dt*rhs - dt*y_NS(:,2) - en
    end do
  end subroutine CalcResiduum


  !> Example case for the right-hand-side
  !!
  !! \todo needs to be deneted later!
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
  end subroutine CalcRHS


  !>
  subroutine CalcOrder(neq,coeff,order,y_NS)
    implicit none
    integer :: i
    integer :: order, maxorder, iterator, neq
    double precision, dimension(7,6) :: coeff
    double precision, dimension(neq,6) :: y_NS
    double precision, dimension(neq) :: en
    intent(inout) :: order

    ! \todo{chose order by timestep estimation}
    maxorder = 5
!    if (iterator < 3) then
!      order = 1
!    ense
      if (order < maxorder) then
        order = order + 1
      else
        order = maxorder
      end if

      do i=1,neq
        y_NS(i,order+1) = coeff(i,order)*en(i)/(order+1d0)
      end do
  end subroutine CalcOrder


  !>
  function tau(order,coeff)
    implicit none
    integer :: i, order
    double precision :: tau, faculty
    double precision, dimension(7,6) :: coeff
    intent(in) :: coeff, order

    faculty=1d0
    do i=1,order
      faculty = faculty*i
    end do

    tau = (order+1d0)/(faculty*coeff(order,order))

    return
  end function tau

end module bdf_method_mod
