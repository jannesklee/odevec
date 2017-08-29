module bdf_method_mod
  implicit none
contains

  subroutine SolveODE_BDF(neq, nvector, y, t, t_stop, dt)
    implicit none
    integer :: neq,nvector, order, iterator
    double precision, dimension(neq) :: y, rhs, en
    double precision, dimension(neq,0:6) :: y_NS
    double precision, dimension(0:6,6) :: coeff
    double precision :: t, t_old, t_stop, dt, r, rtol, atol, start, finish
    double precision, dimension(neq) :: dy, res
    intent(inout) :: y

    ! bdf-matrix with Nordsieck representation
    coeff = reshape((/ double precision :: &
          1.       , 1. , 0.        , 0.       , 0.        , 0.      , 0.         , &
          2./3.    , 1. , 1./3.     , 0.       , 0.        , 0.      , 0.         , &
          6./11.   , 1. , 6./11.    , 1./11.   , 0.        , 0.      , 0.         , &
          12./25.  , 1. , 7./10.    , 1./5.    , 1./50.    , 0.      , 0.         , &
          60./137. , 1. , 225./274. , 85./274. , 15./274.  , 1./274. , 0.         , &
          20./49.  , 1. , 58./63.   , 5./12.   , 25./252.  , 1./84.  , 1./1764. /), &
          (/7, 6/))

    ! load initial condition into Nordsieck array
    order = 1
    y_NS(:,0) = y
    en = 0.0

    ! initial conditions
    iterator = 0
    print *, t, y

    call cpu_time(start)

    ! main solve - solve the linear system
    do while (t <= t_stop)

      ! solve the system
      call SolveLinearSystem(nvector,neq,order,coeff,en,y,y_NS,dt)

      t = t + dt

      print *, t, y

      call CalcOrder(neq,coeff,iterator,order,y_NS)
      iterator = iterator + 1
    end do

    call cpu_time(finish)
    print '("")'
    print '("Time = ",f6.3," seconds.")',finish-start
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
  subroutine SolveLinearSystem(nvector,neq,order,coeff,en_old,y,y_NS,dt)
    implicit none
    integer :: neq, nvector, order
    double precision, dimension(neq) :: dy, res, y, rhs, y_null, rhs_init, den, y_safe, err_weight
    double precision, dimension(neq) :: en_old, en
    double precision, dimension(neq,0:6) :: y_NS
    double precision, dimension(neq,neq) :: L, U
    double precision, dimension(0:6,6) :: coeff
    double precision :: theta, sigma, tol, error, absdy, dt, t, beta, rtol, atol, dt_scale, conv_rate
    logical :: convergence
    integer :: i, j, iterator, maxiter
    intent(inout) :: y_NS, y, dt, en_old, order
    intent(in) :: nvector, neq, coeff


    ! 1. initialization ---------------------------------------------------------------------------!
    ! general variables
    sigma = 0.01d0
    beta = 1.0d0
    theta = HUGE(1d0)
    rtol = 1d-4
    atol = 1d-9
    maxiter = 5

    ! initialize solution arrays
    y_null = y

    ! use the LU matrices from the predictor value
    call GetLU(neq,y,t,beta,dt,L,U)

    ! Calculate initial right hand side
    call CalcRHS(neq,y_NS(:,0),rhs_init)
    y_NS(:,1) = dt*rhs_init

    ! some initializations
    iterator = 0
    den = 0.0

    ! important quantity to measure density
    call CalcErrorWeight(neq,rtol,atol,y,err_weight)

    ! 2. predictor --------------------------------------------------------------------------------!
    ! predictor step - set predicted y_NS, en=0.0
    call PredictSolution(neq,order,y_NS,en)

    ! 3. corrector --------------------------------------------------------------------------------!
    iterator  = 0
    error = HUGE(1d0)

    do while (error > 1d0)!  .AND. iterator < maxiter)
!    do while (iterator < maxiter)
      ! calculates residuum
      call CalcResiduum(nvector,neq,y,y_NS,en,rhs_init,order,dt,beta,res)

      ! calculates the solution for dy with given residuum
      ! \todo this function needs to be enhanced for sparse matrices
      call SolveLU(neq,L,U,res,den)

      ! add correction to solution vector
      en = en + den
      dy = coeff(0,order)*en
      y = y_NS(:,0) + dy

!      print *, en

      ! convergence tests - if failed run again with smaller step size
!      call CheckConvergence(neq,order,coeff,conv_rate,res,err_weight,convergence)

      ! local trucation error this is used for
      error = WeightedNorm(neq,en,err_weight)/tau(order,order,coeff)

      iterator = iterator + 1
    end do

    ! rewrite history array
    call UpdateNordsieck(neq, order, coeff, en, y_NS)

    ! 4. sanity checks & step-size/order control --------------------------------------------------!
    ! calculate step size for current, upper & lower order and use largest one
    ! additionally it changes the order and enlarges the Nordsieck history array if necessary
    !call CalcStepSizeOrder(neq,order,coeff,err_weight,en,en_old,y_NS,dt_scale)

    ! Adjust the Nordsieck history array with new step size & new order
    !call SetStepSize(neq,order,dt_scale,y_NS,dt)
    en_old = en
  end subroutine SolveLinearSystem


  !> Calculates the step-sizes and uses according orders
  subroutine CalcStepSizeOrder(neq,order,coeff,err_weight,en,en_old,y_NS,dt_scale)
    implicit none
    integer :: neq, order, i, dt_maxloc, maxorder
    double precision, dimension(0:6,6) :: coeff
    double precision, dimension(neq) :: en, err_weight, en_old
    double precision, dimension(neq,0:6) :: y_NS
    double precision :: dt_scale, dt_scale_up, dt_scale_down
    intent(in) :: neq, coeff, err_weight, en
    intent(inout) :: order
    intent(out) :: dt_scale

    ! calculate all estimated step sizes
    ! \todo check if order is correct here
    dt_scale_down = (tau(order,order-1,coeff)/WeightedNorm(neq,y_NS(:,order),err_weight))**(1d0/order)
    dt_scale = (tau(order,order,coeff)/WeightedNorm(neq,en,err_weight))**(1d0/(order+1d0))
    dt_scale_up = (tau(order,order+1,coeff)/WeightedNorm(neq,en-en_old,err_weight))**(1d0/(order+2d0))

    print *, dt_scale_down, dt_scale, dt_scale_up

    ! choose largest and search for location of largest
    dt_maxloc = maxloc((/ dt_scale_down, dt_scale, dt_scale_up /),DIM=1)
    dt_scale = maxval((/ dt_scale_down, dt_scale, dt_scale_up /))

    ! set new order
    maxorder = 5
    if (order == maxorder .OR. order == 1) then
      ! do nothing
    else
      if (dt_maxloc == 1) then
        order = order - 1
      else if (dt_maxloc == 2) then
        ! do nothing
      else if (dt_maxloc == 3) then ! enlarge Nordsieck array by 1
        y_NS(:,order+1) = coeff(order,order)*en(:)/(order+1d0)
        order = order + 1
      end if
    end if
  end subroutine

  !> Calculates the error for convergence
  subroutine CheckConvergence(neq,order,coeff,conv_rate,res,weight,convergence)
    implicit none
    integer :: neq, order
    logical :: convergence
    double precision :: conv_rate, border, error_tmp, error, init_val
    double precision, dimension(neq) :: res, weight, y
    double precision, dimension(0:6,6) :: coeff
    intent(in) :: neq, res, weight
    intent(out) :: convergence
    intent(inout) :: conv_rate

    error_tmp = WeightedNorm(neq,res,weight)

    conv_rate = max(0.2d0*conv_rate,error_tmp/error)

    error = error_tmp*min(1d0,1.5d0*conv_rate)

    if (error < tau(order,order,coeff)/(2d0*(order + 2d0))) then
      convergence = .TRUE.
    else
      convergence = .FALSE.
    end if
  end subroutine CheckConvergence


  !>
  function WeightedNorm(neq,en,weight)
    implicit none
    integer :: neq
    double precision, dimension(neq) :: en, weight
    double precision :: WeightedNorm

    WeightedNorm = sqrt(sum(en*en*weight*weight)/neq)
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
  !! When the step-size is changed also the corresponding history array needs to be adjusted. Since
  !! we are using the Nordsieck representation this is easily done by multiplying the track
  !! \f$r,rr,rrr,...\f$ to the Nordsieck array.
  subroutine SetStepSize(neq,order,dt_scale,y_NS,dt)
    implicit none
    integer :: neq,order,j
    double precision :: dt_scale,dtt_scale,dt
    double precision, dimension(neq,0:6) :: y_NS
    intent(inout) :: y_NS,dt
    intent(in) :: order,dt_scale,neq

    dtt_scale = 1.0
    do j = 0,order
      dtt_scale = dt_scale*dtt_scale
      y_NS(:,j) = y_NS(:,j)*dtt_scale
    end do

    dt = dt*dt_scale
  end subroutine SetStepSize


  !>
  subroutine UpdateNordsieck(neq, order, coeff, en, y_NS)
    integer :: neq,i,j, order
    double precision, dimension(neq) :: en
    double precision, dimension(neq,0:6) :: y_NS
    double precision, dimension(0:6,6) :: coeff
    double precision :: h, one, et

    do j = 0,order
      do i = 1,neq
        y_NS(i,j) = y_NS(i,j) + en(i)*coeff(j,order)
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
    double precision, dimension(neq,0:6) :: y_NS
    double precision, dimension(neq) :: en
    intent (in) :: neq,order
    intent (inout) :: y_NS
    intent (out) :: en

    ! this loop effectively solves the Pascal Triangle without multiplications
    do k = 0,order-1
      do j = order,k+1,-1
        y_NS(:,j-1) = y_NS(:,j) + y_NS(:,j-1)
      end do
    end do

    ! set the correction vector to zero
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
  subroutine CalcResiduum(nvector,neq,y,y_NS,en,rhs_init,order,dt,beta,res)
    implicit none
    integer :: neq,nvector,order
    double precision, dimension(neq) :: y, y_null, rhs, a, res, en, rhs_init
    double precision, dimension(neq,0:6) :: y_NS
    double precision :: dt, beta, t
    integer :: i, j
    double precision, dimension(0:6,6) :: coeff           !< coefficients of bdf-formula
    intent (in) :: y, y_NS, nvector, neq, order, dt, beta, en
    intent (out) :: res

    call CalcRHS(neq,y,rhs)

!    res = dt*rhs - dt*y_NS(:,2) - en
    res = dt*rhs - y_NS(:,1) - en
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
  subroutine CalcOrder(neq,coeff,iterator,order,y_NS)
    implicit none
    integer :: i
    integer :: order, maxorder, iterator, neq
    double precision, dimension(0:6,6) :: coeff
    double precision, dimension(neq,0:6) :: y_NS
    double precision, dimension(neq) :: en
    intent(in) :: neq, coeff, iterator
    intent(inout) :: order

    maxorder = 5
    if (order < maxorder) then
      y_NS(:,order+1) = coeff(order,order)*en(:)/(order+1d0)
      order = order + 1
    else
      order = maxorder
    end if
  end subroutine CalcOrder


  !>
  function tau(order,order2,coeff)
    implicit none
    integer :: i, order, order2
    double precision :: tau, faculty
    double precision, dimension(0:6,6) :: coeff
    intent(in) :: coeff, order

    faculty=1d0
    do i=1,order
      faculty = faculty*i
    end do

    tau = (order2+1d0)/(faculty*coeff(order,order))

    return
  end function tau

end module bdf_method_mod
