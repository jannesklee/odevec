module bdf_method_mod
  implicit none

contains

  subroutine SolveODE_BDF(nvector,neq,y,t,t_stop,dt)
    implicit none
    integer :: neq,nvector, order, iterator
    double precision, dimension(nvector,neq) :: y, en
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision, dimension(0:6,6) :: coeff
    double precision :: t, t_stop, dt, start, finish
    intent(in) :: nvector, neq, t_stop
    intent(inout) :: y, t, dt

    ! bdf-matrix with Nordsieck representation
    coeff = reshape((/  &
          1d0        , 1d0 , 0d0         , 0d0        , 0d0         , 0d0       , 0d0         , &
          2d0/3d0    , 1d0 , 1d0/3d0     , 0d0        , 0d0         , 0d0       , 0d0         , &
          6d0/11d0   , 1d0 , 6d0/11d0    , 1d0/11d0   , 0d0         , 0d0       , 0d0         , &
          12d0/25d0  , 1d0 , 7d0/10d0    , 1d0/5d0    , 1d0/50d0    , 0d0       , 0d0         , &
          60d0/137d0 , 1d0 , 225d0/274d0 , 85d0/274d0 , 15d0/274d0  , 1d0/274d0 , 0d0         , &
          20d0/49d0  , 1d0 , 58d0/63d0   , 5d0/12d0   , 25d0/252d0  , 1d0/84d0  , 1d0/1764d0 /), &
          (/7, 6/))

    ! load initial condition into Nordsieck array
    order = 1
    y_NS(:,:,0) = y(:,:)
    en = 0.0

    ! initial conditions
    iterator = 0
    print *, t, y(100,:)

    call cpu_time(start)

    ! main solve - solve the linear system
    do while (t <= t_stop)

      ! solve the system
      call SolveLinearSystem(nvector,neq,order,coeff,en,y,y_NS,t,dt)

      iterator = iterator + 1
    end do

    call cpu_time(finish)
    print '("")'
    print '("Time = ", f10.8, " seconds.")', finish-start
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
  subroutine SolveLinearSystem(nvector,neq,order,coeff,en_old,y,y_NS,t,dt)
    implicit none
    integer :: neq, nvector, order
    double precision, dimension(nvector,neq) :: dy, res, y, rhs_init, den, err_weight
    double precision, dimension(nvector,neq) :: en_old, en
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision, dimension(nvector,neq,neq) :: L, U
    double precision, dimension(0:6,6) :: coeff
    double precision :: theta, sigma, conv_error, error, dt, t, beta, rtol, atol, dt_scale, conv_rate
    logical :: converged,success,reset
    integer :: conv_iterator, lte_iterator
    intent(inout) :: y_NS, y, dt, en_old, order, t
    intent(in) :: nvector, neq, coeff


    ! 1. initialization ---------------------------------------------------------------------------!
    ! general variables
    sigma = 0.01d0
    beta = 1.0d0
    theta = HUGE(1d0)
    rtol = 1d-10
    atol = 1d-20

    ! use the LU matrices from the predictor value
    call GetLU(nvector,neq,y,t,beta,dt,L,U)

    ! Calculate initial right hand side
    call CalcRHS(nvector,neq,t,y_NS(:,:,0),rhs_init)
    y_NS(:,:,1) = dt*rhs_init(:,:)

    ! some initializations
    den = 0.0d0
    conv_rate = 0.7d0
    lte_iterator = 0
    conv_iterator  = 0
    conv_error = 1d20

    ! important quantity to measure density
    call CalcErrorWeight(nvector,neq,rtol,atol,y,err_weight)

    ! 2. predictor --------------------------------------------------------------------------------!
    success = .FALSE.
    predictor: do while(.not. success)
      converged=.FALSE.
      reset=.FALSE.

      ! predictor step - set predicted y_NS, en=0.0
      call PredictSolution(nvector,neq,order,y_NS,en)
      y(:,:) = y_NS(:,:,0)
      ! 3. corrector ------------------------------------------------------------------------------!
      corrector: do while (conv_error > tau(order,order,coeff)/(2d0*(order+2d0)))
        ! calculates residuum
        call CalcResiduum(nvector,neq,y,y_NS,en,dt,res)

        ! calculates the solution for dy with given residuum
        ! \todo this function needs to be enhanced for sparse matrices
        call SolveLU(nvector,neq,L,U,res,den)

        ! add correction to solution vector
        en = en + den
        dy = coeff(0,order)*en
        y = y_NS(:,:,0) + dy

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25 times the step-size
        call CheckConvergence(nvector,neq,order,conv_iterator,coeff,conv_rate,den,conv_error,err_weight,converged,reset)
        conv_iterator = conv_iterator + 1
        if (reset) then
          dt_scale = 0.25
          call ResetSystem(nvector,neq,order,dt_scale,y_NS,dt)
          cycle predictor
        end if
      end do corrector

      ! local truncation error test:
      ! checks if solution is good enough and rewind everything if necessary
      error = WeightedNorm(nvector,neq,en,err_weight)/tau(order,order,coeff)
      if (error > 1d0) then
        dt_scale = 0.2
        call ResetSystem(nvector,neq,order,dt_scale,y_NS,dt)
        lte_iterator = lte_iterator + 1
        cycle predictor
      else
        success = .TRUE.
      end if
    end do predictor

    ! rewrite history array
    call UpdateNordsieck(nvector,neq, order, coeff, en, y_NS)

   ! write the result to the terminal or elsewhere
    t = t + dt
    print *, t, y(100,:)

    ! 4. sanity checks & step-size/order control --------------------------------------------------!
    ! calculate step size for current, upper & lower order and use largest one
    ! additionally it changes the order and enlarges the Nordsieck history array if necessary
    call CalcStepSizeOrder(nvector,neq,order,coeff,err_weight,en,en_old,y_NS,dt_scale)

    ! Adjust the Nordsieck history array with new step size & new order
    call SetStepSize(nvector,neq,order,dt_scale,y_NS,dt)
    en_old = en
  end subroutine SolveLinearSystem


  !>
  subroutine ResetSystem(nvector,neq,order,dt_scale,y_NS,dt)
    implicit none
    integer :: nvector,neq,order,j,k,i,l
    double precision :: dt_scale, dt
    double precision, dimension(nvector,neq,0:6) :: y_NS
    intent(in) :: nvector,neq, order, dt_scale
    intent(inout) :: y_NS, dt

    ! set the matrix back to old value
    do k = 0,order-1
      do j = order,k+1,-1
        do l = 1,neq
          do i = 1,nvector
            y_NS(i,l,j-1) = y_NS(i,l,j-1) - y_NS(i,l,j)
          end do
        end do
      end do
    end do

    call SetStepSize(nvector,neq,order,dt_scale,y_NS,dt)
  end subroutine


  !> Calculates the step-sizes and uses according orders
  subroutine CalcStepSizeOrder(nvector,neq,order,coeff,err_weight,en,en_old,y_NS,dt_scale)
    implicit none
    integer :: nvector, neq, order, dt_maxloc, maxorder
    double precision, dimension(0:6,6) :: coeff
    double precision, dimension(nvector,neq) :: en, err_weight, en_old
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision :: dt_scale, dt_scale_up, dt_scale_down, dt_upper_limit
    intent(in) :: nvector, neq, coeff, err_weight, en
    intent(inout) :: order
    intent(out) :: dt_scale

    ! calculate all estimated step sizes
    ! \todo check if order is correct here
    dt_scale_down = 1.d0/(1.3*(WeightedNorm(nvector,neq,y_NS(:,:,order),err_weight)/tau(order,order-1,coeff))**(1d0/order) + 1d-6)
    dt_scale = 1.d0/(1.2*(WeightedNorm(nvector,neq,en,err_weight)/tau(order,order,coeff))**(1d0/(order+1d0)) + 1d-6)
    dt_scale_up = 1.d0/(1.4*(WeightedNorm(nvector,neq,en-en_old,err_weight)/tau(order,order+1,coeff))**(1d0/(order+2d0)) + 1d-6)

    ! choose largest and search for location of largest
    dt_maxloc = maxloc((/ dt_scale_down, dt_scale, dt_scale_up /),DIM=1)
    dt_scale = maxval((/ dt_scale_down, dt_scale, dt_scale_up /))

    ! maximum increasement in step-size
    dt_upper_limit = 10d0
    dt_scale = min(dt_upper_limit,dt_scale)

    ! set new order
    maxorder = 5
    if (dt_maxloc == 1) then
      if (order /= 1) then
        order = order - 1
      else
        ! do nothing
      end if
    else if (dt_maxloc == 2) then
      ! do nothing
    else if (dt_maxloc == 3) then ! enlarge Nordsieck array by 1
      if (order /= maxorder) then
        y_NS(:,:,order+1) = coeff(order,order)*en(:,:)/(order+1d0)
        order = order + 1
      else
        ! do nothing
      end if
    end if
  end subroutine


  !> Calculates the error for convergence
  subroutine CheckConvergence(nvector,neq,order,iterator,coeff,conv_rate,den,error,weight,converged,reset)
    implicit none
    integer :: nvector, neq, order, iterator
    logical :: converged, reset
    double precision :: conv_rate, error_tmp, error
    double precision, dimension(nvector,neq) :: den, weight
    double precision, dimension(0:6,6) :: coeff
    intent(in) :: neq, den, weight, iterator, order,coeff
    intent(out) :: converged
    intent(inout) :: conv_rate, error, reset

    error_tmp = WeightedNorm(nvector,neq,den,weight)
    conv_rate = max(0.2d0*conv_rate,error_tmp/error)
    error = error_tmp*min(1d0,1.5d0*conv_rate)

    if (iterator == 2 .and. error > 2d0) then
      reset = .TRUE.
    else if (iterator > 2 .and. error > tau(order,order,coeff)/(2d0*(order + 2d0))) then
      reset = .TRUE.
    else if (iterator > 3) then
      reset = .TRUE.
    else
      reset = .FALSE.
    end if
  end subroutine CheckConvergence


  !> Weighted Norm
  function WeightedNorm(nvector,neq,en,weight)
    implicit none
    integer :: nvector,neq
    double precision, dimension(nvector,neq) :: en, weight
    double precision :: WeightedNorm
    intent(in) :: en, weight

    WeightedNorm = sqrt(sum(en*en*weight*weight)/(neq*nvector))
  end function WeightedNorm


  !>
  subroutine CalcErrorWeight(nvector,neq,rtol,atol,y,err_weight)
    implicit none
    integer :: neq, nvector, i, j
    double precision, dimension(nvector,neq) :: y,err_weight
    double precision :: rtol, atol
    intent(in) :: neq, rtol, atol, y
    intent(out) :: err_weight

    ! \todo this has to be different most probably for y
    do j=1,neq
      do i=1,nvector
        err_weight(i,j) = 1./(rtol*abs(y(i,j)) + atol)
      end do
    end do
  end subroutine CalcErrorWeight


  !> Sets the step size.
  !!
  !! When the step-size is changed also the corresponding history array needs to be adjusted. Since
  !! we are using the Nordsieck representation this is easily done by multiplying the track
  !! \f$r,rr,rrr,...\f$ to the Nordsieck array.
  subroutine SetStepSize(nvector,neq,order,dt_scale,y_NS,dt)
    implicit none
    integer :: nvector,neq,order,i,j,k
    double precision :: dt_scale,dtt_scale,dt
    double precision, dimension(nvector,neq,0:6) :: y_NS
    intent(inout) :: y_NS,dt
    intent(in) :: order,dt_scale,neq

    dtt_scale = 1.0
    do k = 1,order
      dtt_scale = dt_scale*dtt_scale
      do j = 1,neq
        do i = 1,nvector
          y_NS(i,j,k) = y_NS(i,j,k)*dtt_scale
        end do
      end do
    end do

    dt = dt*dt_scale
  end subroutine SetStepSize


  !>
  subroutine UpdateNordsieck(nvector,neq, order, coeff, en, y_NS)
    integer :: nvector,neq,i,j,k,order
    double precision, dimension(nvector,neq) :: en
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision, dimension(0:6,6) :: coeff
    intent(in) :: nvector,neq, order, coeff, en
    intent(inout) :: y_NS

    do k = 0,order
      do j = 1,neq
        do i = 1,nvector
          y_NS(i,j,k) = y_NS(i,j,k) + en(i,j)*coeff(k,order)
        end do
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
  subroutine PredictSolution(nvector,neq,order,y_NS,en)
    implicit none
    integer :: i,j,k,l,order,neq,nvector
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision, dimension(nvector,neq) :: en
    intent (in) :: neq,order,nvector
    intent (inout) :: y_NS
    intent (out) :: en

    ! this loop effectively solves the Pascal Triangle without multiplications
    do l = 0,order-1
      do k = order,l+1,-1
        do j = 1,neq
          do i = 1,nvector
            y_NS(i,j,k-1) = y_NS(i,j,k) + y_NS(i,j,k-1)
          end do
        end do
      end do
    end do

    ! set the correction vector to zero
    en = 0.0
  end subroutine PredictSolution


  !> Solve the System \f$ LUx=r \f$
  !!
  !! \todo this is until now for dense systems. it should be changed to sparse matrices
  subroutine SolveLU(nvector,neq,L,U,res,dy)
    implicit none
    integer :: neq,nvector
    double precision, dimension(nvector,neq,neq) :: L, U
    double precision, dimension(nvector,neq) :: dy, res
    double precision, dimension(nvector,neq) :: y
    intent(in) :: res, L, U, neq
    intent(out) :: dy
    integer :: i,j,k

    y = 0.0
    ! lower system
    do k = 1,neq,1
      do j = 1,k-1
        do i = 1,nvector
          y(i,k) = y(i,k) + L(i,k,j)*y(i,j)
        end do
      end do
      y(:,k) = res(:,k) - y(:,k)
      y(:,k) = y(:,k)/L(:,k,k)
    end do

    ! upper system
    dy = 0.0
    do k = neq,1,-1
      do j = k+1,neq
        do i = 1,nvector
          dy(i,k) = dy(i,k) + U(i,k,j)*dy(i,j)
        end do
      end do
      dy(:,k) = y(:,k) - dy(:,k)
      dy(:,k) = dy(:,k)/U(:,k,k)
    end do
  end subroutine SolveLU


  !> Get the Matrices L and U either from A or given
  !!
  !! \todo The Python preprocessor needs to include here the code
  subroutine GetLU(nvector,neq,y,t,beta,dt,L,U)
    implicit none
    integer :: nvector,neq, i
    double precision :: t, beta, dt
    double precision, dimension(nvector,neq) :: y
!    double precision, dimension(3,3) :: L_tmp,U_tmp
    double precision, dimension(nvector,neq,neq) :: L,U
    intent (in) :: neq, y, t, beta, dt
    intent (out) :: L, U

    !> \todo this is most probably not performant yet. Find a better solution
    do i=1,nvector
      L(i,1:3,1:3) = transpose(reshape( &
        (/ 1.0d0, 0.0d0, 0.0d0, &
         0.4d0*beta*dt/(0.04d0*beta*dt+1d0), 1.0d0, 0.0d0 , &
         0.0d0, -60000000.0d0*beta*dt*y(i,2)/(-400.0d0*beta*beta*dt*dt*y(i,3)/(0.04d0*beta*dt + 1d0) &
         - beta*dt*(-60000000.0d0*y(i,2) - 10000.0d0*y(i,3)) + 1d0), 1.0d0 /), &
         (/3,3/)))

      U(i,1:3,1:3) = transpose(reshape( &
        (/ 0.04d0*beta*dt + i, -10000.0*beta*dt*y(i,3), -10000.0*beta*dt*y(i,2) , & ! end row 1
         0d0, -400.0*beta**2*dt**2*y(i,3)/(0.04*beta*dt + 1) - beta*dt*(-60000000.0*y(i,2) - 10000.0*y(i,3)) + i, &
         -400.0*beta**2*dt**2*y(i,2)/(0.04*beta*dt + 1) + 10000.0*beta*dt*y(i,2), &! end row 2
         0d0, 0d0, 60000000.0*beta*dt*y(i,2)*(-400.0*beta**2*dt**2*y(i,2)/(0.04*beta*dt + 1) &
          + 10000.0*beta*dt*y(i,2))/(-400.0*beta**2*dt**2*y(i,3)/(0.04*beta*dt + 1) -  &
          beta*dt*(-60000000.0*y(i,2) - 10000.0*y(i,3)) + 1) + 1 /), & ! end row 3
         (/3,3/)))

!      L(i,:,:) = L_tmp(:,:)
!      U(i,:,:) = U_tmp(:,:)
    end do

  end subroutine GetLU


  !> Calculates the residuum
  !!
  !! Calculates
  !! \f[
  !!    r = y_n - \sum_{i=1}^{order} \alpha_i * y_{n-i} - h \beta_0 f(t_n, y_n).
  !! \f]
  !! For a converging system \f$ r \approx 0 \f$ should become true.
  subroutine CalcResiduum(nvector,neq,y,y_NS,en,dt,res)
    implicit none
    integer :: i,j,neq,nvector
    double precision, dimension(nvector,neq) :: y, rhs, res, en
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision :: dt, t
    intent (in) :: y, y_NS, nvector, neq, dt, en
    intent (out) :: res

    call CalcRHS(nvector,neq,t,y,rhs)

    do j=1,neq
      do i=1,nvector
        res(i,j) = dt*rhs(i,j) - y_NS(i,j,1) - en(i,j)
      end do
    end do
  end subroutine CalcResiduum


  !> Example case for the right-hand-side
  !!
  !! \todo needs to be done external later!
  !! This is Robertson's example
  subroutine CalcRHS(nvector,neq,t,y,rhs)
    implicit none
    integer :: i, neq, nvector
    double precision, dimension(nvector,neq) :: y
    double precision, dimension(nvector,neq) :: rhs
    double precision :: t
    intent (in) :: y, neq, t, nvector
    intent (out) :: rhs

    do i=1,nvector
      rhs(i,1) = -0.04d0*y(i,1) + 1d4*y(i,2)*y(i,3)
      rhs(i,2) = 0.04d0*y(i,1) - 3d7*y(i,2)*y(i,2) - 1d4*y(i,2)*y(i,3)
      rhs(i,3) = 3d7*y(i,2)*y(i,2)
    end do
  end subroutine CalcRHS


  !>
!  subroutine CalcOrder(neq,coeff,iterator,order,y_NS)
!    implicit none
!    integer :: order, maxorder, iterator, neq
!    double precision, dimension(0:6,6) :: coeff
!    double precision, dimension(neq,0:6) :: y_NS
!    double precision, dimension(neq) :: en
!    intent(in) :: neq, coeff, iterator
!    intent(inout) :: order
!
!    maxorder = 5
!    if (order < maxorder) then
!      y_NS(:,order+1) = coeff(order,order)*en(:)/(order+1d0)
!      order = order + 1
!    else
!      order = maxorder
!    end if
!  end subroutine CalcOrder


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
