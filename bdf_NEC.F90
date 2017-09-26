!> module that provides a cell-vectorized ODE-solver based on the BDF-scheme
module bdf_method_mod
  implicit none

  integer, parameter :: nvector=256                       !> vector length
  integer, parameter :: neq=3                             !> number equations
  integer, parameter :: maxorder=5                        !> maximum order
  integer            :: iterator                          !> iterator
  integer            :: order                             !> current order

  double precision :: rtol                                !> relative tolerance
  double precision :: atol                                !> absolute tolerance
  double precision :: t                                   !> current time
  double precision :: dt                                  !> current timestep
  double precision, pointer, dimension(:,:) :: y          !> current solution
                                                          !> array | often y_NS(:,:,0)
  double precision, pointer, dimension(:,:) :: en         !> correction vector
  double precision, pointer, dimension(:,:) :: en_old     !> old corr. vector
  double precision, pointer, dimension(:,:) :: den        !> corrector of en
  double precision, pointer, dimension(:,:) :: den_tmp    !> buffer array
  double precision, pointer, dimension(:,:) :: inv_weight_2 !> error weight
  double precision, pointer, dimension(:,:) :: rhs        !> right-hand-side
  double precision, pointer, dimension(:,:) :: res        !> residuum
  double precision, pointer, dimension(:,:) :: coeff      !> coefficient matrix
  double precision, pointer, dimension(:,:) :: tautable   !> solution table for
                                                          !> precomputation
  double precision, pointer, dimension(:,:,:) :: L,U      !> lower, upper triangular matrix
  double precision, pointer, dimension(:,:,:) :: y_NS     !> Nordsieck history array

  ! two different interfaces, whether it is called with y or y_NS
  interface WeightedNorm
    module procedure WeightedNorm1, WeightedNorm2
  end interface

contains


  !> Allocates all memory for bdf-solver
  subroutine InitBDF()
    implicit none
    integer          :: err, i, j

    ! allocate fields
    allocate( &
              y(nvector,neq), &
              en(nvector,neq), &
              en_old(nvector,neq), &
              den(nvector,neq), &
              den_tmp(nvector,neq), &
              inv_weight_2(nvector,neq), &
              rhs(nvector,neq), &
              res(nvector,neq), &
              coeff(0:6,6), &
              tautable(maxorder+1,0:maxorder+1), &
              L(nvector,neq,neq), &
              U(nvector,neq,neq), &
              y_NS(nvector,neq,0:maxorder+1), &
              STAT=err)
    if (err.ne.0) then
      print *, "Memory allocation error. BDF-Solver could not be intialized."
      stop
    end if

    ! bdf-coefficient-matrix with Nordsieck representation
    coeff = reshape((/  &
    1d0       ,1d0,0d0        ,0d0       ,0d0        ,0d0      ,0d0         , &
    2d0/3d0   ,1d0,1d0/3d0    ,0d0       ,0d0        ,0d0      ,0d0         , &
    6d0/11d0  ,1d0,6d0/11d0   ,1d0/11d0  ,0d0        ,0d0      ,0d0         , &
    12d0/25d0 ,1d0,7d0/10d0   ,1d0/5d0   ,1d0/50d0   ,0d0      ,0d0         , &
    60d0/137d0,1d0,225d0/274d0,85d0/274d0,15d0/274d0 ,1d0/274d0,0d0         , &
    20d0/49d0 ,1d0,58d0/63d0  ,5d0/12d0  ,25d0/252d0 ,1d0/84d0 ,1d0/1764d0 /),&
    (/7, 6/))

    ! precompute all tau and save them in a table (better performance)
    do i=1,maxorder+1
      do j=0,maxorder+1
        tautable(i,j) = tau(i,j,coeff)
      end do
    end do
  end subroutine


  !> deallocates all memory of the bdf-solver
  subroutine CloseBDF()
    implicit none
    integer :: err

    deallocate(y,en,en_old,den,den_tmp,inv_weight_2, &
               rhs,res,coeff,tautable,L,U,y_NS, &
               stat=err)

    if (err.ne.0) then
      print *, "Memory deallocation error. BDF-Solver was not properly closed."
      stop
    end if
  end subroutine


  !> Main subroutine to be called from outside
  subroutine SolveODE_BDF(t,dt,t_stop,y)
    implicit none
    double precision  :: t, t_stop, dt
    double precision, dimension(nvector,neq) :: y
    intent(in)        :: t_stop
    intent(inout)     :: y, t

    ! todo implement routine in order to calculate initial dt

    ! some initializations
    order       = 1
    y(:,:)      = y(:,:)
    y_NS(:,:,0) = y(:,:)
    iterator    = 0

    ! main solve - solve the linear system
    do while (t <= t_stop)
      ! solve the system
      call SolveLinearSystem(t,dt,y)

      iterator = iterator + 1
    end do
  end subroutine


  !> Solves a linear system for a given timestep
  subroutine SolveLinearSystem(t,dt,y)
    implicit none
    double precision, dimension(nvector,neq) :: y
    double precision :: conv_error, error, dt, t, dt_scale, conv_rate
    logical          :: success,reset
    integer          :: conv_iterator, lte_iterator
    intent(inout)    :: y, dt, t

    ! 1. initialization -------------------------------------------------------!
    ! use the LU matrices from the predictor value
    call GetLU(coeff(1,order),y,dt,L,U)

    ! Calculate initial right hand side
    call CalcRHS(y,rhs)
    y_NS(:,:,1) = dt*rhs(:,:)

    ! some initializations
    den           = 0.0d0
    conv_rate     = 0.7d0
    conv_error    = 1d20
    lte_iterator  = 0
    conv_iterator = 0

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(rtol,atol,y,inv_weight_2)


    ! 2. predictor ------------------------------------------------------------!
    success = .FALSE.
    predictor: do while(.not. success)
      reset = .FALSE.

      ! predictor step - set predicted y_NS, en=0.0
      call PredictSolution(y_NS,en)
      y(:,:) = y_NS(:,:,0)

      ! 3. corrector ----------------------------------------------------------!
      corrector: do while (conv_error > tautable(order,order)/ &
                                        (2d0*(order+2d0)))
        ! calculates residuum
        call CalcResiduum(y,dt,res)

        ! calculates the solution for dy with given residuum
        ! \todo  function needs to be enhanced for sparse matrices
        call SolveLU(L,U,res,den)

        ! add correction to solution vector
        en(:,:) = en(:,:) + den(:,:)
        y(:,:) = y_NS(:,:,0) + coeff(0,order)*en

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25
        ! times the step-size
        call CheckConvergence(conv_iterator,conv_rate,den, &
                              conv_error,inv_weight_2,reset)
        conv_iterator = conv_iterator + 1
        if (reset) then
          dt_scale = 0.25
          call ResetSystem(dt_scale,dt,y_NS)
          cycle predictor
        end if
      end do corrector

      ! local truncation error test:
      ! checks if solution is good enough and, similar to the convergence
      ! test, rerun from predictor with adjusted step size
      error = WeightedNorm(en,inv_weight_2)/ &
              tautable(order,order)
      if (error > 1d0) then
        dt_scale = 0.2
        call ResetSystem(dt_scale,dt,y_NS)
        lte_iterator = lte_iterator + 1
        cycle predictor
      else
        success = .TRUE.
      end if
    end do predictor

    ! after a successfull run rewrite the history array
    call UpdateNordsieck(en, y_NS)

    ! advance in time
    t = t + dt
!    print *, t, y(1,:)

    ! 4. step-size/order control ----------------------------------------------!
    ! calc. step size for current order+(0,-1,+1) -> use largest for next step
    call CalcStepSizeOrder(dt_scale,order)

    ! Adjust the Nordsieck history array with new step size & new order
    call SetStepSize(dt_scale,y_NS,dt)
    en_old = en
  end subroutine SolveLinearSystem


  !> Resets the whole system to the state where the predictor starts
  subroutine ResetSystem(dt_scale,dt,y_NS)
    implicit none
    integer          :: j,k,i,l
    double precision :: dt_scale, dt
    double precision, dimension(nvector,neq,0:6) :: y_NS
    intent(in)       :: dt_scale
    intent(inout)    :: y_NS, dt

    ! set the matrix back to old value
    do k = 0,order-1
      do j = order,k+1,-1
!cdir collapse
        do l = 1,neq
          do i = 1,nvector
            y_NS(i,l,j-1) = y_NS(i,l,j-1) - y_NS(i,l,j)
          end do
        end do
      end do
    end do

    call SetStepSize(dt_scale,y_NS,dt)
  end subroutine


  !> Calculates the step-sizes and uses according orders
  subroutine CalcStepSizeOrder(dt_scale,order)
    implicit none
    integer          :: dt_maxloc, order
    double precision :: dt_scale, dt_scale_up, dt_scale_down, dt_upper_limit
    intent(out)      :: dt_scale
    intent(inout)    :: order

    ! calculate all estimated step sizes
    ! todo err_weight_2 = err_weight*err_weight
    dt_scale_down = 1.d0/(1.3*(WeightedNorm(y_NS,inv_weight_2,order)/ &
                          tautable(order,order-1))**(1d0/order) + 1d-6)
    dt_scale      = 1.d0/(1.2*(WeightedNorm(en,inv_weight_2)/ &
                          tautable(order,order))**(1d0/(order+1d0)) + 1d-6)
    dt_scale_up   = 1.d0/(1.4*(WeightedNorm(en-en_old,inv_weight_2)/ &
                          tautable(order,order+1))**(1d0/(order+2d0)) + 1d-6)

    ! choose largest and search for location of largest
    dt_maxloc = maxloc((/ dt_scale_down, dt_scale, dt_scale_up /),DIM=1)
    dt_scale = maxval((/ dt_scale_down, dt_scale, dt_scale_up /))

    ! maximum increasement in step-size
    dt_upper_limit = 10d0
    dt_scale = min(dt_upper_limit,dt_scale)

    ! set new order
    if (dt_maxloc == 1) then
      if (order /= 1) then
        order = order - 1
      else
        ! do nothing
      end if
    else if (dt_maxloc == 2) then
      ! do nothing
    else if (dt_maxloc == 3) then
      if (order /= maxorder) then
        y_NS(:,:,order+1) = coeff(order,order)*en(:,:)/(order+1d0)
        order = order + 1
      else
        ! do nothing
      end if
    end if
  end subroutine


  !> Calculates the error for convergence
  subroutine CheckConvergence(iterator,conv_rate,den,error,inv_weight_2,reset)
    implicit none
    integer          :: iterator
    logical          :: reset
    double precision :: conv_rate, error_tmp, error
    double precision, dimension(nvector,neq) :: den, inv_weight_2
    intent(in)       :: den, inv_weight_2, iterator
    intent(inout)    :: conv_rate, error, reset

    error_tmp = WeightedNorm(den,inv_weight_2)
    conv_rate = max(0.2d0*conv_rate,error_tmp/error)
    error = error_tmp*min(1d0,1.5d0*conv_rate)

    if (iterator == 2 .and. error > 2d0) then
      reset = .TRUE.
    else if (iterator > 2 .and.  &
            error > tautable(order,order)/(2d0*(order + 2d0))) then
      reset = .TRUE.
    else if (iterator > 3) then
      reset = .TRUE.
    else
      reset = .FALSE.
    end if
  end subroutine CheckConvergence


  !> Calculates the weighted norm (two different interfaces for performance)
  pure function WeightedNorm1(en,inv_weight_2)
    implicit none
    integer :: i, j
    double precision :: WeightedNorm1
    double precision, dimension(nvector,neq) :: en, inv_weight_2
    intent(in)       :: en, inv_weight_2

    WeightedNorm1 = sqrt(sum(en(:,:)*en(:,:)*inv_weight_2(:,:))/ &
                         (neq*nvector))
  end function WeightedNorm1
  pure function WeightedNorm2(en,inv_weight_2,column)
    implicit none
    integer          :: column
    double precision :: WeightedNorm2
    double precision, dimension(nvector,neq,0:maxorder+1) :: en
    double precision, dimension(nvector,neq) :: inv_weight_2
    intent(in)       :: en, inv_weight_2, column

    WeightedNorm2 = sqrt(sum(en(:,:,column)*en(:,:,column)*inv_weight_2(:,:))/ &
                         (neq*nvector))
  end function WeightedNorm2


  !> Calculates the error-weight of y for given tolerances
  subroutine CalcErrorWeightInvSquared(rtol,atol,y,inv_weight_2)
    implicit none
    integer          :: i, j
    double precision :: rtol, atol
    double precision, dimension(nvector,neq) :: y, inv_weight_2
    intent(in)       :: rtol, atol, y
    intent(out)      :: inv_weight_2

    inv_weight_2(:,:) = 1./(rtol*abs(y(:,:)) + atol)
    inv_weight_2(:,:) = inv_weight_2(:,:)*inv_weight_2(:,:)
  end subroutine CalcErrorWeightInvSquared


  !> Sets the step size and scales the Nordsieck-array accordingly
  subroutine SetStepSize(dt_scale,y_NS,dt)
    implicit none
    integer          :: i,j,k
    double precision :: dt_scale,dtt_scale,dt
    double precision, dimension(nvector,neq,0:maxorder+1) :: y_NS
    intent(inout)    :: y_NS,dt
    intent(in)       :: dt_scale

    dtt_scale = 1.0
    do k = 1,order
      dtt_scale = dt_scale*dtt_scale
!cdir collapse
      do j = 1,neq
        do i = 1,nvector
          y_NS(i,j,k) = y_NS(i,j,k)*dtt_scale
        end do
      end do
    end do

    dt = dt*dt_scale
  end subroutine SetStepSize


  !> Sets the new Nordsieck array for a given correction vector
  subroutine UpdateNordsieck(en,y_NS)
    implicit none
    integer        :: i,j,k
    double precision, dimension(nvector,neq)                    :: en
    double precision, dimension(nvector,neq,0:maxorder+1)  :: y_NS
    intent(in)     :: en
    intent(inout)  :: y_NS

    do k = 0,order
!cdir collapse
      do j = 1,neq
        do i = 1,nvector
          y_NS(i,j,k) = y_NS(i,j,k) + en(i,j)*coeff(k,order)
        end do
      end do
    end do
  end subroutine UpdateNordsieck


  !> Predicts the solution
  subroutine PredictSolution(y_NS,en)
    implicit none
    integer         :: i,j,k,l
    double precision, dimension(nvector,neq,0:6) :: y_NS
    double precision, dimension(nvector,neq)     :: en
    intent (inout)  :: y_NS
    intent (out)    :: en

    !  loop effectively solves the Pascal Triangle without multiplications
    do l = 0,order-1
      do k = order,l+1,-1
!cdir collapse
        do j = 1,neq
          do i = 1,nvector
            y_NS(i,j,k-1) = y_NS(i,j,k) + y_NS(i,j,k-1)
          end do
        end do
      end do
    end do

    ! set the correction vector to zero
    en(:,:) = 0.0
  end subroutine PredictSolution


  !> Solve the System \f$ LU x=r \f$ with residuum r (res) and return x (den)
  subroutine SolveLU(L,U,res,den)
    implicit none
    integer        :: i,j,k
    double precision, dimension(nvector,neq,neq) :: L, U
    double precision, dimension(nvector,neq)          :: res, den
    intent(in)     :: L, U, res
    intent(out)    :: den

!cdir collapse
    do j=1,neq
      do i=1,nvector
        den_tmp(i,j) = 0.0
      end do
    end do
    ! lower system
    do k = 1,neq,1
      do j = 1,k-1
!cdir nodep
        do i = 1,nvector
          den_tmp(i,k) = den_tmp(i,k) + L(i,k,j)*den_tmp(i,j)
        end do
      end do
!cdir nodep
      do i = 1,nvector
        den_tmp(i,k) = res(i,k) - den_tmp(i,k)
        den_tmp(i,k) = den_tmp(i,k)/L(i,k,k)
      end do
    end do

    ! upper system
!cdir collapse
    do j=1,neq
      do i=1,nvector
        den(i,j) = 0.0
      end do
    end do
    do k = neq,1,-1
      do j = k+1,neq
!cdir nodep
        do i = 1,nvector
          den(i,k) = den(i,k) + U(i,k,j)*den(i,j)
        end do
      end do
!cdir nodep
      do i = 1,nvector
        den(i,k) = den_tmp(i,k) - den(i,k)
        den(i,k) = den(i,k)/U(i,k,k)
      end do
    end do
  end subroutine SolveLU


  !> Provides the Matrices L and U
  !!
  !! \todo The Python preprocessor needs to include the code here
  subroutine GetLU(beta,y,dt,L,U)
    implicit none
    integer          :: i
    double precision :: beta, dt
    double precision, dimension(nvector,neq) :: y
    double precision, dimension(nvector,neq,neq) :: L,U
    intent (in)      :: y, beta, dt
    intent (out)     :: L, U

!cdir nodep
    do i=1,nvector

#ODEVEC_L

#ODEVEC_U

    end do
  end subroutine GetLU


  !> Calculates the residuum
  subroutine CalcResiduum(y,dt,res)
    implicit none
    integer           :: i,j
    double precision  :: dt
    double precision, dimension(nvector,neq) :: y, res
    intent (in)       :: y, dt
    intent (out)      :: res

    call CalcRHS(y,rhs)

!cdir collapse
    do j=1,neq
      do i=1,nvector
        res(i,j) = dt*rhs(i,j) - y_NS(i,j,1) - en(i,j)
      end do
    end do
  end subroutine CalcResiduum


  !> Example case for the right-hand-side
  !!
  !! \todo needs to be done external later!
  !! This is for Robertson's example
  subroutine CalcRHS(y,rhs)
    implicit none
    integer         :: i
    double precision, dimension(nvector,neq) :: y
    double precision, dimension(nvector,neq) :: rhs
    intent (in)     :: y
    intent (out)    :: rhs

!cdir nodep
    do i=1,nvector
#ODEVEC_RHS
    end do
  end subroutine CalcRHS


  !> taufunction needed in several tests
  function tau(order,order2,coeff)
    implicit none
    integer :: i, order, order2
    double precision :: tau, faculty
    double precision, dimension(0:6,6) :: coeff
    intent(in) :: coeff, order

    faculty=1d0
!cdir shortloop
    do i=1,order
      faculty = faculty*i
    end do

    tau = (order2+1d0)/(faculty*coeff(order,order))

    return
  end function tau

end module bdf_method_mod
