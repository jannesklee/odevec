module bdf_method_mod
  implicit none

  type bdf_type
    integer :: nvector                                      !> vector length
    integer :: neq                                          !> number equations
    integer :: order                                        !> current order
    integer :: maxorder=5                                   !> maximum order (smaller value possible)
    integer :: iterator                                     !> iterator
    double precision :: t                                   !> current time
    double precision :: dt                                  !> current timestep
    double precision :: rtol                                !> relative tolerance
    double precision :: atol                                !> absolute tolerance
    double precision, pointer, dimension(:,:) :: y          !> current solution array | often y_NS(:,:,0)
    double precision, pointer, dimension(:,:) :: en         !> correction vector
    double precision, pointer, dimension(:,:) :: en_old     !> old correction vector
    double precision, pointer, dimension(:,:) :: den        !> corrector to correction vector
    double precision, pointer, dimension(:,:) :: den_tmp    !> buffer array
    double precision, pointer, dimension(:,:) :: err_weight !> error weight
    double precision, pointer, dimension(:,:) :: rhs        !> right-hand-side
    double precision, pointer, dimension(:,:) :: res        !> residuum in linear equation
    double precision, pointer, dimension(:,:) :: coeff      !> coefficient matrix
    double precision, pointer, dimension(:,:) :: tautable   !> solution table of an often evaluated function
    double precision, pointer, dimension(:,:,:) :: L,U      !> lower, upper triangular matrix
    double precision, pointer, dimension(:,:,:) :: y_NS     !> Nordsieck history array
  end type bdf_type


contains


  !> Allocates all memory for bdf-solver
  subroutine InitBDF(this,nvector,neq,rtol,atol)
    implicit none
    type(bdf_type)   :: this
    integer          :: nvector, neq
    integer          :: err, i, j
    double precision :: rtol, atol
    intent(in)       :: nvector, neq
    intent(out)      :: this

    this%nvector = nvector
    this%neq     = neq
    this%rtol    = rtol
    this%atol    = atol

    ! allocate fields
    allocate( &
              this%y(this%nvector,this%neq), &
              this%en(this%nvector,this%neq), &
              this%en_old(this%nvector,this%neq), &
              this%den(this%nvector,this%neq), &
              this%den_tmp(this%nvector,this%neq), &
              this%err_weight(this%nvector,this%neq), &
              this%rhs(this%nvector,this%neq), &
              this%res(this%nvector,this%neq), &
              this%coeff(0:6,6), &
              this%tautable(this%maxorder,0:this%maxorder), &
              this%L(this%nvector,this%neq,this%neq), &
              this%U(this%nvector,this%neq,this%neq), &
              this%y_NS(this%nvector,this%neq,0:this%maxorder+1), &
              STAT=err)
    if (err.ne.0) then
      print *, "Memory allocation error. BDF-Solver could not be intialized."
      stop
    end if

    ! bdf-coefficient-matrix in Nordsieck representation
    this%coeff = reshape((/  &
          1d0        , 1d0 , 0d0         , 0d0        , 0d0         , 0d0       , 0d0         , &
          2d0/3d0    , 1d0 , 1d0/3d0     , 0d0        , 0d0         , 0d0       , 0d0         , &
          6d0/11d0   , 1d0 , 6d0/11d0    , 1d0/11d0   , 0d0         , 0d0       , 0d0         , &
          12d0/25d0  , 1d0 , 7d0/10d0    , 1d0/5d0    , 1d0/50d0    , 0d0       , 0d0         , &
          60d0/137d0 , 1d0 , 225d0/274d0 , 85d0/274d0 , 15d0/274d0  , 1d0/274d0 , 0d0         , &
          20d0/49d0  , 1d0 , 58d0/63d0   , 5d0/12d0   , 25d0/252d0  , 1d0/84d0  , 1d0/1764d0 /), &
          (/7, 6/))

    ! precompute all tau and save them in a table (better performance)
    do i=1,this%maxorder
      do j=0,this%maxorder
        this%tautable(i,j) = tau(i,j,this%coeff)
      end do
    end do
  end subroutine


  !> deallocates all memory of the bdf-solver
  subroutine CloseBDF(this)
    implicit none
    type(bdf_type) :: this
    integer :: err

    deallocate(this%y,this%en,this%en_old,this%den,this%den_tmp,this%err_weight, &
               this%rhs,this%res,this%coeff,this%tautable,this%L,this%U,this%y_NS, &
               stat=err)

    if (err.ne.0) then
      print *, "Memory deallocation error. BDF-Solver was not properly closed."
      stop
    end if
  end subroutine


  !> Main subroutine to be called from outside
  subroutine SolveODE_BDF(this,t,dt,t_stop,y)
    implicit none
    type(bdf_type)    :: this
    double precision  :: t, t_stop, dt, rtol, atol
    double precision, dimension(this%nvector,this%neq) :: y
!    double precision  :: finish, start
    intent(in)        :: t_stop
    intent(inout)     :: y, t


    ! todo implement routine in order to calculate initial dt
    ! todo y needs to be either always run on this%y or it needs to be passed all the time

    ! some initializations
    this%order       = 1
    this%y(:,:)      = y(:,:)
    this%y_NS(:,:,0) = this%y(:,:)
    this%iterator    = 0

    ! initial conditions
    print *, t, y(1,:)

!    call cpu_time(start)

    ! main solve - solve the linear system
    do while (t <= t_stop)
      ! solve the system
      call SolveLinearSystem(this,t,dt,y)

      this%iterator = this%iterator + 1
    end do

!    call cpu_time(finish)
!    print '("")'
!    print '("Time = ", f6.3, " seconds.")', finish-start
  end subroutine


  !> Solves a linear system for a given timestep
  subroutine SolveLinearSystem(this,t,dt,y)
    implicit none
    type(bdf_type)   :: this
    double precision, dimension(this%nvector,this%neq) :: y
    double precision :: conv_error, error, dt, t, dt_scale, conv_rate
    logical          :: success,reset
    integer          :: conv_iterator, lte_iterator
    intent(inout)    :: y, dt, t, this

    ! 1. initialization ---------------------------------------------------------------------------!
    ! use the LU matrices from the predictor value
    ! todo: check if this%coeff(1,this%order)=beta or this%coeff(this%order,1)
    call GetLU(this,this%coeff(1,this%order),y,dt,this%L,this%U)

    ! Calculate initial right hand side
    call CalcRHS(this,this%y_NS(:,:,0),this%rhs)
    this%y_NS(:,:,1) = dt*this%rhs(:,:)

    ! some initializations
    this%den      = 0.0d0
    conv_rate     = 0.7d0
    conv_error    = 1d20
    lte_iterator  = 0
    conv_iterator = 0

    ! important quantity to measure density
    call CalcErrorWeight(this,this%rtol,this%atol,this%y,this%err_weight)

    ! 2. predictor --------------------------------------------------------------------------------!
    success = .FALSE.
    predictor: do while(.not. success)
      reset=.FALSE.

      ! predictor step - set predicted y_NS, en=0.0
      call PredictSolution(this,this%y_NS,this%en)
      y(:,:) = this%y_NS(:,:,0)
      ! 3. corrector ------------------------------------------------------------------------------!
      corrector: do while (conv_error > this%tautable(this%order,this%order)/(2d0*(this%order+2d0)))
        ! calculates residuum
        call CalcResiduum(this,y,dt,this%res)

        ! calculates the solution for dy with given residuum
        ! \todo this function needs to be enhanced for sparse matrices
        call SolveLU(this,this%L,this%U,this%res,this%den)

        ! add correction to solution vector
        this%en = this%en + this%den
        y = this%y_NS(:,:,0) + this%coeff(0,this%order)*this%en

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25 times the step-size
        call CheckConvergence(this,conv_iterator,conv_rate,this%den,conv_error,this%err_weight,reset)
        conv_iterator = conv_iterator + 1
        if (reset) then
          dt_scale = 0.25
          call ResetSystem(this,dt_scale,dt,this%y_NS)
          cycle predictor
        end if
      end do corrector

      ! local truncation error test:
      ! checks if solution is good enough and rewind everything if necessary
      error = WeightedNorm(this,this%en,this%err_weight)/this%tautable(this%order,this%order)
      if (error > 1d0) then
        dt_scale = 0.2
        call ResetSystem(this,dt_scale,dt,this%y_NS)
        lte_iterator = lte_iterator + 1
        cycle predictor
      else
        success = .TRUE.
      end if
    end do predictor

    ! rewrite history array
    call UpdateNordsieck(this, this%en, this%y_NS)

   ! write the result to the terminal or elsewhere
    t = t + dt
    print *, t, y(1,:)

    ! 4. sanity checks & step-size/order control --------------------------------------------------!
    ! calculate step size for current, upper & lower order and use largest one
    ! additionally it changes the order and enlarges the Nordsieck history array if necessary
    call CalcStepSizeOrder(this,this%order,this%err_weight,this%en,this%en_old,this%y_NS,dt_scale)

    ! Adjust the Nordsieck history array with new step size & new order
    call SetStepSize(this,dt_scale,this%y_NS,dt)
    this%en_old = this%en
  end subroutine SolveLinearSystem


  !>
  subroutine ResetSystem(this,dt_scale,dt,y_NS)
    implicit none
    type(bdf_type)   :: this
    integer          :: j,k,i,l
    double precision :: dt_scale, dt
    double precision, dimension(this%nvector,this%neq,0:6) :: y_NS
    intent(in)       :: dt_scale
    intent(inout)    :: y_NS, dt, this

    ! set the matrix back to old value
    do k = 0,this%order-1
      do j = this%order,k+1,-1
        do l = 1,this%neq
!cdir nodep
          do i = 1,this%nvector
            y_NS(i,l,j-1) = y_NS(i,l,j-1) - y_NS(i,l,j)
          end do
        end do
      end do
    end do

    call SetStepSize(this,dt_scale,y_NS,dt)
  end subroutine


  !> Calculates the step-sizes and uses according orders
  subroutine CalcStepSizeOrder(this,order,err_weight,en,en_old,y_NS,dt_scale)
    implicit none
    type(bdf_type) :: this
    integer :: order, dt_maxloc
    double precision, dimension(this%nvector,this%neq) :: en, err_weight, en_old
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: y_NS
    double precision  :: dt_scale, dt_scale_up, dt_scale_down, dt_upper_limit
    intent(in)        :: err_weight, en
    intent(inout)     :: order, this
    intent(out)       :: dt_scale

    ! calculate all estimated step sizes
    ! \todo check if order is correct here
    dt_scale_down = 1.d0/(1.3*(WeightedNorm(this,y_NS(:,:,order),err_weight)/this%tautable(order,order-1))**(1d0/order) + 1d-6)
    dt_scale = 1.d0/(1.2*(WeightedNorm(this,en,err_weight)/this%tautable(order,order))**(1d0/(order+1d0)) + 1d-6)
    dt_scale_up = 1.d0/(1.4*(WeightedNorm(this,en-en_old,err_weight)/this%tautable(order,order+1))**(1d0/(order+2d0)) + 1d-6)

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
    else if (dt_maxloc == 3) then ! enlarge Nordsieck array by 1
      if (order /= this%maxorder) then
        y_NS(:,:,order+1) = this%coeff(order,order)*en(:,:)/(order+1d0)
        order = order + 1
      else
        ! do nothing
      end if
    end if
  end subroutine


  !> Calculates the error for convergence
  subroutine CheckConvergence(this,iterator,conv_rate,den,error,weight,reset)
    implicit none
    type(bdf_type) :: this
    integer :: iterator
    logical :: reset
    double precision :: conv_rate, error_tmp, error
    double precision, dimension(this%nvector,this%neq) :: den, weight
    intent(in) :: den, weight, iterator, this
    intent(inout) :: conv_rate, error, reset

    error_tmp = WeightedNorm(this,den,weight)
    conv_rate = max(0.2d0*conv_rate,error_tmp/error)
    error = error_tmp*min(1d0,1.5d0*conv_rate)

    if (iterator == 2 .and. error > 2d0) then
      reset = .TRUE.
    else if (iterator > 2 .and. error > this%tautable(this%order,this%order)/(2d0*(this%order + 2d0))) then
      reset = .TRUE.
    else if (iterator > 3) then
      reset = .TRUE.
    else
      reset = .FALSE.
    end if
  end subroutine CheckConvergence


  !> Weighted Norm
  function WeightedNorm(this,en,weight)
    implicit none
    type(bdf_type) :: this
    double precision, dimension(this%nvector,this%neq) :: en, weight
    double precision :: WeightedNorm
    intent(in) :: en, weight

    WeightedNorm = sqrt(sum(en*en*weight*weight)/(this%neq*this%nvector))
  end function WeightedNorm


  !>
  subroutine CalcErrorWeight(this,rtol,atol,y,err_weight)
    implicit none
    type(bdf_type) :: this
    integer :: i, j
    double precision, dimension(this%nvector,this%neq) :: y,err_weight
    double precision  :: rtol, atol
    intent(in)        :: rtol, atol, y
    intent(out)       :: err_weight

    ! \todo this has to be different most probably for y
!cdir collapse
    do j=1,this%neq
      do i=1,this%nvector
        err_weight(i,j) = 1./(rtol*abs(y(i,j)) + atol)
      end do
    end do
  end subroutine CalcErrorWeight


  !> Sets the step size.
  !!
  !! When the step-size is changed also the corresponding history array needs to be adjusted. Since
  !! we are using the Nordsieck representation this is easily done by multiplying the track
  !! \f$r,rr,rrr,...\f$ to the Nordsieck array.
  subroutine SetStepSize(this,dt_scale,y_NS,dt)
    implicit none
    type(bdf_type)   :: this
    integer          :: i,j,k
    double precision :: dt_scale,dtt_scale,dt
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: y_NS
    intent(inout)    :: y_NS,dt,this
    intent(in)       :: dt_scale

    dtt_scale = 1.0
    do k = 1,this%order
      dtt_scale = dt_scale*dtt_scale
!cdir collapse
      do j = 1,this%neq
        do i = 1,this%nvector
          y_NS(i,j,k) = y_NS(i,j,k)*dtt_scale
        end do
      end do
    end do

    dt = dt*dt_scale
  end subroutine SetStepSize


  !>
  subroutine UpdateNordsieck(this,en,y_NS)
    implicit none
    type(bdf_type) :: this
    integer        :: i,j,k
    double precision, dimension(this%nvector,this%neq)                    :: en
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1)  :: y_NS
    intent(in)     :: en
    intent(inout)  :: y_NS, this

    do k = 0,this%order
!cdir collapse
      do j = 1,this%neq
        do i = 1,this%nvector
          y_NS(i,j,k) = y_NS(i,j,k) + en(i,j)*this%coeff(k,this%order)
        end do
      end do
    end do
  end subroutine UpdateNordsieck


  !> Predicts the solution
  subroutine PredictSolution(this,y_NS,en)
    implicit none
    type(bdf_type)  :: this
    integer         :: i,j,k,l
    double precision, dimension(this%nvector,this%neq,0:6) :: y_NS
    double precision, dimension(this%nvector,this%neq)     :: en
    intent (in)     :: this
    intent (inout)  :: y_NS
    intent (out)    :: en

    ! this loop effectively solves the Pascal Triangle without multiplications
    do l = 0,this%order-1
      do k = this%order,l+1,-1
!cdir collapse
        do j = 1,this%neq
          do i = 1,this%nvector
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
  subroutine SolveLU(this,L,U,res,den)
    implicit none
    type(bdf_type) :: this
    double precision, dimension(this%nvector,this%neq,this%neq) :: L, U
    double precision, dimension(this%nvector,this%neq)          :: res, den
    integer        :: i,j,k
    intent(in)     :: L, U, res
    intent(out)    :: den

    this%den_tmp(:,:) = 0.0
    ! lower system
    do k = 1,this%neq,1
      do j = 1,k-1
!cdir nodep
        do i = 1,this%nvector
          this%den_tmp(i,k) = this%den_tmp(i,k) + L(i,k,j)*this%den_tmp(i,j)
        end do
      end do
!cdir nodep
      this%den_tmp(:,k) = res(:,k) - this%den_tmp(:,k)
      this%den_tmp(:,k) = this%den_tmp(:,k)/L(:,k,k)
    end do

    ! upper system
    den(:,:) = 0.0
    do k = this%neq,1,-1
      do j = k+1,this%neq
!cdir nodep
        do i = 1,this%nvector
          den(i,k) = den(i,k) + U(i,k,j)*den(i,j)
        end do
      end do
      den(:,k) = this%den_tmp(:,k) - den(:,k)
      den(:,k) = den(:,k)/U(:,k,k)
    end do
  end subroutine SolveLU


  !> Get the Matrices L and U either from A or given
  !!
  !! \todo The Python preprocessor needs to include the code here
  subroutine GetLU(this,beta,y,dt,L,U)
    implicit none
    type(bdf_type) :: this
    integer :: i
    double precision :: beta, dt
    double precision, dimension(this%nvector,this%neq) :: y
    double precision, dimension(this%nvector,this%neq,this%neq) :: L,U
    intent (in) :: y, beta, dt
    intent (out) :: L, U

!cdir nodep
    do i=1,this%nvector
      L(i,1,1) = 1d0
      L(i,1,2) = 0d0
      L(i,1,3) = 0d0
      L(i,2,1) = 0.4d0*beta*dt/(0.04d0*beta*dt+1d0)
      L(i,2,2) = 1d0
      L(i,2,3) = 0d0
      L(i,3,1) = 0d0
      L(i,3,2) = - 60000000.0d0*beta*dt*y(i,2)/(-400.0d0*beta*beta*dt*dt*y(i,3)/ &
                  (0.04d0*beta*dt + 1d0) - beta*dt*(-60000000.0d0*y(i,2) - 10000.0d0*y(i,3)) + 1d0)
      L(i,3,3) = 1d0
      U(i,1,1) = 0.04d0*beta*dt + 1d0
      U(i,1,2) = -10000.0*beta*dt*y(i,3)
      U(i,1,3) = -10000.0*beta*dt*y(i,2)
      U(i,2,1) = 0d0
      U(i,2,2) = -400.0*beta**2*dt**2*y(i,3)/(0.04*beta*dt + 1) - beta*dt* &
                 (-60000000.0*y(i,2) - 10000.0*y(i,3)) + 1d0
      U(i,2,3) = -400.0*beta**2*dt**2*y(i,2)/(0.04*beta*dt + 1) + 10000.0*beta*dt*y(i,2)
      U(i,3,1) = 0d0
      U(i,3,2) = 0d0
      U(i,3,3) = 60000000.0*beta*dt*y(i,2)*(-400.0*beta**2*dt**2*y(i,2)/(0.04*beta*dt + 1d0) &
                + 10000.0*beta*dt*y(i,2))/(-400.0*beta**2*dt**2*y(i,3)/(0.04*beta*dt + 1d0) &
                - beta*dt*(-60000000.0*y(i,2) - 10000.0*y(i,3)) + 1d0) + 1d0
    end do
  end subroutine GetLU


  !> Calculates the residuum
  !!
  !! Calculates
  !! \f[
  !!    r = y_n - \sum_{i=1}^{order} \alpha_i * y_{n-i} - h \beta_0 f(t_n, y_n).
  !! \f]
  !! For a converging system \f$ r \approx 0 \f$ should become true.
  subroutine CalcResiduum(this,y,dt,res)
    implicit none
    type(bdf_type)    :: this
    integer           :: i,j
    double precision  :: dt
    double precision, dimension(this%nvector,this%neq) :: y, res
    intent (in)       :: y, dt
    intent (out)      :: res

    call CalcRHS(this,y,this%rhs)

!cdir collapse
    do j=1,this%neq
      do i=1,this%nvector
        res(i,j) = dt*this%rhs(i,j) - this%y_NS(i,j,1) - this%en(i,j)
      end do
    end do
  end subroutine CalcResiduum


  !> Example case for the right-hand-side
  !!
  !! \todo needs to be done external later!
  !! This is Robertson's example
  subroutine CalcRHS(this,y,rhs)
    implicit none
    type(bdf_type)  :: this
    integer         :: i
    double precision, dimension(this%nvector,this%neq) :: y
    double precision, dimension(this%nvector,this%neq) :: rhs
    intent (in)     :: y
    intent (out)    :: rhs

!cdir nodep
    do i=1,this%nvector
      rhs(i,1) = -0.04d0*y(i,1) + 1d4*y(i,2)*y(i,3)
      rhs(i,2) = 0.04d0*y(i,1) - 3d7*y(i,2)*y(i,2) - 1d4*y(i,2)*y(i,3)
      rhs(i,3) = 3d7*y(i,2)*y(i,2)
    end do
  end subroutine CalcRHS

  !>
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
