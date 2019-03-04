!> module that provides a cell-vectorized ODE-solver based on the BDF-scheme
module odevec_main
  implicit none

  type :: odevec
#ODEVEC_VECTORLENGTH
#ODEVEC_EQUATIONS
#ODEVEC_REACTIONS
#ODEVEC_MAXORDER
    integer :: iterator                                     !> iterator
    integer :: order                                        !> current order
!    integer :: successes, necessary_successes
    logical :: UpdateJac

    double precision :: rtol                                !> relative tolerance
    double precision :: atol                                !> absolute tolerance
    double precision :: t                                   !> current time
    double precision :: dt                                  !> current timestep
    double precision :: dt_min                              !> minimal timestep
    double precision, pointer, dimension(:,:) :: y          !> current solution
                                                            !> array | often y_NS(:,:,0)
    integer         , pointer, dimension(:)   :: Piv        !> pivoting vector
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
    double precision, pointer, dimension(:,:,:) :: LU       !> LU matrix, also jacobian
                                                            !! shortly (in-situ replace)
    double precision, pointer, dimension(:,:,:) :: y_NS     !> Nordsieck history array
#ODEVEC_LU_PRESENT

  end type

  ! two different interfaces, whether it is called with y or y_NS
  interface WeightedNorm
    module procedure WeightedNorm1, WeightedNorm2
  end interface

contains


  !> Allocates all memory for bdf-solver
  subroutine InitOdevec(this)
    implicit none
    type(odevec) :: this
    integer      :: err, i, j

    ! allocate fields
    allocate( &
              this%y(this%nvector,this%neq), &
              this%Piv(this%neq), &
              this%en(this%nvector,this%neq), &
              this%en_old(this%nvector,this%neq), &
              this%den(this%nvector,this%neq), &
              this%den_tmp(this%nvector,this%neq), &
              this%inv_weight_2(this%nvector,this%neq), &
              this%rhs(this%nvector,this%neq), &
              this%res(this%nvector,this%neq), &
              this%coeff(6,0:6), &
              this%tautable(this%maxorder+1,0:this%maxorder+1), &
              this%LU(this%nvector,this%neq,this%neq), &
              this%y_NS(this%nvector,this%neq,0:this%maxorder+1), &
              STAT=err)
    if (err.ne.0) then
      print *, "ODEVEC: Memory allocation error. OdeVec could not be intialized."
      stop
    end if

    ! bdf-coefficient-matrix with Nordsieck representation
    this%coeff = TRANSPOSE(reshape((/  &
    1d0       ,1d0, 0d0        ,0d0       ,0d0        ,0d0      ,0d0         , &
    2d0/3d0   ,1d0, 1d0/3d0    ,0d0       ,0d0        ,0d0      ,0d0         , &
    6d0/11d0  ,1d0, 6d0/11d0   ,1d0/11d0  ,0d0        ,0d0      ,0d0         , &
    12d0/25d0 ,1d0, 7d0/10d0   ,1d0/5d0   ,1d0/50d0   ,0d0      ,0d0         , &
    60d0/137d0,1d0, 225d0/274d0,85d0/274d0,15d0/274d0 ,1d0/274d0,0d0         , &
    20d0/49d0 ,1d0, 58d0/63d0  ,5d0/12d0  ,25d0/252d0 ,1d0/84d0 ,1d0/1764d0 /),&
    (/7, 6/)))

    ! precompute all tau and save them in a table (better performance)
    do i=1,this%maxorder+1
      do j=0,this%maxorder+1
        this%tautable(i,j) = tau(i,j,this%coeff)
      end do
    end do

#ODEVEC_DT_MIN

  end subroutine


  !> deallocates all memory of the bdf-solver
  subroutine CloseOdevec(this)
    implicit none
    type(odevec) :: this
    integer :: err

    deallocate(this%y,this%en,this%en_old,this%den,this%den_tmp,this%inv_weight_2, &
               this%rhs,this%res,this%coeff,this%tautable,this%LU, &
               this%y_NS, &
               stat=err)

    if (err.ne.0) then
      print *, "ODEVEC: Memory deallocation error. OdeVec was not properly closed."
      stop
    end if
  end subroutine


  !> Main subroutine to be called from outside
  subroutine SolveODE(this,time,dt,t_stop,y,GetRHS,GetJac,GetLU)
    implicit none
    type(odevec) :: this
    external          :: GetRHS,GetJac,GetLU
    double precision  :: time, t_stop, dt
    double precision, dimension(this%nvector,this%neq) :: y
    intent(in)        :: t_stop
    intent(inout)     :: y, time, dt

    ! some initializations
    this%order       = 1
    this%y_NS(:,:,0) = y(:,:)
    this%iterator    = 0

    ! main solve - solve the linear system
    do while (time < t_stop)
      ! solve the system
      call SolveLinearSystem(this,time,dt,y,GetRHS,GetJac,GetLU)

      this%iterator = this%iterator + 1
    end do
    call InterpolateSolution(this,time,dt,t_stop,this%order,this%y_NS,y)
  end subroutine

  !> Solves a linear system for a given timestep
  subroutine SolveLinearSystem(this,t,dt,y,GetRHS,GetJac,GetLU)
    implicit none
    type(odevec) :: this
    external          :: GetRHS,GetJac,GetLU
    double precision, dimension(this%nvector,this%neq) :: y
    double precision :: dt, t, dt_scale
    integer          :: i,j,k
    double precision :: conv_error, conv_rate, error
    integer          :: conv_iterator, conv_failed, lte_iterator
    logical          :: success, reset
    intent(inout)    :: y, dt, t

    ! 1. initialization -------------------------------------------------------!
    ! use the LU matrices from the predictor value

    ! Calculate initial right hand side
    call GetRHS(this,y,this%rhs)
    this%y_NS(:,:,1) = dt*this%rhs(:,:)

    ! some initializations
    this%den      = 0.0d0
    conv_rate     = 0.7d0
    conv_error    = HUGE(1d0)
    lte_iterator  = 0
    conv_iterator = 0
    conv_failed   = 0
    this%UpdateJac = .TRUE.
!    this%successes     = 0
!    this%necessary_successes = this%order + 1

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(this,this%rtol,this%atol,y,this%inv_weight_2)

    ! 2. predictor ------------------------------------------------------------!
    success = .FALSE.
    predictor: do while(.not. success)
      reset = .FALSE.

      ! predictor step - set predicted y_NS, en=0.0
      call PredictSolution(this,this%y_NS,this%en)
      y(:,:) = this%y_NS(:,:,0)

      if (this%UpdateJac) then
        if (this%LU_PRESENT) then
          call GetLU(this,this%coeff(this%order,0),y,dt,this%LU,this%Piv)
          this%Piv = [(i, i=1, this%neq)]
        else
          call GetJac(this,this%coeff(this%order,0),y,dt,this%LU)
          call LUDecompose(this,this%LU,this%Piv)
        end if
        this%UpdateJac = .FALSE.
      end if

      conv_iterator = 0
      ! 3. corrector ----------------------------------------------------------!
      corrector: do while ((conv_error > this%tautable(this%order,this%order)/ &
                                        (2d0*(this%order+2d0))) .or. &
                           (conv_iterator .eq. 0))

        ! calculates residuum
        call CalcResiduum(this,GetRHS,y,dt,this%res)

        ! calculates the solution for dy with given residuum
        call SolveLU(this,this%LU,this%Piv,this%res,this%den)

        ! add correction to solution vector
        this%en(:,:) = this%en(:,:) + this%den(:,:)
        y(:,:) = this%y_NS(:,:,0) + this%coeff(this%order,0)*this%en

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25
        ! times the step-size
        call CheckConvergence(this,conv_iterator,conv_rate,this%den, &
                              conv_error,this%inv_weight_2,reset)

        conv_iterator = conv_iterator + 1
        if (reset) then
          dt_scale = 0.25
          conv_failed = conv_failed + 1
          if (conv_failed .gt. 10) then
            print *, "ODEVEC: Convergence failed! Abortion after more than 10 iterations."
            stop
          end if
          conv_iterator = 0
          call ResetSystem(this,dt_scale,dt,this%y_NS)
          if (dt .lt. this%dt_min) then
            print *, "ODEVEC: Convergence failed! Abortion, because timestep too small."
            stop
          end if
          this%UpdateJac = .TRUE.
          cycle predictor
        end if
      end do corrector
      conv_failed = 0

      ! local truncation error test:
      ! checks if solution is good enough and, similar to the convergence
      ! test, rerun from predictor with adjusted step size
      error = WeightedNorm(this,this%en,this%inv_weight_2)/ &
              this%tautable(this%order,this%order)
      if (error > 1d0) then
        dt_scale = 0.2
        call ResetSystem(this,dt_scale,dt,this%y_NS)
        !if (dt .lt. this%dt_min) stop
        lte_iterator = lte_iterator + 1
        this%UpdateJac = .TRUE.
        cycle predictor
      else
        success = .TRUE.
!        this%successes = this%successes + 1
      end if
    end do predictor

    ! after a successfull run rewrite the history array
    call UpdateNordsieck(this,this%en,this%y_NS)

    ! advance in time
    t = t + dt

!    if (this%successes .ge. this%necessary_successes) then
    ! 4. step-size/order control ----------------------------------------------!
    ! calc. step size for order+(0,-1,+1) -> use largest for next step
    call CalcStepSizeOrder(this,dt_scale,this%order)

    ! Adjust the Nordsieck history array with new step size & new order
    call SetStepSize(this,dt_scale,this%y_NS,dt)

!    this%necessary_successes = this%order + 1
!    this%successes = 0
!    end if
    this%en_old = this%en
  end subroutine SolveLinearSystem

  !> Resets the whole system to the state where the predictor starts
  subroutine ResetSystem(this,dt_scale,dt,y_NS)
    implicit none
    type(odevec)     :: this
    integer          :: j,k,i,l
    double precision :: dt_scale, dt
    double precision, dimension(this%nvector,this%neq,0:6) :: y_NS
    intent(in)       :: dt_scale
    intent(inout)    :: y_NS, dt

    ! set the matrix back to old value
    do k = 0,this%order-1
      do j = this%order,k+1,-1
!NEC$ collapse
        do l = 1,this%neq
          do i = 1,this%nvector
            y_NS(i,l,j-1) = y_NS(i,l,j-1) - y_NS(i,l,j)
          end do
        end do
      end do
    end do

    call SetStepSize(this,dt_scale,y_NS,dt)
  end subroutine


  !> Calculates the step-sizes and uses according orders
  subroutine CalcStepSizeOrder(this,dt_scale,order)
    implicit none
    type(odevec)     :: this
    integer          :: dt_maxloc, order
    double precision :: dt_scale, dt_scale_up, dt_scale_down, dt_upper_limit
    intent(out)      :: dt_scale
    intent(inout)    :: order,this

    ! calculate all estimated step sizes
    ! todo err_weight_2 = err_weight*err_weight
    dt_scale_down = 1.d0/(1.3*(WeightedNorm(this,this%y_NS,this%inv_weight_2,order)/ &
                          this%tautable(order,order-1))**(1d0/order) + 1d-6)
    dt_scale      = 1.d0/(1.2*(WeightedNorm(this,this%en,this%inv_weight_2)/ &
                          this%tautable(order,order))**(1d0/(order+1d0)) + 1d-6)
    dt_scale_up   = 1.d0/(1.4*(WeightedNorm(this,this%en-this%en_old,this%inv_weight_2)/ &
                          this%tautable(order,order+1))**(1d0/(order+2d0)) + 1d-6)

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
      if (order /= this%maxorder) then
        this%y_NS(:,:,order+1) = this%coeff(order,order)*this%en(:,:)/(order+1d0)
        order = order + 1
      else
        ! do nothing
      end if
    end if
  end subroutine


  !> Calculates the error for convergence
  subroutine CheckConvergence(this,iterator,conv_rate,den,conv_error,inv_weight_2,reset)
    implicit none
    type(odevec)     :: this
    integer          :: iterator
    logical          :: reset
    double precision :: conv_rate, conv_rate2, conv_error_tmp, conv_error
    double precision, dimension(this%nvector,this%neq) :: den, inv_weight_2
    intent(in)       :: den, inv_weight_2, iterator
    intent(inout)    :: conv_rate, conv_error, reset

    conv_error_tmp = WeightedNorm(this,den,inv_weight_2)
    conv_rate2 = conv_error_tmp/conv_error
    conv_rate = max(0.2d0*conv_rate,conv_rate2)
    conv_error = conv_error_tmp*min(1d0,1.5d0*conv_rate)

    if (iterator == 2 .and. conv_rate2 > 2d0) then
      reset = .TRUE.
    else if (iterator > 2 .and.  &
            conv_error > this%tautable(this%order,this%order)/(2d0*(this%order + 2d0))) then
      reset = .TRUE.
    else if (iterator > 3) then
      reset = .TRUE.
    else
      reset = .FALSE.
    end if
  end subroutine CheckConvergence


  !> Calculates the weighted norm (two different interfaces for performance)
  pure function WeightedNorm1(this,en,inv_weight_2)
    implicit none
    type(odevec)     :: this
    double precision :: WeightedNorm1
    double precision, dimension(this%nvector,this%neq) :: en, inv_weight_2
    intent(in)       :: en,inv_weight_2,this

    WeightedNorm1 = sqrt(sum(en(:,:)*en(:,:)*inv_weight_2(:,:))/ &
                         (this%neq*this%nvector))
  end function WeightedNorm1
  pure function WeightedNorm2(this,en,inv_weight_2,column)
    implicit none
    type(odevec)     :: this
    integer          :: column
    double precision :: WeightedNorm2
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: en
    double precision, dimension(this%nvector,this%neq) :: inv_weight_2
    intent(in)       :: en,inv_weight_2,column,this

    WeightedNorm2 = sqrt(sum(en(:,:,column)*en(:,:,column)*inv_weight_2(:,:))/ &
                         (this%neq*this%nvector))
  end function WeightedNorm2


  !> Calculates the error-weight of y for given tolerances
  subroutine CalcErrorWeightInvSquared(this,rtol,atol,y,inv_weight_2)
    implicit none
    type(odevec)     :: this
    double precision :: rtol, atol
    double precision, dimension(this%nvector,this%neq) :: y, inv_weight_2
    intent(in)       :: rtol, atol, y
    intent(out)      :: inv_weight_2

    inv_weight_2(:,:) = 1./(rtol*abs(y(:,:)) + atol)
    inv_weight_2(:,:) = inv_weight_2(:,:)*inv_weight_2(:,:)
  end subroutine CalcErrorWeightInvSquared


  !> Sets the step size and scales the Nordsieck-array accordingly
  subroutine SetStepSize(this,dt_scale,y_NS,dt)
    implicit none
    type(odevec)     :: this
    integer          :: i,j,k
    double precision :: dt_scale,dtt_scale,dt
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: y_NS
    intent(in)       :: dt_scale,this
    intent(inout)    :: y_NS,dt

    dtt_scale = 1.0
    do k = 1,this%order
      dtt_scale = dt_scale*dtt_scale
!NEC$ collapse
      do j = 1,this%neq
        do i = 1,this%nvector
          y_NS(i,j,k) = y_NS(i,j,k)*dtt_scale
        end do
      end do
    end do

    dt = dt*dt_scale
  end subroutine SetStepSize


  !> Sets the new Nordsieck array for a given correction vector
  subroutine UpdateNordsieck(this,en,y_NS)
    implicit none
    type(odevec)   :: this
    integer        :: i,j,k
    double precision, dimension(this%nvector,this%neq)               :: en
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1)  :: y_NS
    intent(in)     :: en
    intent(inout)  :: y_NS

    do k = 0,this%order
!NEC$ collapse
      do j = 1,this%neq
        do i = 1,this%nvector
          y_NS(i,j,k) = y_NS(i,j,k) + en(i,j)*this%coeff(this%order,k)
        end do
      end do
    end do
  end subroutine UpdateNordsieck


  !> Predicts the solution
  subroutine PredictSolution(this,y_NS,en)
    implicit none
    type(odevec)    :: this
    double precision, dimension(this%nvector,this%neq,0:6) :: y_NS
    double precision, dimension(this%nvector,this%neq)     :: en
    integer         :: i,j,k,l
    intent (inout)  :: y_NS
    intent (out)    :: en

    !  loop effectively solves the Pascal Triangle without multiplications
    do l = 0,this%order-1
      do k = this%order,l+1,-1
!NEC$ collapse
        do j = 1,this%neq
          do i = 1,this%nvector
            y_NS(i,j,k-1) = y_NS(i,j,k) + y_NS(i,j,k-1)
          end do
        end do
      end do
    end do

    ! set the correction vector to zero
    en(:,:) = 0.0
  end subroutine PredictSolution


  !> Solve the System \f$ LU x=r \f$ with residuum r (res) and return x (den)
  subroutine SolveLU(this,LU,Piv,res,den)
    implicit none
    type(odevec)   :: this
    double precision, dimension(this%nvector,this%neq,this%neq) :: LU
    double precision, dimension(this%nvector,this%neq)          :: res, den
    integer         , dimension(this%neq)                       :: Piv
    integer        :: i,j,k
    intent(in)     :: LU,Piv,res
    intent(out)    :: den

!NEC$ collapse
    do j=1,this%neq
      do i=1,this%nvector
        this%den_tmp(i,j) = 0.0
      end do
    end do
    ! lower system
    do k = 1,this%neq,1
      do j = 1,k-1
!NEC$ ivdep
        do i = 1,this%nvector
          this%den_tmp(i,k) = this%den_tmp(i,k) + LU(i,Piv(k),j)*this%den_tmp(i,j)
        end do
      end do
!NEC$ ivdep
      do i = 1,this%nvector
        this%den_tmp(i,k) = res(i,k) - this%den_tmp(i,k)
        this%den_tmp(i,k) = this%den_tmp(i,k)!/LU(i,Piv(k),k)
      end do
    end do

    ! upper system
!NEC$ collapse
    do j=1,this%neq
      do i=1,this%nvector
        den(i,j) = 0.0
      end do
    end do
    do k = this%neq,1,-1
      do j = k+1,this%neq
!NEC$ ivdep
        do i = 1,this%nvector
          den(i,k) = den(i,k) + LU(i,Piv(k),j)*den(i,j)
        end do
      end do
!NEC$ ivdep
      do i = 1,this%nvector
        den(i,k) = this%den_tmp(i,k) - den(i,k)
        den(i,k) = den(i,k)/LU(i,Piv(k),k)
      end do
    end do
  end subroutine SolveLU


  !> Get L and U from Jacobian (in-place transformation)
  !!
  !! Based on dgetrf from LAPACK. LU-decomposition with partial pivoting. Taken
  !! and modified from https://rosettacode.org/wiki/LU_decomposition#Fortran.
  subroutine LUDecompose(this,A,P)
    implicit none
    type(odevec)   :: this
    integer        :: i, j, k, m, kmax
    double precision, dimension(this%nvector,this%neq,this%neq) :: A
    integer,          dimension(this%neq) :: P
    integer,          dimension(2) :: maxloc_ij
    intent (inout) :: A
    intent (out)   :: p

    ! initialize P
    P = [(i, i=1, this%neq)]
    do k = 1,this%neq-1
      maxloc_ij(:) = maxloc(abs(A(:,P(k:),k)))
      kmax = maxloc_ij(2) + k - 1
      if (kmax /= k ) P([k, kmax]) = P([kmax, k])
    end do

    do k = 1,this%neq-1
      do j = k+1,this%neq
        do i = 1,this%nvector
          A(i,P(j),k) = A(i,P(j),k) / A(i,P(k),k)
        end do
      end do
      do j = k+1,this%neq
        do m = k+1,this%neq
!NEC$ ivdep
          do i = 1,this%nvector
            A(i,P(m),j) = A(i,P(m),j) - A(i,P(m),k) * A(i,P(k),j)
          end do
        end do
      end do
    end do

  end subroutine LUDecompose


  !> Calculates the residuum
  subroutine CalcResiduum(this,GetRHS,y,dt,res)
    implicit none
    type(odevec)      :: this
    external          :: GetRHS
    integer           :: i,j
    double precision  :: dt
    double precision, dimension(this%nvector,this%neq) :: y, res
    intent (in)       :: y, dt
    intent (out)      :: res

    call GetRHS(this,y,this%rhs)

!NEC$ collapse
    do j=1,this%neq
      do i=1,this%nvector
        res(i,j) = dt*this%rhs(i,j) - this%y_NS(i,j,1) - this%en(i,j)
      end do
    end do
  end subroutine CalcResiduum


  !> Interpolates solution at output time
  !!
  !! If an output is requested because \$ t > t_{\mathrm{out}} \$ the solution
  !! is interpolated at the chosen output position. The interpolation is done
  !! by Taylor series expansion, where Horner's rule is applied for efficiency.
  subroutine InterpolateSolution(this,t,dt,t_out,order,y_NS,y)
    implicit none
    type(odevec)    :: this
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: y_NS
    double precision, dimension(this%nvector,this%neq) :: y
    double precision  :: t, t_out, rel_t, dt
    integer           :: k, k_dot, order
    intent(in)        :: dt, t_out, order, y_NS
    intent(inout)     :: t, y

    ! Horner's rule
    y(:,:) = y_NS(:,:,order)
    rel_t = (t_out - t)/(dt)
    do k = 1,order
      k_dot = order - k
      y(:,:) = y_NS(:,:,k_dot) + rel_t*y(:,:)
    end do
    t = t_out
  end subroutine


  !> taufunction needed in several tests
  function tau(order,order2,coeff)
    implicit none
    integer :: i, order, order2
    double precision :: tau, faculty
    double precision, dimension(0:6,6) :: coeff
    intent(in) :: coeff, order

    faculty=1d0
!NEC$ shortloop
    do i=1,order
      faculty = faculty*i
    end do

    tau = (order2+1d0)/(faculty*coeff(order,order))

    return
  end function tau

end module odevec_main
