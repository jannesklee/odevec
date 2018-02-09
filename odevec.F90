!> module that provides a cell-vectorized ODE-solver based on the BDF-scheme
module bdf_method_mod
  implicit none

  type :: odevec
    integer :: nvector=6                         !> vector length
    integer :: neq=12                            !> number equations
    integer :: maxorder=5                        !> maximum order
    integer :: nrea = 38 ! delete again
    integer :: iterator                          !> iterator
    integer :: order                             !> current order

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
    double precision, pointer, dimension(:,:,:) :: jac      !> jacobian
    double precision, pointer, dimension(:,:,:) :: L,U      !> lower, upper triangular matrix
    double precision, pointer, dimension(:,:,:) :: y_NS     !> Nordsieck history array

#ODEVEC_LU_PRESENT

  end type

  ! two different interfaces, whether it is called with y or y_NS
  interface WeightedNorm
    module procedure WeightedNorm1, WeightedNorm2
  end interface

contains


  !> Allocates all memory for bdf-solver
  subroutine InitBDF(this)
    implicit none
    type(odevec) :: this
    integer      :: err, i, j

    ! allocate fields
    allocate( &
              this%y(this%nvector,this%neq), &
              this%en(this%nvector,this%neq), &
              this%en_old(this%nvector,this%neq), &
              this%den(this%nvector,this%neq), &
              this%den_tmp(this%nvector,this%neq), &
              this%inv_weight_2(this%nvector,this%neq), &
              this%rhs(this%nvector,this%neq), &
              this%res(this%nvector,this%neq), &
              this%coeff(0:6,6), &
              this%tautable(this%maxorder+1,0:this%maxorder+1), &
              this%jac(this%nvector,this%neq,this%neq), &
              this%L(this%nvector,this%neq,this%neq), &
              this%U(this%nvector,this%neq,this%neq), &
              this%y_NS(this%nvector,this%neq,0:this%maxorder+1), &
              STAT=err)
    if (err.ne.0) then
      print *, "Memory allocation error. BDF-Solver could not be intialized."
      stop
    end if

    ! bdf-coefficient-matrix with Nordsieck representation
    this%coeff = reshape((/  &
    1d0       ,1d0,0d0        ,0d0       ,0d0        ,0d0      ,0d0         , &
    2d0/3d0   ,1d0,1d0/3d0    ,0d0       ,0d0        ,0d0      ,0d0         , &
    6d0/11d0  ,1d0,6d0/11d0   ,1d0/11d0  ,0d0        ,0d0      ,0d0         , &
    12d0/25d0 ,1d0,7d0/10d0   ,1d0/5d0   ,1d0/50d0   ,0d0      ,0d0         , &
    60d0/137d0,1d0,225d0/274d0,85d0/274d0,15d0/274d0 ,1d0/274d0,0d0         , &
    20d0/49d0 ,1d0,58d0/63d0  ,5d0/12d0  ,25d0/252d0 ,1d0/84d0 ,1d0/1764d0 /),&
    (/7, 6/))

    ! precompute all tau and save them in a table (better performance)
    do i=1,this%maxorder+1
      do j=0,this%maxorder+1
        this%tautable(i,j) = tau(i,j,this%coeff)
      end do
    end do
  end subroutine


  !> deallocates all memory of the bdf-solver
  subroutine CloseBDF(this)
    implicit none
    type(odevec) :: this
    integer :: err

    deallocate(this%y,this%en,this%en_old,this%den,this%den_tmp,this%inv_weight_2, &
               this%rhs,this%res,this%coeff,this%tautable,this%jac, &
               this%L,this%U,this%y_NS, &
               stat=err)

    if (err.ne.0) then
      print *, "Memory deallocation error. BDF-Solver was not properly closed."
      stop
    end if
  end subroutine


  !> Main subroutine to be called from outside
  subroutine SolveODE_BDF(this,time,dt,t_stop,y)
    implicit none
    type(odevec) :: this
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
      call SolveLinearSystem(this,time,dt,y)

      this%iterator = this%iterator + 1
    end do
    call InterpolateSolution(this,time,dt,t_stop,this%order,this%y_NS,y)
  end subroutine


  !> Solves a linear system for a given timestep
  subroutine SolveLinearSystem(this,t,dt,y)
    implicit none
    type(odevec) :: this
    double precision, dimension(this%nvector,this%neq) :: y
    double precision :: dt, t, dt_scale
    double precision :: conv_error, conv_rate, error
    integer          :: conv_iterator, lte_iterator
    logical          :: success, reset
    intent(inout)    :: y, dt, t

    ! 1. initialization -------------------------------------------------------!
    ! use the LU matrices from the predictor value

    if (this%LU_PRESENT) then
      call GetLU(this,this%coeff(1,this%order),y,dt,this%L,this%U)
    else
      call GetJac(this,this%coeff(1,this%order),y,dt,this%jac)
      call LUDecompose(this,this%jac,this%L,this%U)
    end if

    ! Calculate initial right hand side
    call CalcRHS(this,y,this%rhs)
    this%y_NS(:,:,1) = dt*this%rhs(:,:)


    ! some initializations
    this%den           = 0.0d0
    conv_rate     = 0.7d0
    conv_error    = 1d20
    lte_iterator  = 0
    conv_iterator = 0

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(this,this%rtol,this%atol,y,this%inv_weight_2)

    ! 2. predictor ------------------------------------------------------------!
    success = .FALSE.
    predictor: do while(.not. success)
      reset = .FALSE.

      ! predictor step - set predicted y_NS, en=0.0
      call PredictSolution(this,this%y_NS,this%en)
      y(:,:) = this%y_NS(:,:,0)

      ! 3. corrector ----------------------------------------------------------!
      corrector: do while (conv_error > this%tautable(this%order,this%order)/ &
                                        (2d0*(this%order+2d0)))
        ! calculates residuum
        call CalcResiduum(this,y,dt,this%res)

        ! calculates the solution for dy with given residuum
        call SolveLU(this,this%L,this%U,this%res,this%den)

        ! add correction to solution vector
        this%en(:,:) = this%en(:,:) + this%den(:,:)
        y(:,:) = this%y_NS(:,:,0) + this%coeff(0,this%order)*this%en

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25
        ! times the step-size
        call CheckConvergence(this,conv_iterator,conv_rate,this%den, &
                              conv_error,this%inv_weight_2,reset)
        conv_iterator = conv_iterator + 1
        if (reset) then
          dt_scale = 0.25
          call ResetSystem(this,dt_scale,dt,this%y_NS)
          cycle predictor
        end if
      end do corrector

      ! local truncation error test:
      ! checks if solution is good enough and, similar to the convergence
      ! test, rerun from predictor with adjusted step size
      error = WeightedNorm(this,this%en,this%inv_weight_2)/ &
              this%tautable(this%order,this%order)
      if (error > 1d0) then
        dt_scale = 0.2
        call ResetSystem(this,dt_scale,dt,this%y_NS)
        lte_iterator = lte_iterator + 1
        cycle predictor
      else
        success = .TRUE.
      end if
    end do predictor

    ! after a successfull run rewrite the history array
    call UpdateNordsieck(this,this%en,this%y_NS)

    ! advance in time
    t = t + dt

    ! 4. step-size/order control ----------------------------------------------!
    ! calc. step size for current order+(0,-1,+1) -> use largest for next step
    call CalcStepSizeOrder(this,dt_scale,this%order)

    ! Adjust the Nordsieck history array with new step size & new order
    call SetStepSize(this,dt_scale,this%y_NS,dt)
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
!cdir collapse
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
  subroutine CheckConvergence(this,iterator,conv_rate,den,error,inv_weight_2,reset)
    implicit none
    type(odevec)     :: this
    integer          :: iterator
    logical          :: reset
    double precision :: conv_rate, error_tmp, error
    double precision, dimension(this%nvector,this%neq) :: den, inv_weight_2
    intent(in)       :: den, inv_weight_2, iterator
    intent(inout)    :: conv_rate, error, reset

    error_tmp = WeightedNorm(this,den,inv_weight_2)
    conv_rate = max(0.2d0*conv_rate,error_tmp/error)
    error = error_tmp*min(1d0,1.5d0*conv_rate)

    if (iterator == 2 .and. error > 2d0) then
      reset = .TRUE.
    else if (iterator > 2 .and.  &
            error > this%tautable(this%order,this%order)/(2d0*(this%order + 2d0))) then
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
!cdir collapse
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
    type(odevec)    :: this
    double precision, dimension(this%nvector,this%neq,0:6) :: y_NS
    double precision, dimension(this%nvector,this%neq)     :: en
    integer         :: i,j,k,l
    intent (inout)  :: y_NS
    intent (out)    :: en

    !  loop effectively solves the Pascal Triangle without multiplications
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
    en(:,:) = 0.0
  end subroutine PredictSolution


  !> Solve the System \f$ LU x=r \f$ with residuum r (res) and return x (den)
  subroutine SolveLU(this,L,U,res,den)
    implicit none
    type(odevec)   :: this
    double precision, dimension(this%nvector,this%neq,this%neq) :: L, U
    double precision, dimension(this%nvector,this%neq)          :: res, den
    integer        :: i,j,k
    intent(in)     :: L, U, res
    intent(out)    :: den

!cdir collapse
    do j=1,this%neq
      do i=1,this%nvector
        this%den_tmp(i,j) = 0.0
      end do
    end do
    ! lower system
    do k = 1,this%neq,1
      do j = 1,k-1
!cdir nodep
        do i = 1,this%nvector
          this%den_tmp(i,k) = this%den_tmp(i,k) + L(i,k,j)*this%den_tmp(i,j)
        end do
      end do
!cdir nodep
      do i = 1,this%nvector
        this%den_tmp(i,k) = res(i,k) - this%den_tmp(i,k)
        this%den_tmp(i,k) = this%den_tmp(i,k)/L(i,k,k)
      end do
    end do

    ! upper system
!cdir collapse
    do j=1,this%neq
      do i=1,this%nvector
        den(i,j) = 0.0
      end do
    end do
    do k = this%neq,1,-1
      do j = k+1,this%neq
!cdir nodep
        do i = 1,this%nvector
          den(i,k) = den(i,k) + U(i,k,j)*den(i,j)
        end do
      end do
!cdir nodep
      do i = 1,this%nvector
        den(i,k) = this%den_tmp(i,k) - den(i,k)
        den(i,k) = den(i,k)/U(i,k,k)
      end do
    end do
  end subroutine SolveLU


  !> Provides NOT the jacobian J but P = 1 - dt*beta*J
  subroutine GetJac(this,beta,y,dt,jac)
    implicit none
    type(odevec)     :: this
    double precision :: beta, dt
    double precision, dimension(this%nvector,this%neq) :: y
    double precision, dimension(this%nvector,this%neq,this%neq) :: jac
    double precision, dimension(this%nvector,this%nrea) :: k
    integer          :: i
    intent (in)      :: y, beta, dt
    intent (out)     :: jac

    ! \todo remove again
    k = 0.0d0
    k = coe(this)

!cdir nodep
    do i=1,this%nvector

#ODEVEC_JAC

    end do

  end subroutine GetJac


  !> Get L and U from Jacobian
  subroutine LUDecompose(this,A,L,U)
    implicit none
    type(odevec)   :: this
    integer        :: i, j, k, m
    double precision, dimension(this%nvector,this%neq,this%neq) :: A
    double precision, dimension(this%nvector,this%neq,this%neq) :: L, U
    intent (inout) :: A
    intent (out)   :: L, U

    ! calculate L
    do k=1,this%neq-1
!cdir nodep
      do j=k+1,this%neq
        do i=1,this%nvector
          L(i,j,k) = A(i,j,k)/A(i,k,k)
        end do
!cdir nodep
        do m=k+1,this%neq
          do i=1,this%nvector
            A(i,j,m) = A(i,j,m)-L(i,j,k)*A(i,k,m)
          end do
        end do
      end do
    end do

    do j=1,this%neq
      do i=1,this%nvector
        L(i,j,j) = 1.0
      end do
    end do

    ! calculate U
    do k=1,this%neq
!cdir nodep
      do j=1,k
        do i=1,this%nvector
          U(i,j,k) = A(i,j,k)
        end do
      end do
    end do
  end subroutine LUDecompose


  !> Provides the Matrices L and U
  subroutine GetLU(this,beta,y,dt,L,U)
    implicit none
    type(odevec), intent(in) :: this
    double precision, dimension(this%nvector,this%neq), &
         intent(in)  :: y
    double precision, dimension(this%nvector,this%neq,this%neq), &
         intent(out) :: L,U
    double precision :: beta, dt
    integer          :: i

!cdir nodep
    do i=1,this%nvector

#ODEVEC_L

#ODEVEC_U

    end do
  end subroutine GetLU


  !> Calculates the residuum
  subroutine CalcResiduum(this,y,dt,res)
    implicit none
    type(odevec)      :: this
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
  !! This is for Robertson's example
  subroutine CalcRHS(this,y,rhs)
    implicit none
    type(odevec)    :: this
    integer         :: i
    double precision, dimension(this%nvector,this%neq) :: y
    double precision, dimension(this%nvector,this%neq) :: rhs
    double precision, dimension(this%nvector,this%nrea) :: k
    intent (in)     :: y
    intent (out)    :: rhs

    k = 0.0d0
    k = coe(this)

!cdir nodep
    do i=1,this%nvector
#ODEVEC_RHS
    end do
  end subroutine CalcRHS

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


  !compute reaction rates cm^3(n-1)/s
  function coe(this)
    implicit none
    type(odevec)    :: this
    real*8::coe(this%nvector,this%nrea)
    real*8::Tgas(this%nvector)
    real*8::Te(this%nvector)
    real*8::lnTe(this%nvector)
    real*8::T32(this%nvector)
    real*8::invT(this%nvector)
    real*8::invTe(this%nvector)
    real*8::sqrTgas(this%nvector)
    real*8::log10Te(this%nvector)
    real*8::small(this%nvector),nmax(this%nvector)
    real*8::vec_k(this%nvector,this%nrea)
    real*8::Tmax,Tmin
    integer::i

    !Tgas is in K

!    Tmin=100d0

   ! replaced with a proper value according to the environment
    do i=1,this%nvector
      nmax(i) = max(maxval(this%y(i,1:this%neq)),1d0)
    end do
    small = 1d-40/(nmax*nmax*nmax)

    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
    lnTe = log(Te) !ln of Te (#)
    T32 = Tgas*0.0033333333333333335 !Tgas/(300 K) (#)
    invT = 1.d0/Tgas !inverse of T (1/K)
    invTe = 1.d0/Te !inverse of T (1/eV)
    sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)

    do i=1,this%nvector
      vec_k(i,:) = small(i) !inizialize coefficients
    end do

    !H + E -> H+ + E + E
    vec_k(:,1) = small + (exp(-32.71396786d0+13.5365560d0&
        *lnTe-5.73932875d0*(lnTe**2)+1.56315498d0&
        *(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2&
        *(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4&
        *(lnTe**7)-2.03914985d-6*(lnTe**8)))

    !H+ + E -> H
    where(Tgas.LE.5.5d3)
      vec_k(:,2) = small + (3.92d-13&
          *invTe**0.6353d0)
    end where

    !H+ + E -> H
    where(Tgas.GT.5.5d3)
      vec_k(:,3) = small + (exp(-28.61303380689232d0-0.7241125657826851d0&
          *lnTe-0.02026044731984691d0*lnTe**2-0.002380861877349834d0&
          *lnTe**3-0.0003212605213188796d0&
          *lnTe**4-0.00001421502914054107d0&
          *lnTe**5+4.989108920299513d-6*lnTe**6+5.755614137575758d-7&
          *lnTe**7-1.856767039775261d-8*lnTe**8-3.071135243196595d-9&
          *lnTe**9))
    end where

    !HE + E -> HE+ + E + E
    vec_k(:,4) = small + (dexp(-44.09864886d0+23.91596563d0&
        *lnTe-10.7532302d0*(lnTe**2)+3.05803875d0&
        *(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2&
        *(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4&
        *(lnTe**7)-3.64916141d-6*(lnTe**8)))

    !HE+ + E -> HE
    where(Tgas.LE.9.28d3)
      vec_k(:,5) = small + (3.92d-13&
          *invTe**0.6353d0)
    end where

    !HE+ + E -> HE
    where(Tgas.GT.9.28d3)
      vec_k(:,6) = small + (1.54d-9&
          *(1.d0+0.3d0&
          /exp(8.099328789667d0&
          *invTe))&
          /(exp(40.49664394833662d0*invTe)&
          *Te**1.5d0)+3.92d-13&
          /Te**0.6353d0)
    end where

    !HE+ + E -> HE++ + E + E
    vec_k(:,7) = small + (exp(-68.71040990212001d0+43.93347632635d0&
        *lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0&
        *lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0&
        *lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0&
        *lnTe**7-3.165581065665d-6*lnTe**8))

    !HE++ + E -> HE+
    vec_k(:,8) = small + (3.36d-10/sqrTgas/(Tgas&
        /1.d3)**0.2d0/(1+(Tgas/1.d6)**0.7d0))

    !H + E -> H-
    vec_k(:,9) = small + (6.77d-15*Te**0.8779d0)

    !H- + H -> H2 + E
    where(Tgas.LT.1160d0)
      vec_k(:,10) = small + (1.43d-9)
    end where

    !H- + H -> H2 + E
    where(Tgas.GT.1160d0)
      vec_k(:,11) = small + (exp(-20.06913897587003d0+0.2289800603272916d0&
          *lnTe+0.03599837721023835d0*lnTe**2-0.004555120027032095d0&
          *lnTe**3-0.0003105115447124016d0&
          *lnTe**4+0.0001073294010367247d0*lnTe**5-8.36671960467864d-6&
          *lnTe**6+2.238306228891639d-7*lnTe**7))
    end where

    !H + H+ -> H2+
    where(Tgas.LE.6.7d3)
      vec_k(:,12) = small + (1.85d-23&
          *Tgas**1.8d0)
    end where

    !H + H+ -> H2+
    where(Tgas.GT.6.7d3)
      vec_k(:,13) = small + (5.81d-16&
          *(Tgas&
          /5.62d4)**(-0.6657d0*log10(Tgas/5.62d4)))
    end where

    !H2+ + H -> H2 + H+
    vec_k(:,14) = small + (6.0d-10)

    !H2 + H+ -> H2+ + H
    where(Tgas.GT.3.48d3)
      vec_k(:,15) = small + (exp(-24.24914687731536d0+3.400824447095291d0&
          *lnTe-3.898003964650152d0*lnTe**2+2.045587822403071d0&
          *lnTe**3-0.5416182856220388d0*lnTe**4+0.0841077503763412d0&
          *lnTe**5-0.007879026154483455d0&
          *lnTe**6+0.0004138398421504563d0*lnTe**7-9.36345888928611d-6&
          *lnTe**8))
    end where

    !H2 + E -> H + H + E
    vec_k(:,16) = small + (5.6d-11&
        *exp(-102124.d0*invT)*Tgas**0.5d0)

    !H2 + H -> H + H + H
    vec_k(:,17) = small + (1.0670825d-10&
        *Te**2.012d0*exp(-4.463d0*invTe)&
        /(1.d0+0.2472d0*Te)**3.512d0)

    !H- + E -> H + E + E
    vec_k(:,18) = small + (exp(-18.01849334273d0+2.360852208681d0&
        *lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0&
        *lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0&
        *lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0&
        *lnTe**7-2.631285809207d-6*lnTe**8))

    !H- + H -> H + H + E
    where(Tgas.LE.1.16d3)
      vec_k(:,19) = small + (2.56d-9&
          *Te**1.78186d0)
    end where

    !H- + H -> H + H + E
    where(Tgas.GT.1.16d3)
      vec_k(:,20) = small + (exp(-20.37260896533324d0+1.139449335841631d0&
          *lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0&
          *lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0&
          *lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0&
          *lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8&
          *lnTe**9))
    end where

    !H- + H+ -> H + H
    vec_k(:,21) = small + (6.5d-9/sqrt(Te))

    !H- + H+ -> H2+ + E
    vec_k(:,22) = small + (1.d-8*Tgas**(-0.4d0))

    !H2+ + E -> H + H
    where(Tgas.LE.6.17d2)
      vec_k(:,23) = small + (1.d-8)
    end where

    !H2+ + E -> H + H
    where(Tgas.GT.6.17d2)
      vec_k(:,24) = small + (1.32d-6&
          *Tgas**(-0.76d0))
    end where

    !H2+ + H- -> H + H2
    vec_k(:,25) = small + (5.d-7*sqrt(1.d2*invT))

    !H + H + H -> H2 + H
    where(Tgas.LE.3d2)
      vec_k(:,26) = small + (1.3d-32&
          *(T32)**(-0.38d0))
    end where

    !H + H + H -> H2 + H
    where(Tgas.GT.3d2)
      vec_k(:,27) = small + (1.3d-32&
          *(T32)**(-1.00d0))
    end where

    !H2 + H + H -> H2 + H2
    where(Tgas.LE.3d2)
      vec_k(:,28) = small + (1.3d-32&
          *(T32)**(-0.38d0) &
          / 8.d0)
    end where

    !H2 + H + H -> H2 + H2
    where(Tgas.GT.3d2)
      vec_k(:,29) = small + (1.3d-32&
          *(T32)**(-1.00d0) &
          / 8.d0)
    end where

    !H+ + D -> H + D+
    where(Tgas.GE.5d1)
      vec_k(:,30) = small + (2.00d-10&
          *Tgas**(0.402d0)*exp(-37.1d0*invT)-3.31d-17&
          *Tgas**(1.48d0))
    end where

    !H + D+ -> H+ + D
    where(Tgas.GE.5d1)
      vec_k(:,31) = small + (2.06d-10&
          *Tgas**(0.396)*exp(-33.d0*invT)+2.03d-9&
          *Tgas**(-0.332))
    end where

    !H2 + D+ -> HD + H+
    vec_k(:,32) = small + (1.d-9*(0.417+0.846&
        *log10(Tgas)-0.137*(log10(Tgas))**2))

    !HD + H+ -> H2 + D+
    vec_k(:,33) = small + (1.0d-9 *exp(-4.57d2&
        *invT))

    !H2 + D -> HD + H
    where(Tgas.LE.2.d3)
      vec_k(:,34) = small + (10**(-56.4737+5.88886&
          *log10(Tgas)+7.19692*(log10(Tgas))**2+2.25069&
          *(log10(Tgas))**3-2.16903*(log10(Tgas))**4+0.317887&
          *(log10(Tgas))**5))
    end where

    !H2 + D -> HD + H
    where(Tgas.GT.2d3)
      vec_k(:,35) = small + (3.17d-10&
          *exp(-5207.*invT))
    end where

    !HD + H -> H2 + D
    where(Tgas.GT.2d2)
      vec_k(:,36) = small + (5.25d-11&
          *exp(-4430.*invT+1.739d5*(invT)**2))
    end where

    !D + H- -> HD + E
    vec_k(:,37) = small + (1.5d-9*(T32)**(-0.1d0))

    !D+ + E -> D
    vec_k(:,38) = small + (3.6d-12&
        *(Tgas&
        /300)**(-0.75d0))

    coe(:,:) = vec_k(:,:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(vec_k(:,i)<0d0.or.vec_k(:,i)>kmax) print *,i,vec_k(:,i)
    !   end do
    !end if

  end function coe





!    !Tgas is in K
!
!    !------------------------------------!!!!!!!!!!
!    Tgas(:) = 1d3
!    !------------------------------------!!!!!!!!!!
!
!    !maxn initialization can be removed and small can be
!    !replaced with a proper value according to the environment
!!    nmax = max(maxval(n(1:vec_nvector,1:nmols)),1d0)
!    small = 1d-40!/(nmax*nmax*nmax)
!
!    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
!    lnTe = log(Te) !ln of Te (#)
!    T32 = Tgas*0.0033333333333333335 !Tgas/(300 K) (#)
!    invT = 1.d0/Tgas !inverse of T (1/K)
!    invTe = 1.d0/Te !inverse of T (1/eV)
!    sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)
!    log10Te = log10(Tgas)
!
!
!    do i=1,this%nvector
!      vec_k(i,:) = small(i) !initialize coefficients
!    end do
!
!    !H + E -> H+ + E + E
!    vec_k(:,1) = small + (exp(-32.71396786d0+13.5365560d0&
!        *lnTe-5.73932875d0*(lnTe**2)+1.56315498d0&
!        *(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2&
!        *(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4&
!        *(lnTe**7)-2.03914985d-6*(lnTe**8)))
!
!    !H+ + E -> H
!    where(Tgas.LE.5.5d3)
!      vec_k(:,2) = small + (3.92d-13&
!          *invTe**0.6353d0)
!    end where
!
!    !H+ + E -> H
!    where(Tgas.GT.5.5d3)
!      vec_k(:,3) = small + (exp(-28.61303380689232d0-0.7241125657826851d0&
!          *lnTe-0.02026044731984691d0*lnTe**2-0.002380861877349834d0&
!          *lnTe**3-0.0003212605213188796d0&
!          *lnTe**4-0.00001421502914054107d0&
!          *lnTe**5+4.989108920299513d-6*lnTe**6+5.755614137575758d-7&
!          *lnTe**7-1.856767039775261d-8*lnTe**8-3.071135243196595d-9&
!          *lnTe**9))
!    end where
!
!    !HE + E -> HE+ + E + E
!    vec_k(:,4) = small + (dexp(-44.09864886d0+23.91596563d0&
!        *lnTe-10.7532302d0*(lnTe**2)+3.05803875d0&
!        *(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2&
!        *(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4&
!        *(lnTe**7)-3.64916141d-6*(lnTe**8)))
!
!    !HE+ + E -> HE
!    where(Tgas.LE.9.28d3)
!      vec_k(:,5) = small + (3.92d-13&
!          *invTe**0.6353d0)
!    end where
!
!    !HE+ + E -> HE
!    where(Tgas.GT.9.28d3)
!      vec_k(:,6) = small + (1.54d-9&
!          *(1.d0+0.3d0&
!          /exp(8.099328789667d0&
!          *invTe))&
!          /(exp(40.49664394833662d0*invTe)&
!          *Te**1.5d0)+3.92d-13&
!          /Te**0.6353d0)
!    end where
!
!!    vec_k(:,7) = small + (exp(-68.71040990212001d0+43.93347632635d0&
!!        *lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0&
!!        *lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0&
!!        *lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0&
!!        *lnTe**7-3.165581065665d-6*lnTe**8))
!    vec_k(:,7) = small + (exp(-68.71040990212001d0+43.93347632635d0&
!        *lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0&
!        *lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0&
!        *lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0&
!        *lnTe**7-3.165581065665d-6*lnTe**8))
!
!    !HE++ + E -> HE+
!    vec_k(:,8) = small + (3.36d-10/sqrTgas/(Tgas&
!        /1.d3)**0.2d0/(1+(Tgas/1.d6)**0.7d0))
!
!    !H + E -> H-
!    vec_k(:,9) = small + (6.77d-15*Te**0.8779d0)
!
!    !H- + H -> H2 + E
!    where(Tgas.LT.1160d0)
!      vec_k(:,10) = small + (1.43d-9)
!    end where
!
!    !H- + H -> H2 + E
!    where(Tgas.GT.1160d0)
!      vec_k(:,11) = small + (exp(-20.06913897587003d0+0.2289800603272916d0&
!          *lnTe+0.03599837721023835d0*lnTe**2-0.004555120027032095d0&
!          *lnTe**3-0.0003105115447124016d0&
!          *lnTe**4+0.0001073294010367247d0*lnTe**5-8.36671960467864d-6&
!          *lnTe**6+2.238306228891639d-7*lnTe**7))
!    end where
!
!    !H + H+ -> H2+
!    where(Tgas.LE.6.7d3)
!      vec_k(:,12) = small + (1.85d-23&
!          *Tgas**1.8d0)
!    end where
!
!    !H + H+ -> H2+
!    where(Tgas.GT.6.7d3)
!      vec_k(:,13) = small + (5.81d-16&
!          *(Tgas&
!          /5.62d4)**(-0.6657d0*(log10Te-log10(5.62d4))))
!    end where
!
!    !H2+ + H -> H2 + H+
!    vec_k(:,14) = small + (6.0d-10)
!
!    !H2 + H+ -> H2+ + H
!    where(Tgas.GT.3.48d3)
!      vec_k(:,15) = small + (exp(-24.24914687731536d0+3.400824447095291d0&
!          *lnTe-3.898003964650152d0*lnTe**2+2.045587822403071d0&
!          *lnTe**3-0.5416182856220388d0*lnTe**4+0.0841077503763412d0&
!          *lnTe**5-0.007879026154483455d0&
!          *lnTe**6+0.0004138398421504563d0*lnTe**7-9.36345888928611d-6&
!          *lnTe**8))
!    end where
!
!    !H2 + E -> H + H + E
!!    vec_k(:,16) = small + (5.6d-11&
!!        *exp(-102124.d0*invT)*Tgas**0.5d0)
!    vec_k(:,16) = small + (5.6d-11&
!        *exp(-102124.d0*invT)*Tgas**0.5d0)
!
!    !H2 + H -> H + H + H
!    vec_k(:,17) = small + (1.0670825d-10&
!        *Te**2.012d0*exp(-4.463d0*invTe)&
!        /(1.d0+0.2472d0*Te)**3.512d0)
!
!    !H- + E -> H + E + E
!    vec_k(:,18) = small + (exp(-18.01849334273d0+2.360852208681d0&
!        *lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0&
!        *lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0&
!        *lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0&
!        *lnTe**7-2.631285809207d-6*lnTe**8))
!
!    !H- + H -> H + H + E
!    where(Tgas.LE.1.16d3)
!      vec_k(:,19) = small + (2.56d-9&
!          *Te(:)**1.78186d0)
!    end where
!
!    !H- + H -> H + H + E
!    where(Tgas.GT.1.16d3)
!      vec_k(:,20) = small + (exp(-20.37260896533324d0+1.139449335841631d0&
!          *lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0&
!          *lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0&
!          *lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0&
!          *lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8&
!          *lnTe**9))
!    end where
!
!    !H- + H+ -> H + H
!    vec_k(:,21) = small + (6.5d-9/sqrt(Te))
!
!    !H- + H+ -> H2+ + E
!    vec_k(:,22) = small + (1.d-8*Tgas**(-0.4d0))
!
!    !H2+ + E -> H + H
!    where(Tgas.LE.6.17d2)
!      vec_k(:,23) = small + (1.d-8)
!    end where
!
!    !H2+ + E -> H + H
!    where(Tgas.GT.6.17d2)
!      vec_k(:,24) = small + (1.32d-6&
!          *Tgas**(-0.76d0))
!    end where
!
!    !H2+ + H- -> H + H2
!    vec_k(:,25) = small + (5.d-7*sqrt(1.d2*invT))
!
!   !H + H + H -> H2 + H
!    where(Tgas.LE.3d2)
!      vec_k(:,26) = small + (1.3d-32&
!          *(T32(:))**(-0.38d0))
!    end where
!
!    !H + H + H -> H2 + H
!    where(Tgas.GT.3d2)
!      vec_k(:,27) = small + (1.3d-32&
!          *(T32)**(-1.00d0))
!    end where
!
!    !H2 + H + H -> H2 + H2
!    where(Tgas.LE.3d2)
!      vec_k(:,28) = small + (1.3d-32&
!          *(T32)**(-0.38d0) &
!          / 8.d0)
!    end where
!
!    !H2 + H + H -> H2 + H2
!    where(Tgas.GT.3d2)
!      vec_k(:,29) = small + (1.3d-32&
!          *(T32)**(-1.00d0) &
!          / 8.d0)
!    end where
!
!    !H+ + D -> H + D+
!    where(Tgas.GE.5d1)
!      vec_k(:,30) = small + (2.00d-10&
!          *Tgas**(0.402d0)*exp(-37.1d0*invT)-3.31d-17&
!          *Tgas**(1.48d0))
!    end where
!
!    !H + D+ -> H+ + D
!    where(Tgas.GE.5d1)
!      vec_k(:,31) = small + (2.06d-10&
!          *Tgas**(0.396)*exp(-33.d0*invT)+2.03d-9&
!          *Tgas**(-0.332))
!    end where
!
!    !H2 + D+ -> HD + H+
!    vec_k(:,32) = small + (1.d-9*(0.417+0.846&
!        *log10Te-0.137*(log10Te)**2))
!
!    !HD + H+ -> H2 + D+
!    vec_k(:,33) = small + (1.0d-9 *exp(-4.57d2&
!        *invT))
!
!   !H2 + D -> HD + H
!    where(Tgas.LE.2.d3)
!      vec_k(:,34) = small + (10**(-56.4737+5.88886&
!          *log10Te+7.19692*(log10Te)**2+2.25069&
!          *(log10Te)**3-2.16903*(log10Te)**4+0.317887&
!          *(log10Te)**5))
!    end where
!
!    !H2 + D -> HD + H
!    where(Tgas.GT.2d3)
!      vec_k(:,35) = small + (3.17d-10&
!          *exp(-5207.*invT))
!    end where
!
!    !HD + H -> H2 + D
!    where(Tgas.GT.2d2)
!      vec_k(:,36) = small + (5.25d-11&
!          *exp(-4430.*invT+1.739d5*(invT)**2))
!    end where
!
!    !D + H- -> HD + E
!    vec_k(:,37) = small + (1.5d-9*(T32)**(-0.1d0))
!
!    !D+ + E -> D
!    vec_k(:,38) = small + (3.6d-12&
!        *(Tgas&
!        /300)**(-0.75d0))
!
!    coe = vec_k !set coefficients to return variable
!
!
!    !!uncomment below to check coefficient values
!    !kmax = 1d0
!    !if(maxval(k)>kmax.or.minval(k)<0d0) then
!    !   print *,"***************"
!    !   do i=1,size(k)
!    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
!    !   end do
!    !end if
!
!  end function coe



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
