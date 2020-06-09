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
!> \brief provides a cell-vectorized ODE-solver based on the BDF-scheme
!!
!>  \attention This file needs to be preprocessed by pre_odevec.py
module odevec_main
  implicit none

  type :: csc_matrix
    double precision, pointer, dimension(:,:) :: sdata
    integer, pointer, dimension(:) :: u_col_start, l_col_start, row_index
  end type

  type :: odevec
#ODEVEC_VECTORLENGTH
#ODEVEC_EQUATIONS
#ODEVEC_REACTIONS
#ODEVEC_MAXORDER
#ODEVEC_NNZ
    integer :: Order                        !> current order
    integer :: SuccessesWithoutUpdate       !> count successful steps
    integer :: LastJacobianUpdate           !> count till last jacobian update
    integer :: FailedErrorTests             !> count failed error tests
    integer :: FailedCorrections            !> count failed convergence tests
    logical :: UpdateJac                    !> whether to update on next step
    logical :: NeglectUpperOrder            !> no allowance for order increase
    logical :: FirstStep                    !> extra handling for first step
    logical :: DoPermute                    !> skips permutation region if false
    logical :: CheckNegatives               !> flag for checking negatives in y
    integer :: ErrorCode                    !> value to return error

    double precision :: RelativeTolerance   !> relative tolerance
    double precision :: AbsoluteTolerance   !> absolute tolerance
    double precision :: MinimalTimestep     !> minimal timestep
    double precision :: ConvergenceRate     !> rate of convergence
    double precision :: JacobianChange      !> change of Jacobian
    double precision :: MaximumStepChange   !> maximal next timestep
    double precision :: OldCoefficient      !> saves old coefficient
    double precision, allocatable, dimension(:,:) :: y          !> current solution
                                                                !> array | often y_NS(:,:,0)
    integer         , allocatable, dimension(:)   :: Perm       !> permutation vector
    double precision, allocatable, dimension(:,:) :: en         !> correction vector
    double precision, allocatable, dimension(:,:) :: en_old     !> old corr. vector
    double precision, allocatable, dimension(:,:) :: den        !> corrector of en
    double precision, allocatable, dimension(:,:) :: den_tmp    !> buffer array
    double precision, allocatable, dimension(:,:) :: inv_weight_2 !> error weight
    double precision, allocatable, dimension(:,:) :: rhs        !> right-hand-side
    double precision, allocatable, dimension(:,:) :: res        !> residuum
    double precision, allocatable, dimension(:,:) :: coeff      !> coefficient matrix
    double precision, allocatable, dimension(:,:) :: tautable   !> solution table for
                                                                !> precomputation
    double precision, allocatable, dimension(:,:,:) :: y_NS     !> Nordsieck history array
#ODEVEC_LU_MATRIX
#ODEVEC_LU_METHOD
#ODEVEC_PACKAGING
  end type

  ! two different interfaces, whether it is called with y or y_NS
  interface WeightedNorm
    module procedure WeightedNorm1, WeightedNorm2
  end interface

  interface LUDecompose
    module procedure LUDecompose_dense, LUDecompose_sparse
  end interface
  interface SolveLU
    module procedure SolveLU_dense, SolveLU_sparse
  end interface
  interface CalcJac
    module procedure CalcJac_dense, CalcJac_sparse
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
              this%Perm(this%neq), &
              this%en(this%nvector,this%neq), &
              this%en_old(this%nvector,this%neq), &
              this%den(this%nvector,this%neq), &
              this%den_tmp(this%nvector,this%neq), &
              this%inv_weight_2(this%nvector,this%neq), &
              this%rhs(this%nvector,this%neq), &
              this%res(this%nvector,this%neq), &
              this%coeff(6,0:6), &
              this%tautable(this%maxorder+1,0:this%maxorder+1), &
              this%y_NS(this%nvector,this%neq,0:this%maxorder+1), &
              STAT=err)
    if (err.ne.0) then
      print *, "ODEVEC: Memory allocation error. OdeVec could not be intialized."
      stop
    end if

#ODEVEC_ALLOCATE_LU

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

#ODEVEC_DT_MIN

    ! permutations
#ODEVEC_PERMUTATIONS
    if(all(this%Perm.eq.(/(i, i=1, this%neq)/))) then
      this%DoPermute = .false.
    else
      this%DoPermute = .true.
    end if

    this%CheckNegatives = .false.

    ! sparsity patterns
#ODEVEC_SET_LU_SPARSITY

  end subroutine


  !> Main subroutine to be called from outside
  subroutine SolveODE(this,time,dt,t_stop,y,GetRHS,GetJac,GetLU,Mask)
    implicit none
    type(odevec) :: this
    external          :: GetRHS,GetJac,GetLU
    double precision  :: time, t_stop, dt
    double precision, dimension(this%nvector,this%neq) :: y
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    logical :: start_solver
    intent(in)        :: t_stop
    intent(inout)     :: y, time, dt

    ! some initializations
    this%FirstStep   = .true.
    this%order       = 1
    this%SuccessesWithoutUpdate = 0
    start_solver = .true.
    this%ConvergenceRate = 0.7d0
    this%JacobianChange = 0.0
    this%OldCoefficient = 1.0
    this%NeglectUpperOrder = .false.
    this%den      = 0.0d0
    this%FailedErrorTests  = 0
    this%FailedCorrections   = 0
    this%UpdateJac = .true.
    this%MaximumStepChange = 10.0
    this%ErrorCode = 0

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(this,this%RelativeTolerance,this%AbsoluteTolerance,y,this%inv_weight_2)

    if (present(Mask)) then
      if (.not.any(Mask)) then
        ! do nothing
        time = t_stop
        start_solver = .false.
      else
        start_solver = .true.
      end if
   end if

   if(start_solver) then
      this%y_NS(:,:,0) = y(:,:)

      ! Calculate initial right hand side and step-size
      call GetRHS(this,y,this%rhs)
      this%y_NS(:,:,1) = dt*this%rhs(:,:)

      call CalcInitialStepSize(this,time,t_stop,y,dt,Mask)

      ! main solve - solve the linear system
      do while (time < t_stop)
        ! solve the system
        call SolveLinearSystem(this,time,dt,y,GetRHS,GetJac,GetLU,Mask)
        if (this%ErrorCode.ne.0) return
      end do

      call InterpolateSolution(this,time,dt,t_stop,this%order,this%y_NS,y,Mask)
    end if
  end subroutine

  subroutine CalcInitialStepSize(this,time,t_stop,y,dt_init,Mask)
    implicit none
    type(odevec) :: this
    double precision :: time, t_stop, dt_init, tol, W0
    double precision, dimension(this%nvector,this%neq) :: W_inv2, y
    logical, intent(in), optional, dimension(this%nvector) :: Mask

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(this,this%RelativeTolerance,this%AbsoluteTolerance,y,this%inv_weight_2)

    tol = this%RelativeTolerance
    W_inv2 = this%inv_weight_2*tol*tol
    W0 = max(abs(time),abs(t_stop))


    if(present(Mask)) then
      dt_init = sqrt(tol/(W0**(-2.0) + maxval(sum(this%rhs(:,:)**2.0*W_inv2(:,:),DIM=2),MASK=Mask)/this%neq))
    else
      dt_init = sqrt(tol/(W0**(-2.0) + maxval(sum(this%rhs(:,:)**2.0*W_inv2(:,:),DIM=2))/this%neq))
    end if

    dt_init = min(dt_init,abs(t_stop-time))
  end subroutine


  !> Solves a linear system for a given timestep
  subroutine SolveLinearSystem(this,t,dt,y,GetRHS,GetJac,GetLU,Mask)
    implicit none
    type(odevec) :: this
    external          :: GetRHS,GetJac,GetLU
    double precision, dimension(this%nvector,this%neq) :: y
    double precision :: dt, t, dt_scale, told
    integer          :: i,j,k
    double precision :: error
    integer          :: CorrectorIterations
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    logical          :: success, ConvergenceFailed, Converged
    intent(inout)    :: y, dt, t
    ! advance in time
    told = t

    ! 2. predictor ------------------------------------------------------------!
    success = .false.
    predictor: do while(.not. success)
      t = t + dt

      ! predictor step - set predicted y_NS, en=0.0
      call PredictSolution(this,this%y_NS,this%en)
      if (present(Mask)) then
        where (spread(Mask,2,this%neq))
          y(:,:) = this%y_NS(:,:,0)
        end where
      else
        y(:,:) = this%y_NS(:,:,0)
      end if

      if (this%UpdateJac) then
        if (this%LUmethod.eq.1) then
          call GetJac(this,this%coeff(this%order,0),y,dt,this%LU)
          call LUDecompose(this,this%LU)
        else if (this%LUmethod.eq.2) then
          call GetLU(this,this%coeff(this%order,0),y,dt,this%LU)
        else if (this%LUmethod.eq.3) then
          call CalcJac(this,this%coeff(this%order,0),GetRHS,y,dt,this%LU,Mask)
          call LUDecompose(this,this%LU)
        end if
        this%UpdateJac = .false.
        this%LastJacobianUpdate = 0
        this%ConvergenceRate = 0.7d0
        this%JacobianChange = 1.0
      end if

      ! 3. corrector ----------------------------------------------------------!
      ConvergenceFailed = .false.
      Converged = .false.
      CorrectorIterations = 0
      corrector: do while (.not.Converged)

        call GetRHS(this,y,this%rhs)

        call CalcResiduum(this,dt,this%rhs,this%y_NS,this%en,this%res)

        ! calculates the solution for dy with given residuum
        if(this%DoPermute) then
          do i=1,this%neq
            this%den(:,i) = this%res(:,this%Perm(i))
          end do
          this%res(:,:) = this%den(:,:)
        end if

        call SolveLU(this,this%LU,this%res,this%den)

        if(this%DoPermute) then
          do i=1,this%neq
            this%res(:,this%Perm(i)) = this%den(:,i)
          end do
          this%den(:,:) = this%res(:,:)
        end if

        ! add correction to solution vector
        this%en(:,:) = this%en(:,:) + this%den(:,:)
        if (present(Mask)) then
          where (spread(Mask,2,this%neq))
            y(:,:) = this%y_NS(:,:,0) + this%coeff(this%order,0)*this%en
          end where
        else
          y(:,:) = this%y_NS(:,:,0) + this%coeff(this%order,0)*this%en
        end if

        if (this%CheckNegatives) then
          if(any(y(:,:).lt.0.0)) then
            this%ErrorCode = 1
          end if
        end if
        if(this%ErrorCode.ne.0) return

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25
        ! times the step-size
        call CheckConvergence(this,this%den,this%inv_weight_2, &
                              CorrectorIterations,this%ConvergenceRate, &
                              ConvergenceFailed,Converged,Mask)

        if (ConvergenceFailed.and..not.Converged) then
          this%MaximumStepChange = 2.0
          dt_scale = 0.25
          t = told
          this%FailedCorrections = this%FailedCorrections + 1
          if (this%FailedCorrections .gt. 10) then
            print *, "ODEVEC: Convergence failed! Abortion after more than 10 iterations."
            stop
          end if
          CorrectorIterations = 0
          call ResetSystem(this,this%y_NS)
          call SetStepSize(this,dt_scale,this%y_NS,dt)
          if (dt .lt. this%MinimalTimestep) then
            print *, "ODEVEC: Convergence failed! Abortion, because timestep too small."
            stop
          end if
          this%UpdateJac = .true.
          cycle predictor
        end if
      end do corrector
      this%FailedCorrections = 0

      ! local truncation error test:
      ! checks if solution is good enough and, similar to the convergence
      ! test, rerun from predictor with adjusted step size
      error = WeightedNorm(this,this%en,this%inv_weight_2,Mask)/ &
              tau(this%order,this%order,this%coeff)
      call CalcErrorWeightInvSquared(this,this%RelativeTolerance,this%AbsoluteTolerance,y,this%inv_weight_2)
      if (error > 1d0) then
        if (this%FailedErrorTests.ge.3) then
          dt_scale = 0.1
          t = told
          this%order = 1

          y(:,:) = this%y_NS(:,:,0)
          call GetRHS(this,y,this%rhs)

          this%y_NS(:,:,1) = dt*this%rhs(:,:)

          this%UpdateJac = .true.
        else
          dt_scale = 0.2
          t = told

          call ResetSystem(this,this%y_NS)

          this%UpdateJac = .true.
          this%SuccessesWithoutUpdate = 0
          this%MaximumStepChange = 2.0
          this%NeglectUpperOrder = .true.
        end if

        call CalcStepSizeOrder(this,dt_scale,this%order,dt,Mask)

        call SetStepSize(this,dt_scale,this%y_NS,dt)

        this%FailedErrorTests = this%FailedErrorTests + 1
        cycle predictor
      else
        this%SuccessesWithoutUpdate = this%SuccessesWithoutUpdate + 1
        this%LastJacobianUpdate = this%LastJacobianUpdate + 1

        ! after a successfull run rewrite the history array
        call UpdateNordsieck(this,this%en,this%y_NS)

        if (this%SuccessesWithoutUpdate.eq.this%order + 1) then
          call CalcStepSizeOrder(this,dt_scale,this%order,dt,Mask)
          call SetStepSize(this,dt_scale,this%y_NS,dt)
        end if

        this%FailedErrorTests = 0

        success = .true.
      end if
    end do predictor

    if (this%LastJacobianUpdate.ge.20) then
      this%UpdateJac = .true.
    else if (abs(this%JacobianChange - 1.0) .gt. 0.3) then
      this%UpdateJac = .true.
    end if

    this%en_old = this%en
  end subroutine SolveLinearSystem


  !> Calculates the Jacobian numerically
  !!
  !! Taken from DRPEPJ in odepack.
  subroutine CalcJac_dense(this,beta,GetRHS,y,dt,jac,Mask)
    implicit none
    type(odevec)     :: this
    integer          :: j,k,i
    external         :: GetRHS
    double precision :: beta,srur, dt,uround
    double precision, dimension(this%nvector) :: deltay, r, fac, r0, ytmp
    double precision, dimension(this%nvector,this%neq) :: y,Drhs
    double precision, dimension(this%nvector,this%neq,this%neq) :: jac
    logical, intent(in), optional, dimension(this%nvector) :: Mask

    uround = 2.2204460492503131E-016
    srur = 1.4901161193847656E-008

    call GetRHS(this, y, this%rhs)
    fac(:) = WeightedNorm(this, this%rhs, this%inv_weight_2, Mask)
    r0 = 1d3*abs(dt)*uround*this%neq*fac
    where (r0 .EQ. 0d0)
      r0 = 1d0
    end where
    do k = 1,this%neq
      ytmp(:) = y(:,k)
      r(:) = MAX(srur*abs(y(:,k)),r0(:)/sqrt(this%inv_weight_2(:,k)))
      y(:,k) = y(:,k) + r(:)
      fac(:) = -dt*beta/r
      call GetRHS(this, y, Drhs)
      do j = 1,this%neq
        jac(:,j,k) = (Drhs(:,j) - this%rhs(:,j))*fac
      end do
      y(:,k) = ytmp(:)
    end do
    do j = 1,this%neq
      jac(:,j,j) = 1 + jac(:,j,j)
    end do

  end subroutine CalcJac_dense

  !> Calculates the Jacobian numerically
  subroutine CalcJac_sparse(this,beta,GetRHS,y,dt,jac,Mask)
    implicit none
    type(odevec)     :: this
    integer          :: j,k,i
    external         :: GetRHS
    double precision :: beta,srur, dt
    double precision, dimension(this%nvector) :: deltay, r, fac, r0
    double precision, dimension(this%nvector,this%neq) :: y,ytmp,Drhs
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    type(csc_matrix) :: jac


    print *, "The Method Flags LUmethod=3 is not supported with sparse matrices, yet. Aborting..."
    stop
  end subroutine CalcJac_sparse

  !> Resets the whole system to the state where the predictor starts
  subroutine ResetSystem(this,y_NS)
    implicit none
    type(odevec)     :: this
    integer          :: i,j,k,l
    double precision, dimension(this%nvector,this%neq,0:6), intent(inout) :: y_NS

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

  end subroutine

  !> Calculates the step-sizes and uses according orders
  subroutine CalcStepSizeOrder(this,dt_scale,order,dt,Mask)
    implicit none
    type(odevec)     :: this
    integer          :: dt_maxloc, order
    double precision :: dt_scale, dt_scale_same, dt_scale_up, dt_scale_down, dt_upper_limit, dt
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    intent(out)      :: dt_scale
    intent(inout)    :: order,this,dt

    ! calculate all estimated step sizes
    ! for lower order
    if (order.eq.1) then
      dt_scale_down = 0.0
    else
      dt_scale_down = 1.d0/(1.3*(WeightedNorm(this,this%y_NS,this%inv_weight_2,order,Mask)/ &
                            tau(order,order-1,this%coeff))**(1d0/order) + 1d-6)
    end if
    ! for same order
    dt_scale_same = 1.d0/(1.2*(WeightedNorm(this,this%en,this%inv_weight_2,Mask)/ &
                          tau(order,order,this%coeff))**(1d0/(order+1d0)) + 1d-6)
    ! for higher order
    if (order.eq.this%maxorder.or.this%NeglectUpperOrder) then
      dt_scale_up = 0d0
    else
      dt_scale_up   = 1.0d0/(1.4*(WeightedNorm(this,this%en-this%en_old,this%inv_weight_2,Mask)/ &
                            tau(order,order+1,this%coeff))**(1d0/(order+2d0)) + 1d-6)
    end if

    ! choose largest and search for location of largest
    dt_maxloc = maxloc((/ dt_scale_down, dt_scale_same, dt_scale_up /),DIM=1)
    dt_scale = maxval((/ dt_scale_down, dt_scale_same, dt_scale_up /))

    ! set new order
    if (dt_scale.ge.1.1) then
      if (dt_maxloc == 1) then
        if (order /= 1) then
          order = order - 1
        else
          print *, "Error in CalcStepSizeOrder: order and step-size not fitting together"
          stop
        end if
      else if (dt_maxloc == 2) then
        dt_scale = dt_scale_same
        order = order
      else if (dt_maxloc == 3) then
        if (order /= this%maxorder) then
          this%y_NS(:,:,order+1) = this%coeff(order,order)*this%en(:,:)/(order+1d0)
          order = order + 1
        else
          print *, "Error in CalcStepSizeOrder: order and step-size not fitting together"
          stop
        end if
      end if
    else
      dt_scale = dt_scale_same
      order = order
    end if

    ! maximum increasement in step-size
    if (this%FirstStep) then
      dt_upper_limit = 1d4
      this%FirstStep = .false.
    else
      dt_upper_limit = 10d0
    end if
    dt_scale = min(dt_upper_limit,dt_scale)

    this%JacobianChange = this%JacobianChange*dt_scale

    this%SuccessesWithoutUpdate = 0
    this%NeglectUpperOrder = .false.
  end subroutine


  !> Calculates the error for convergence
  subroutine CheckConvergence(this,den,inv_weight_2,CorrectorIterations,ConvergenceRate,ConvergenceFailed,Converged,Mask)
    implicit none
    type(odevec)     :: this
    integer          :: CorrectorIterations
    logical          :: ConvergenceFailed, Converged
    double precision :: ConvergenceRate, ConvergenceRate_tmp, ConvergenceCriterion, ConvergenceError_tmp
    double precision, save :: ConvergenceError
    double precision, dimension(this%nvector,this%neq) :: den, inv_weight_2
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    intent(in)       :: den, inv_weight_2
    intent(inout)    :: ConvergenceRate, ConvergenceFailed, CorrectorIterations, Converged

    ConvergenceError_tmp = WeightedNorm(this,den,inv_weight_2,Mask)
    if (CorrectorIterations.ne.0) then
      ConvergenceRate_tmp = ConvergenceError_tmp/ConvergenceError
      ConvergenceRate = max(0.2d0*ConvergenceRate,ConvergenceRate_tmp)
    end if

    ConvergenceCriterion = ConvergenceError_tmp*min(1d0,1.5d0*ConvergenceRate)/ &
                 (tau(this%order,this%order,this%coeff)/(2d0*(this%order + 2d0)))

    ! sanity checks - these are important in extreme cases were NaNs resuls from
    ! the predictor or following routines
    if (any(den/=den).or.any(inv_weight_2/=inv_weight_2)) then
      ConvergenceFailed = .true.
      return
    end if

    if (ConvergenceCriterion.lt.1.0) then
      Converged = .true.
      return
    end if

    if (CorrectorIterations.ge.2 .and. ConvergenceRate.gt.2d0*ConvergenceRate_tmp) then
      ConvergenceFailed = .true.
    else if (CorrectorIterations.eq.3) then
      ConvergenceFailed = .true.
    else
      ConvergenceFailed = .false.
      ConvergenceError = ConvergenceError_tmp
      CorrectorIterations = CorrectorIterations + 1
    end if

  end subroutine CheckConvergence


  !> Calculates the weighted norm (two different interfaces for performance)
  function WeightedNorm1(this,en,inv_weight_2,Mask)
    implicit none
    type(odevec)     :: this
    double precision :: WeightedNorm1
    double precision, dimension(this%nvector,this%neq) :: en, inv_weight_2
    logical, optional, dimension(this%nvector) :: Mask
    intent(in)       :: en,inv_weight_2,this,Mask

    if (present(Mask)) then
    WeightedNorm1 = maxval(sqrt(sum(en(:,:)*en(:,:)*inv_weight_2(:,:),DIM=2)/ &
                         (this%neq)),MASK=Mask)
    else
    WeightedNorm1 = maxval(sqrt(sum(en(:,:)*en(:,:)*inv_weight_2(:,:),DIM=2)/ &
                         (this%neq)))
    end if
  end function WeightedNorm1
  pure function WeightedNorm2(this,en,inv_weight_2,column,Mask)
    implicit none
    type(odevec)     :: this
    integer          :: column
    double precision :: WeightedNorm2
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: en
    double precision, dimension(this%nvector,this%neq) :: inv_weight_2
    logical, optional, dimension(this%nvector) :: Mask
    intent(in)       :: en,inv_weight_2,column,this, Mask

    if (present(Mask)) then
    WeightedNorm2 = maxval(sqrt(sum(en(:,:,column)*en(:,:,column)*inv_weight_2(:,:),DIM=2)/ &
                         (this%neq)),MASK=Mask)
    else
    WeightedNorm2 = maxval(sqrt(sum(en(:,:,column)*en(:,:,column)*inv_weight_2(:,:),DIM=2)/ &
                         (this%neq)))
    end if
  end function WeightedNorm2


  !> Calculates the error-weight of y for given tolerances
  subroutine CalcErrorWeightInvSquared(this,RelativeTolerance,AbsoluteTolerance,y,inv_weight_2)
    implicit none
    type(odevec)     :: this
    double precision :: RelativeTolerance, AbsoluteTolerance
    double precision, dimension(this%nvector,this%neq) :: y, inv_weight_2
    intent(in)       :: RelativeTolerance, AbsoluteTolerance, y
    intent(out)      :: inv_weight_2

    inv_weight_2(:,:) = 1./(RelativeTolerance*abs(y(:,:)) + AbsoluteTolerance)
    inv_weight_2(:,:) = inv_weight_2(:,:)*inv_weight_2(:,:)
  end subroutine CalcErrorWeightInvSquared


  !> Sets the step size and scales the Nordsieck-array accordingly
  subroutine SetStepSize(this,dt_scale,y_NS,dt)
    implicit none
    type(odevec)     :: this
    integer          :: i,j,k
    double precision :: dt_scale,dtt_scale,dt
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: y_NS
    intent(in)       :: dt_scale
    intent(inout)    :: y_NS,dt,this

    dtt_scale = 1.0
    do k = 1,this%order
      dtt_scale = dt_scale*dtt_scale
!NEC$ collapse
      do j = 1,this%neq
!NEC$ ivdep
        do i = 1,this%nvector
          y_NS(i,j,k) = y_NS(i,j,k)*dtt_scale
        end do
      end do
    end do

    dt = dt*dt_scale
    this%JacobianChange = this%JacobianChange*dt_scale
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
!NEC$ ivdep
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
!NEC$ ivdep
          do i = 1,this%nvector
            y_NS(i,j,k-1) = y_NS(i,j,k) + y_NS(i,j,k-1)
          end do
        end do
      end do
    end do

    ! set the correction vector to zero
    en(:,:) = 0.0
  end subroutine PredictSolution


  !> Solve the DENSE System \f$ LU x=r \f$ with residuum r (res) and return x (den)
  subroutine SolveLU_dense(this,LU,res,den)
    implicit none
    type(odevec)   :: this
    double precision, dimension(this%nvector,this%neq,this%neq) :: LU
    double precision, dimension(this%nvector,this%neq)          :: res, den
    integer        :: i,j,k
    intent(in)     :: LU,res
    intent(inout)    :: den

!NEC$ collapse
    do j=1,this%neq
!NEC$ ivdep
      do i=1,this%nvector
        this%den_tmp(i,j) = 0.0
      end do
    end do
    ! lower system
    do k = 1,this%neq
      do j = 1,k-1
!NEC$ ivdep
        do i = 1,this%nvector
          this%den_tmp(i,k) = this%den_tmp(i,k) + LU(i,k,j)*this%den_tmp(i,j)
        end do
      end do
!NEC$ ivdep
      do i = 1,this%nvector
        this%den_tmp(i,k) = res(i,k) - this%den_tmp(i,k)
        this%den_tmp(i,k) = this%den_tmp(i,k)
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
          den(i,k) = den(i,k) + LU(i,k,j)*den(i,j)
        end do
      end do
!NEC$ ivdep
      do i = 1,this%nvector
        den(i,k) = this%den_tmp(i,k) - den(i,k)
        den(i,k) = den(i,k)/LU(i,k,k)
      end do
    end do

  end subroutine SolveLU_dense


  !> Solve the SPARSE System \f$ LU x=r \f$ with residuum r (res) and return x (den)
  subroutine SolveLU_sparse(this,LU,res,den)
    implicit none
    type(odevec) :: this
    type(csc_matrix) :: LU
    double precision, dimension(this%nvector,this%neq) :: res, den
    double precision, dimension(this%nvector) :: mult
    integer :: i,j,k,kk
    intent(in) :: LU,res
    intent(out) :: den

    den = res
    do k = 1,this%neq
!NEC$ ivdep
      do i = 1,this%nvector
        mult(i) = den(i,k)/LU%sdata(i,LU%l_col_start(k))
        den(i,k) = mult(i)
      end do
      do kk = LU%l_col_start(k) + 1, LU%u_col_start(k+1) - 1
!NEC$ ivdep
        do i = 1,this%nvector
          den(i,LU%row_index(kk)) = den(i,LU%row_index(kk)) - mult(i)*LU%sdata(i,kk)
        end do
      end do
    end do

    do k = this%neq,1,-1
!NEC$ ivdep
      do i = 1,this%nvector
        mult(i) = den(i,k)
      end do
      do kk = LU%u_col_start(k), LU%l_col_start(k) - 1
!NEC$ ivdep
        do i = 1, this%nvector
          den(i,LU%row_index(kk)) = den(i,LU%row_index(kk)) - mult(i)*LU%sdata(i,kk)
        end do
      end do
    end do

  end subroutine SolveLU_sparse


  !> Get LU from Jacobian for the dense case (in-place transformation)
  !!
  !! Based on dgetrf from LAPACK. LU-decomposition with partial pivoting. Taken
  !! and modified from https://rosettacode.org/wiki/LU_decomposition#Fortran.
  subroutine LUDecompose_dense(this,A)
    implicit none
    type(odevec)   :: this
    integer        :: i, j, k, m, kmax
    double precision, dimension(this%nvector,this%neq,this%neq) :: A
    integer,          dimension(1) :: maxloc_ij
    intent (inout) :: A

    do k = 1,this%neq-1
      do j = k+1,this%neq
!NEC$ ivdep
        do i = 1,this%nvector
          A(i,j,k) = A(i,j,k) / A(i,k,k)
        end do
      end do
      do j = k+1,this%neq
        do m = k+1,this%neq
!NEC$ ivdep
          do i = 1,this%nvector
            A(i,m,j) = A(i,m,j) - A(i,m,k) * A(i,k,j)
          end do
        end do
      end do
    end do
  end subroutine LUDecompose_dense


  !> Get LU with sparse computations
  !!
  !! Based on Duff, Erisman and Reid (2017): "Direct Methods for Sparse Matrices"
  !! Especially, ch. 10.3, pp. 210
  subroutine LUDecompose_sparse(this,A)
    implicit none
    type(odevec) :: this
    integer        :: i, j, k, jj, kk
    type(csc_matrix) :: A
    double precision, dimension(this%nvector,this%neq) :: w !< temporary column
    double precision, dimension(this%nvector) :: alpha

    do k=1,this%neq
      ! scatter
      do kk = A%u_col_start(k), A%u_col_start(k+1) - 1
!NEC$ ivdep
        do i = 1,this%nvector
          w(i,A%row_index(kk)) =  A%sdata(i,kk)
        end do
      end do

      ! calculate LU
      do kk = A%u_col_start(k), A%l_col_start(k) - 1
        j = A%row_index(kk)
!NEC$ ivdep
        do i = 1,this%nvector
          alpha(i) = w(i,j)/A%sdata(i,A%l_col_start(j))
          w(i,j) = alpha(i)
        end do
        do jj = A%l_col_start(j)+1,A%u_col_start(j+1)-1
!NEC$ ivdep
          do i = 1,this%nvector
            w(i,A%row_index(jj)) = w(i,A%row_index(jj)) - alpha(i)*A%sdata(i,jj)
          end do
        end do
      end do

      ! gather
      do kk = A%u_col_start(k), A%u_col_start(k+1) - 1
!NEC$ ivdep
        do i = 1,this%nvector
          A%sdata(i,kk) = w(i,A%row_index(kk))
        end do
      end do
    end do

  end subroutine LUDecompose_sparse


  !> Calculates the residuum
  subroutine CalcResiduum(this,dt,rhs,y_NS,en,res)
    implicit none
    type(odevec)      :: this
    integer           :: i,j
    double precision, intent(in) :: dt
    double precision, dimension(this%nvector,this%neq), intent(in) :: rhs, en
    double precision, dimension(this%nvector,this%neq,0:6), intent(in) :: y_NS
    double precision, dimension(this%nvector,this%neq), intent(out) :: res

!NEC$ collapse
    do j=1,this%neq
      do i=1,this%nvector
        res(i,j) = dt*rhs(i,j) - (y_NS(i,j,1) + en(i,j))
      end do
    end do
  end subroutine CalcResiduum


  !> Interpolates solution at output time
  !!
  !! If an output is requested because \$ t > t_{\mathrm{out}} \$ the solution
  !! is interpolated at the chosen output position. The interpolation is done
  !! by Taylor series expansion, where Horner's rule is applied for efficiency.
  subroutine InterpolateSolution(this,t,dt,t_out,order,y_NS,y,Mask)
    implicit none
    type(odevec)    :: this
    double precision, dimension(this%nvector,this%neq,0:this%maxorder+1) :: y_NS
    double precision, dimension(this%nvector,this%neq) :: y
    logical, optional, dimension(this%nvector) :: Mask
    double precision  :: t, t_out, rel_t, dt
    integer           :: k, k_dot, order
    intent(in)        :: dt, t_out, order, y_NS
    intent(inout)     :: t, y

    ! Horner's rule
    if (present(Mask)) then
      where (spread(Mask,2,this%neq))
        y(:,:) = y_NS(:,:,order)
      end where
    else
      y(:,:) = y_NS(:,:,order)
    end if
    rel_t = (t_out - t)/(dt)
    do k = 1,order
      k_dot = order - k
      if (present(Mask)) then
        where (spread(Mask,2,this%neq))
          y(:,:) = y_NS(:,:,k_dot) + rel_t*y(:,:)
        end where
      else
        y(:,:) = y_NS(:,:,k_dot) + rel_t*y(:,:)
      end if
    end do

    t = t_out
  end subroutine


  !> taufunction needed in several tests
  function tau(order,order2,coeff)
    implicit none
    integer :: i, order, order2
    double precision :: tau, faculty
    double precision, dimension(1:6,0:6) :: coeff
    intent(in) :: coeff, order

    faculty=1d0
!NEC$ shortloop
    do i=1,order
      faculty = faculty*i
    end do

    if (order-1.eq.order2) then
      tau = (order2+1d0)/(faculty)
    else
      tau = (order2+1d0)/(faculty*coeff(order,order))
    end if

    return
  end function tau


  !> deallocates all memory of the bdf-solver
  subroutine CloseOdevec(this)
    implicit none
    type(odevec) :: this
    integer :: err

    deallocate(this%y,this%en,this%en_old,this%den,this%den_tmp,this%inv_weight_2, &
               this%rhs,this%Perm,this%res,this%coeff,this%tautable, &
               this%y_NS, &
               stat=err)

#ODEVEC_DEALLOCATE_LU

    if (err.ne.0) then
      print *, "ODEVEC: Memory deallocation error. OdeVec was not properly closed."
      stop
    end if
  end subroutine




end module odevec_main
