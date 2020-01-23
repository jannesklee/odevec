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
    integer :: iterator                                     !> iterator
    integer :: order                                        !> current order
    integer :: dtorder_count                                !> last step-size/order change
    logical :: UpdateJac
    integer :: UpdateJac_count
    logical :: FirstStep
    logical :: update_dtorder
    logical :: check_negatives = .false.

    double precision :: rtol                                !> relative tolerance
    double precision :: atol                                !> absolute tolerance
    double precision :: dt_min                              !> minimal timestep
    double precision, pointer, dimension(:,:) :: y          !> current solution
                                                            !> array | often y_NS(:,:,0)
    integer         , pointer, dimension(:,:) :: Piv        !> pivoting vector
    integer         , pointer, dimension(:)   :: Perm       !> permutation vector
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
    double precision, pointer, dimension(:,:,:) :: y_NS     !> Nordsieck history array
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
              this%Piv(this%nvector,this%neq), &
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
    this%iterator    = 0
    this%dtorder_count = 0
    this%update_dtorder = .true.
    this%UpdateJac_count = 0
    start_solver = .true.

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
      if (present(Mask)) then
        CALL CalcInitialStepSize(this,time,t_stop,y,dt,Mask)
      else
        CALL CalcInitialStepSize(this,time,t_stop,y,dt)
      end if
      this%y_NS(:,:,1) = dt*this%rhs(:,:)

      ! main solve - solve the linear system
      do while (time < t_stop)
        ! solve the system
        if (present(Mask)) then
          call SolveLinearSystem(this,time,dt,y,GetRHS,GetJac,GetLU,Mask)
        else
          call SolveLinearSystem(this,time,dt,y,GetRHS,GetJac,GetLU)
        end if

        this%iterator = this%iterator + 1
      end do
      if (present(Mask)) then
        call InterpolateSolution(this,time,dt,t_stop,this%order,this%y_NS,y,Mask)
      else
        call InterpolateSolution(this,time,dt,t_stop,this%order,this%y_NS,y)
      end if
    end if
  end subroutine

  subroutine CalcInitialStepSize(this,time,t_stop,y,dt_init,Mask)
    implicit none
    type(odevec) :: this
    double precision :: time, t_stop, dt_init, tol, W0
    double precision, dimension(this%nvector,this%neq) :: W_inv2, y
    logical, intent(in), optional, dimension(this%nvector) :: Mask

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(this,this%rtol,this%atol,y,this%inv_weight_2)

    tol = this%rtol
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
    double precision :: conv_error, conv_rate, ConvergenceCriterion, error
    integer          :: CorrectorIterations, conv_failed, lte_iterator
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    logical          :: success, ConvergenceFailed
    intent(inout)    :: y, dt, t

    ! 1. initialization -------------------------------------------------------!
    ! use the LU matrices from the predictor value
    ! some initializations
    this%den      = 0.0d0
    conv_rate     = 0.7d0
    lte_iterator  = 0
    conv_failed   = 0
    this%UpdateJac = .true.

    ! needed in order to calculate the weighted norm
    call CalcErrorWeightInvSquared(this,this%rtol,this%atol,y,this%inv_weight_2)

    ! advance in time
    told = t

    ! 2. predictor ------------------------------------------------------------!
    success = .false.
    predictor: do while(.not. success)
      ConvergenceFailed = .false.
      CorrectorIterations = 0

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
          call LUDecompose(this,this%LU,this%Piv)
        else if (this%LUmethod.eq.2) then
          call GetLU(this,this%coeff(this%order,0),y,dt,this%LU)
        else if (this%LUmethod.eq.3) then
          if (present(Mask)) then
            call CalcJac(this,this%coeff(this%order,0),GetRHS,y,dt,this%LU,Mask)
          else
            call CalcJac(this,this%coeff(this%order,0),GetRHS,y,dt,this%LU)
          end if
          call LUDecompose(this,this%LU,this%Piv)
        end if
        this%UpdateJac = .false.
        this%UpdateJac_count = 0
        conv_rate     = 0.7d0
      end if

      ! 3. corrector ----------------------------------------------------------!
      corrector: do while (ConvergenceCriterion.gt.1.0 .or. CorrectorIterations.eq.0)


        ! calculates residuum
        call CalcResiduum(this,GetRHS,y,dt,this%res)

        ! calculates the solution for dy with given residuum
        do i=1,this%neq
          this%den(:,i) = this%res(:,this%Perm(i))
        end do
        this%res(:,:) = this%den(:,:)

        call SolveLU(this,this%LU,this%Piv,this%res,this%den)

        do i=1,this%neq
          this%res(:,this%Perm(i)) = this%den(:,i)
        end do
        this%den(:,:) = this%res(:,:)


        ! add correction to solution vector
        this%en(:,:) = this%en(:,:) + this%den(:,:)
        if (present(Mask)) then
          where (spread(Mask,2,this%neq))
            y(:,:) = this%y_NS(:,:,0) + this%coeff(this%order,0)*this%en
          end where
        else
          y(:,:) = this%y_NS(:,:,0) + this%coeff(this%order,0)*this%en
        end if

        ! convergence test:
        ! if fail reset and run again starting at predictor step with 0.25
        ! times the step-size
        if (present(Mask)) then
          call CheckConvergence(this,CorrectorIterations,conv_rate,this%den, &
                                conv_error,this%inv_weight_2,ConvergenceFailed,ConvergenceCriterion,Mask)
        else
          call CheckConvergence(this,CorrectorIterations,conv_rate,this%den, &
                                conv_error,this%inv_weight_2,ConvergenceFailed,ConvergenceCriterion)
        end if

        if (this%check_negatives) then
          if(ANY(y.lt.-TINY(y))) ConvergenceFailed=.true.
        end if
        if(any(y(:,:)/=y(:,:))) ConvergenceFailed=.true.

        if (ConvergenceFailed) then
          dt_scale = 0.25
          t = told
          conv_failed = conv_failed + 1
          if (conv_failed .gt. 10) then
            print *, "ODEVEC: Convergence failed! Abortion after more than 10 iterations."
            stop
          end if
          CorrectorIterations = 0
          call ResetSystem(this,dt_scale,dt,this%y_NS)
          if (dt .lt. this%dt_min) then
            print *, "ODEVEC: Convergence failed! Abortion, because timestep too small."
            stop
          end if
          this%UpdateJac = .true.
          cycle predictor
        end if
      end do corrector
      conv_failed = 0

      ! local truncation error test:
      ! checks if solution is good enough and, similar to the convergence
      ! test, rerun from predictor with adjusted step size
      if (present(Mask)) then
        error = WeightedNorm(this,this%en,this%inv_weight_2,Mask)/ &
                tau(this%order,this%order,this%coeff)
      else
        error = WeightedNorm(this,this%en,this%inv_weight_2)/ &
                tau(this%order,this%order,this%coeff)
      end if
      call CalcErrorWeightInvSquared(this,this%rtol,this%atol,y,this%inv_weight_2)
      if (error > 1d0) then
        dt_scale = 0.2
        t = told
        call ResetSystem(this,dt_scale,dt,this%y_NS)
        lte_iterator = lte_iterator + 1
        this%UpdateJac = .true.
        this%update_dtorder = .true.
        this%dtorder_count = 0
        cycle predictor
      else
        this%dtorder_count = this%dtorder_count + 1
        this%UpdateJac_count = this%UpdateJac_count + 1
        success = .true.
      end if
    end do predictor

    ! after a successfull run rewrite the history array
    call UpdateNordsieck(this,this%en,this%y_NS)

    if (this%UpdateJac_count.ge.20) then
      this%UpdateJac = .true.
    end if

    ! preparation for next step
    ! 4. step-size/order control ----------------------------------------------!
    ! calc. step size for order+(0,-1,+1) -> use largest for next step
    if (this%update_dtorder.or.(this%dtorder_count.eq.this%order + 1)) then
      if (present(Mask)) then
        call CalcStepSizeOrder(this,dt_scale,this%order,dt,Mask)
      else
        call CalcStepSizeOrder(this,dt_scale,this%order,dt)
      end if
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
    double precision, dimension(this%nvector) :: deltay, r, fac, r0
    double precision, dimension(this%nvector,this%neq) :: y,ytmp,Drhs
    double precision, dimension(this%nvector,this%neq,this%neq) :: jac
    logical, intent(in), optional, dimension(this%nvector) :: Mask

    uround = 2.2204460492503131E-016
    srur = 1.4901161193847656E-008

    call GetRHS(this, y, this%rhs)
    if (present(Mask)) then
      fac(:) = WeightedNorm(this, this%rhs, this%inv_weight_2, Mask)
    else
      fac(:) = WeightedNorm(this, this%rhs, this%inv_weight_2)
    end if
    r0 = 1d3*abs(dt)*uround*this%neq*fac
    where (r0 .EQ. 0d0)
      r0 = 1d0
    end where
    do k = 1,this%neq
      ytmp(:,k) = y(:,k)
      r(:) = MAX(srur*abs(y(:,k)),r0(:)/this%inv_weight_2(:,k))
      y(:,k) = y(:,k) + r(:)
      fac(:) = -dt*beta/r
      call GetRHS(this, y, Drhs)
      do j = 1,this%neq
        jac(:,j,k) = (Drhs(:,j) - this%rhs(:,j))*fac
      end do
      y(:,k) = ytmp(:,k)
    end do
    do j = 1,this%neq
      jac(:,j,j) = 1 + jac(:,j,j)
    end do

  end subroutine CalcJac_dense

  !> Calculates the Jacobian numerically
  subroutine CalcJac_sparse(this,beta,GetRHS,y,dt,jac)
    implicit none
    type(odevec)     :: this
    integer          :: j,k,i
    external         :: GetRHS
    double precision :: beta,srur, dt
    double precision, dimension(this%nvector) :: deltay, r, fac, r0
    double precision, dimension(this%nvector,this%neq) :: y,ytmp,Drhs
    type(csc_matrix) :: jac


    print *, "The Method Flags LUmethod=3 is not supported with sparse matrices, yet. Aborting..."
    stop
  end subroutine CalcJac_sparse

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
    if (present(Mask)) then
      dt_scale_down = 1.d0/(1.3*(WeightedNorm(this,this%y_NS,this%inv_weight_2,order,Mask)/ &
                            tau(order,order-1,this%coeff))**(1d0/order) + 1d-6)
    else
      dt_scale_down = 1.d0/(1.3*(WeightedNorm(this,this%y_NS,this%inv_weight_2,order)/ &
                            tau(order,order-1,this%coeff))**(1d0/order) + 1d-6)
    end if
    end if
    ! for same order
    if (present(Mask)) then
      dt_scale_same = 1.d0/(1.2*(WeightedNorm(this,this%en,this%inv_weight_2,Mask)/ &
                            tau(order,order,this%coeff))**(1d0/(order+1d0)) + 1d-6)
    else
      dt_scale_same = 1.d0/(1.2*(WeightedNorm(this,this%en,this%inv_weight_2)/ &
                            tau(order,order,this%coeff))**(1d0/(order+1d0)) + 1d-6)
    end if
    ! for higher order
    if (order.eq.this%maxorder) then
      dt_scale_up = 0d0
    else
    if (present(Mask)) then
      dt_scale_up   = 1.0d0/(1.4*(WeightedNorm(this,this%en-this%en_old,this%inv_weight_2,Mask)/ &
                            tau(order,order+1,this%coeff))**(1d0/(order+2d0)) + 1d-6)
    else
      dt_scale_up   = 1.0d0/(1.4*(WeightedNorm(this,this%en-this%en_old,this%inv_weight_2)/ &
                            tau(order,order+1,this%coeff))**(1d0/(order+2d0)) + 1d-6)
    end if
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


    ! Adjust the Nordsieck history array with new step size & new order
    call SetStepSize(this,dt_scale,this%y_NS,dt)

    this%update_dtorder = .false.
    this%dtorder_count = 0
  end subroutine


  !> Calculates the error for convergence
  subroutine CheckConvergence(this,CorrectorIterations,conv_rate,den,conv_error,inv_weight_2,ConvergenceFailed,ConvergenceCriterion,Mask)
    implicit none
    type(odevec)     :: this
    integer          :: CorrectorIterations
    logical          :: ConvergenceFailed
    double precision :: conv_rate, conv_rate2, conv_error_tmp, conv_error, ConvergenceCriterion
    double precision, dimension(this%nvector,this%neq) :: den, inv_weight_2
    logical, intent(in), optional, dimension(this%nvector) :: Mask
    intent(in)       :: den, inv_weight_2
    intent(inout)    :: conv_rate, conv_error, ConvergenceFailed, CorrectorIterations, ConvergenceCriterion

    if (present(Mask)) then
      conv_error_tmp = WeightedNorm(this,den,inv_weight_2,Mask)
    else
      conv_error_tmp = WeightedNorm(this,den,inv_weight_2)
    end if
    if (CorrectorIterations.ne.0) then
      conv_rate2 = conv_error_tmp/conv_error
      conv_rate = max(0.2d0*conv_rate,conv_rate2)
    end if
    ConvergenceCriterion = conv_error_tmp*min(1d0,1.5d0*conv_rate)/ &
                 (tau(this%order,this%order,this%coeff)/(2d0*(this%order + 2d0)))


    if (CorrectorIterations.ge.2 .and. conv_rate.gt.2d0*conv_rate2) then
      ConvergenceFailed = .true.
    else if (CorrectorIterations.eq.3) then
      ConvergenceFailed = .true.
    else
      ConvergenceFailed = .false.
    end if

    conv_error = conv_error_tmp
    CorrectorIterations = CorrectorIterations + 1
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
!NEC$ ivdep
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
  subroutine SolveLU_dense(this,LU,Piv,res,den)
    implicit none
    type(odevec)   :: this
    double precision, dimension(this%nvector,this%neq,this%neq) :: LU
    double precision, dimension(this%nvector,this%neq)          :: res, den
    integer         , dimension(this%nvector,this%neq)          :: Piv
    integer        :: i,j,k
    intent(in)     :: LU,Piv,res
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
!          this%den_tmp(i,k) = this%den_tmp(i,k) + LU(i,Piv(i,k),j)*this%den_tmp(i,j)
        end do
      end do
!NEC$ ivdep
      do i = 1,this%nvector
        this%den_tmp(i,k) = res(i,k) - this%den_tmp(i,k)
        this%den_tmp(i,k) = this%den_tmp(i,k)
!        this%den_tmp(i,k) = res(i,Piv(i,k)) - this%den_tmp(i,k)
!        this%den_tmp(i,k) = this%den_tmp(i,k)
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
!          den(i,k) = den(i,k) + LU(i,Piv(i,k),j)*den(i,j)
        end do
      end do
!NEC$ ivdep
      do i = 1,this%nvector
        den(i,k) = this%den_tmp(i,k) - den(i,k)
        den(i,k) = den(i,k)/LU(i,k,k)
!        den(i,k) = this%den_tmp(i,k) - den(i,k)
!        den(i,k) = den(i,k)/LU(i,Piv(i,k),k)
      end do
    end do

  end subroutine SolveLU_dense


  !> Solve the SPARSE System \f$ LU x=r \f$ with residuum r (res) and return x (den)
  subroutine SolveLU_sparse(this,LU,Piv,res,den)
    implicit none
    type(odevec) :: this
    type(csc_matrix) :: LU
    double precision, dimension(this%nvector,this%neq) :: res, den
    integer, dimension(this%neq) :: Piv
    double precision, dimension(this%nvector) :: mult
    integer :: i,j,k,kk
    intent(in) :: LU,Piv,res
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
  subroutine LUDecompose_dense(this,A,P)
    implicit none
    type(odevec)   :: this
    integer        :: i, j, k, m, kmax
    double precision, dimension(this%nvector,this%neq,this%neq) :: A
    integer,          dimension(this%nvector,this%neq) :: P
    integer,          dimension(1) :: maxloc_ij!,kmax
    intent (inout) :: A
    intent (out)   :: p

    ! initialize P
!    do i=1,this%nvector
!      P(i,:) = [(j, j=1, this%neq)]
!    end do
!    do k = 1,this%neq-1
!      do i=1,this%nvector
!        maxloc_ij = maxloc(abs(A(i,P(i,k:),k)))
!        kmax = maxloc_ij(1) + k - 1
!        if (kmax /= k ) then
!          P(i,[k, kmax]) = P(i,[kmax, k])
!        end if
!      end do
!    end do

    do k = 1,this%neq-1
      do j = k+1,this%neq
!NEC$ ivdep
        do i = 1,this%nvector
!          A(i,P(i,j),k) = A(i,P(i,j),k) / A(i,P(i,k),k)
          A(i,j,k) = A(i,j,k) / A(i,k,k)
        end do
      end do
      do j = k+1,this%neq
        do m = k+1,this%neq
!NEC$ ivdep
          do i = 1,this%nvector
!            A(i,P(i,m),j) = A(i,P(i,m),j) - A(i,P(i,m),k) * A(i,P(i,k),j)
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
  subroutine LUDecompose_sparse(this,A,P)
    implicit none
    type(odevec) :: this
    integer        :: i, j, k, jj, kk
    type(csc_matrix) :: A
    double precision, dimension(this%nvector,this%neq) :: w !< temporary column
    double precision, dimension(this%nvector) :: alpha
    integer,          dimension(this%neq) :: P

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
               this%rhs,this%Piv,this%Perm,this%res,this%coeff,this%tautable, &
               this%y_NS, &
               stat=err)

#ODEVEC_DEALLOCATE_LU

    if (err.ne.0) then
      print *, "ODEVEC: Memory deallocation error. OdeVec was not properly closed."
      stop
    end if
  end subroutine




end module odevec_main
