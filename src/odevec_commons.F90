module odevec_commons
  use odevec_main
  implicit none

  contains

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
#ifdef HAVE_COEFFICIENTS
    k = 0.0d0
    k = coe(this)
#endif

!cdir nodep
    do i=1,this%nvector

#ODEVEC_JAC

    end do

  end subroutine GetJac


  !> Provides the Matrices L and U
  subroutine GetLU(this,beta,y,dt,L,U)
    implicit none
    type(odevec), intent(in) :: this
    double precision, dimension(this%nvector,this%neq), &
         intent(in)  :: y
    double precision, dimension(this%nvector,this%neq,this%neq), &
         intent(out) :: L,U
    double precision, dimension(this%nvector,this%nrea) :: k
    double precision :: beta, dt
    integer          :: i

#ifdef HAVE_COEFFICIENTS
    k = 0.0d0
    k = coe(this)
#endif

!cdir nodep
    do i=1,this%nvector

#ODEVEC_L

#ODEVEC_U

    end do
  end subroutine GetLU


  !> Example case for the right-hand-side
  !!
  !! This is for Robertson's example
  subroutine GetRHS(this,y,rhs)
    implicit none
    type(odevec)    :: this
    integer         :: i
    double precision, dimension(this%nvector,this%neq) :: y
    double precision, dimension(this%nvector,this%neq) :: rhs
    double precision, dimension(this%nvector,this%nrea) :: k
    intent (in)     :: y
    intent (out)    :: rhs

#ifdef HAVE_COEFFIENTS
    k = 0.0d0
    k = coe(this
#endif

!cdir nodep
    do i=1,this%nvector
#ODEVEC_RHS
    end do
  end subroutine GetRHS


#ifdef HAVE_COEFFICIENTS
  ! copied from krome as a test
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

    Tgas=1d5

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
#endif

end module odevec_commons


