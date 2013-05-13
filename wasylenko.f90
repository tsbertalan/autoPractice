!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   wasylenko :     Wasylenko's network of simplified HH-style neurons
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

!    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
    DOUBLE PRECISION, dimension(NDIM/4) :: VTC, hTC, VRE, hRE, dVTCdt, dhTCdt, dVREdt, dhREdt, theSum
    INTEGER N, i, k, omega, theIndex
    DOUBLE PRECISION s, minf, Hinf, sinf, tauinf

    DOUBLE PRECISION ths, sigs, thh, sigh, thm, sigm, thtau, sigtau, tau1, tau2
    DOUBLE PRECISION gLTC, VLTC, VsynTC, gLRE, VLRE, epRE, VsynRE, gCa, Cm, v

!# bifurcation parameters
    double precision epTC, gTC, gRE, VCa

    ths = -20
    sigs = 2
    thh = -79
    sigh = -5
    thm = -65
    sigm = 7.8
    thtau = -65
    sigtau = 4
    tau1 = 1
    tau2 = 80

    gLTC = .01
    VLTC = -75
    VsynTC = 0
    gLRE = .2
    VLRE = -80
    epRE = 2
    VsynRE = -80
    gCa = 1
    Cm = 1

    N = NDIM / 4

    VTC = U(0*N:1*N)
    hTC = U(1*N:2*N)
    VRE = U(2*N:3*N)
    hRE = U(3*N:4*N)

    s=PAR(1)
    epTC = 1 + 2 * s
    gTC = 0.03 + 0.07 * s
    gRE = 0.1 + 0.2 * s
    VCa = 120 - 30 * s

    omega = 6

    theSum = 0d0
    do i=1,N
        do k=-omega,omega+1
            theIndex = i + k
            if (theIndex >= N) then
             theIndex = theIndex - N
            endif
            theSum(i) = theSum(i) + sinf(VTC(theIndex))
    !                theSum(i) = theSum(i) + 1
        end do
    end do

    !      sinfVTC = sinf(VTC)
    do i=1,N
      dVTCdt(i) = (-gLTC * (VTC(i) - VLTC) &
                  -gCa * (minf(VTC(i))) ** 3 * hTC(i) * (VTC(i) - VCa) &
                  -gTC * sinf(VRE(i)) * (VTC(i) - VsynRE) &
                  ) / Cm
      dhTCdt(i) = epTC * (hinf(VTC(i)) - hTC(i)) / tauinf(VTC(i))
      dVREdt(i) = (-gLRE * (VRE(i) - VLRE) &
                  -gCa * (minf(VRE(i))) ** 3 * hRE(i) * (VRE(i) - VCa) &
                  -gRE / (2*omega+1) * theSum(i) * (VRE(i) - VsynTC) &
                 ) / Cm
      dhREdt(i) = epRE * (hinf(VRE(i)) - hRE(i)) / tauinf(VRE(i))
    end do

    do i=1,N
        F(0*N + i) = dVTCdt(i)
        F(1*N + i) = dhTRdt(i)
        F(2*N + i) = dVREdt(i)
        F(3*N + i) = dhREdt(i)
    end do

!        # Wasylenko 2010, Sec. 3.2:
!        #  First, we started with a perturbation in the vTC from the (spatially
!        #  uniform) rest state. This generated both left- and right-traveling
!        #  waves, but by transiently changing the boundary conditions it was
!        #  possible to eliminate the left-travelingwave.
!        # Here is that transient modification. It's possible that I could do it better:
!I DON'T KNOW HOW TO DO THIS IN FORTRAN FOR AUTO
!        if t<600:
!            dVTCdt[0:omega] = 0
!            dhTCdt[0:omega] = 0
!            dVREdt[0:omega] = 0
!            dhREdt[0:omega] = 0

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

       PAR(1)=0.5

        U(0*NDIM/4:1*NDIM/4) = -49  ! VTC_0
        U(1*NDIM/4:2*NDIM/4) = 0    ! hTC_0
        U(2*NDIM/4:3*NDIM/4) = -78  ! VRE_0
        U(3*NDIM/4:4*NDIM/4) = 0.5  ! hRE_0
        U(6+3) = 0    ! a small perturbation in VTC_0

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

      DOUBLE PRECISION function minf(v)
      DOUBLE PRECISION, INTENT(IN) :: v
        minf = 1 / (1 + exp(-(v - thm) / sigm))
      END FUNCTION
    
      DOUBLE PRECISION function hinf(v)
      DOUBLE PRECISION, INTENT(IN) :: v
        hinf = 1 / (1 + exp(-(v - thh) / sigh))
      END FUNCTION
    
      DOUBLE PRECISION function sinf(v)
      DOUBLE PRECISION, INTENT(IN) :: v
        sinf = 1 / (1 + exp(-(v - ths) / sigs))
      END FUNCTION
     
      DOUBLE PRECISION function tauinf(v)
      DOUBLE PRECISION, INTENT(IN) :: v
        tauinf = tau1 + (tau2 - tau1) / (1 + exp(-(v - thtau) / sigtau))
      END FUNCTION


