!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   wasylenko :     Wasylenko's network of simplified HH-style neurons
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

!      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION VTC(NDIM/4), hTC(NDIM/4), VRE(NDIM/4), hRE(NDIM/4), s
      DOUBLE PRECISION sinfVTC(NDIM/4), dVTCdt(NDIM/4), dhTCdt(NDIM/4), dVREdt(NDIM/4), dhREdt(NDIM/4)
      INTEGER N
!      DOUBLE PRECISION, dimension(NDIM/4) :: minf_(NDIM/4), Hinf_(NDIM/4), sinf_(NDIM/4), tauinf_(NDIM/4), theSum_(NDIM/4)
      DOUBLE PRECISION, dimension(NDIM/4) :: minf_, Hinf_, sinf_, tauinf_, theSum_

real this, sigs, thh, sigh, thm, sigm, thtau, sigtau, tau1, tau2
real gLTC, VLTC, VsynTC, gLRE, VLRE, epRE, VsynRE, gCa, Cm

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

      sinfVTC = sinf(VTC)

      dVTCdt = (-gLTC * (VTC - VLTC) &
                  -gCa * (minf(VTC)) ** 3 * hTC * (VTC - VCa) &
                  -gTC * sinf(VRE) * (VTC - VsynRE) &
                 ) / Cm
      dhTCdt = epTC * (hinf(VTC) - hTC) / tauinf(VTC)
      dVREdt = (-gLRE * (VRE - VLRE) &
                  -gCa * (minf(VRE)) ** 3 * hRE * (VRE - VCa) &
                  -gRE / (2*omega+1) * theSum(sinfVTC) * (VRE - VsynTC) &
                 ) / Cm
      dhREdt = epRE * (hinf(VRE) - hRE) / tauinf(VRE)

      U(0*N:1*N) = dVTCdt
      U(1*N:2*N) = dhTRdt
      U(2*N:3*N) = dVREdt
      U(3*N:4*N) = dhREdt

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

        U(0*N:1*N) = -49  ! VTC_0
        U(1*N:2*N) = 0    ! hTC_0
        U(2*N:3*N) = -78  ! VRE_0
        U(3*N:4*N) = 0.5  ! hRE_0
        U(omega+3) = 0    ! a small perturbation in VTC_0

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

      DOUBLE PRECISION function theSum(sinfVTC)
      DOUBLE PRECISION, INTENT(IN) :: sinfVTC(*)
      INTEGER i, k, N, theIndex, theShape(*)
      DOUBLE PRECISION, DIMENSION(*) :: theSum
!        '''Used for computing the sum of inputs to each of the
!        RE neurons. Takes as argument a precomputed vector of s_inf(v_i)
!        values for each of the TC neurons' voltages, v_i.'''
        theShape = shape(sinfVTC)
        N = theShape(1)
        theSum = 0d0
        do i=1,N
            do k=-omega,omega+1
                theIndex = i + k
                if (theIndex >= N) then
                 theIndex = theIndex - N
                endif
                theSum(i) = theSum(i) + sinfVTC(theIndex)
            end do
        end do
        END SUBROUTINE

