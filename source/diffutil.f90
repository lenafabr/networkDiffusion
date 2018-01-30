MODULE DIFFUTIL
  ! utilities for sampling diffusive motion of particles

  IMPLICIT NONE
CONTAINS

  SUBROUTINE SAMPLEFPT2ABS(X0,L,D,FPT,WHICHLEAVE)
    ! sample the first passage time for a particle starting at position X0
    ! between two absorbing boundaries (at 0 and L)
    ! particle diffusivity D
    ! WHICHLEAVE = 0 if leaving on the left, 1 if leaving on the right
    USE GENUTIL, ONLY : PI
    USE mt19937, ONLY : GRND
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X0, L, D
    DOUBLE PRECISION, INTENT(OUT) :: FPT
    INTEGER, INTENT(OUT) :: WHICHLEAVE

    ! hard-coded maximum indices to avoid unrealistically large loops
    INTEGER, PARAMETER :: MAXN = 500, MAXK = 5000
    ! list of N values for sine series
    INTEGER :: NLIST(MAXN), NSET, KLIST(MAXK)
    ! coefficients for sine series
    DOUBLE PRECISION :: COEFF(MAXN)
    ! precalculate a coefficient for figuring out how
    ! many N indices are needed for sin series convergence
    DOUBLE PRECISION :: COEFFCHECKN, COEFFCHECKK

    ! nondimensionalized initial position
    DOUBLE PRECISION :: X0L

    ! for testing purposes only
    INTEGER :: TC
    DOUBLE PRECISION :: T, H, TMP
    
    ! non-dimensionalize length units by L, time units by L^2/D
    X0L = X0/L

    ! decide in what direction to leave   
    TMP = GRND()      
    IF (TMP < X0L) THEN
       WHICHLEAVE = 1
    ELSE
       WHICHLEAVE = 0
    ENDIF    

    
    ! nset = for sine series, how many coefficients have already been calculated
    NSET = 1

    ! precalculated coefficient for establishing how many
    ! n indices are necessary for convergence
    COEFFCHECKN =  SQRT(-LOG(EPSILON(1D0))/PI**2*10)
    ! coefficient for how many k indices needed for convergence
    ! with image method summation
    COEFFCHECKK = SQRT(-LOG(EPSILON(1D0)))
    
    ! ----------
    ! For testing only: plot the cumulative FPT function
    ! ----------
    ! OPEN(UNIT=88,FILE='test.out')
    ! DO TC = 1,100
    !    T = 10**(-2.5+3*DBLE(TC)/100)       
    !    H = CUMFPT2ABS(T)
    !    WRITE(88,*) TC, T, H
    !    PRINT*, TC, T, H
    ! ENDDO
    ! CLOSE(88)
    

    ! sample a FPT, using inversion of cumulative function
    TMP = GRND()
    !PRINT*, 'TESTX2:', TMP
    ! look for a solution with t>0
    FPT=FSOLVE_INCR(CUMFPT2ABS,TMP,0D0,INTERPTYPE=3)

    ! rescale to actual units
    FPT = FPT*L**2/D
  CONTAINS
    DOUBLE PRECISION FUNCTION CUMFPT2ABS(T)
      ! calculate the cumulative first passage time probability
      ! up to time t
      ! for a particle diffusing between two absorbing boundaries (in 1D)
      ! length units nondimensionalized by L, time units by L^2/D
      ! SAVECOEFF: optional argument, whether to reset saved coefficients (resets by default)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: T
      INTEGER :: NLIM, N
      INTEGER ::  K, KLIM
      DOUBLE PRECISION :: SQT
      LOGICAL :: USESINE = .TRUE.

      IF (T.LE.0) THEN
         CUMFPT2ABS = 0D0
         RETURN
      ENDIF
     
      ! decide whether to use method of images or sine series      
      SQT = SQRT(T)     
      ! n index necessary for convergence
      NLIM = INT(COEFFCHECKN/SQT)           
      KLIM = INT((COEFFCHECKK*SQT+X0-1)/2)*3

      ! Use at least 3 k indices
      KLIM = MAX(KLIM, 3)

      ! decide whether to use sine series
      USESINE = NLIM.LT.KLIM*2

      IF (USESINE) THEN
         ! use sine series summation (converges well for large t)

         IF (NLIM.GT.MAXN) THEN
            PRINT*, 'ERROR IN CUMFPT2ABS: NLIM out of bounds', NLIM, MAXN, T
            stop 1
         ELSEIF (NLIM.GE.NSET) THEN
            ! update those coefficients that haven't
            ! yet been calculated during this sampling run
            NLIST(NSET:NLIM) = (/(N, N = NSET,NLIM)/) 
            IF (WHICHLEAVE.EQ.1) THEN ! absorb to right boundary
               COEFF(NSET:NLIM) = (-1)**NLIST(NSET:NLIM)&
                    & *SIN(NLIST(NSET:NLIM)*PI*X0L)/NLIST(NSET:NLIM)/X0L*2/PI
            ELSE ! absorb to left boundary 
               COEFF(NSET:NLIM) = (-1)**NLIST(NSET:NLIM) &
                    & *SIN(NLIST(NSET:NLIM)*PI*(1-X0L))/NLIST(NSET:NLIM)/(1-X0L)*2/PI
            ENDIF           
            NSET=NLIM
         ENDIF
         
         ! probability particle has hit before this time
         CUMFPT2ABS = 1 + SUM(COEFF(1:NLIM)*EXP(-(PI**2*T)*NLIST(1:NLIM)**2))
      ELSE
         ! use image method summation (converges well for small t)

         IF (KLIM.GT.MAXK) THEN
            PRINT*, 'ERROR IN CUMFPT2ABS: KLIM out of bounds', KLIM, MAXK, T
            stop 1
         ENDIF

         KLIST(1:KLIM+1) = (/(2*K+1, K=0,KLIM)/)
        
         IF (WHICHLEAVE.EQ.1) THEN ! absorb to right boundary
            CUMFPT2ABS = 0D0
            DO K = 1,KLIM+1 
               CUMFPT2ABS = CUMFPT2ABS + (DERF((KLIST(K)+X0L)/SQT/2)  - DERF((KLIST(K)-X0L)/SQT/2))
            ENDDO
            CUMFPT2ABS = CUMFPT2ABS*1/X0L
         ELSE ! absorb to left boundary
            CUMFPT2ABS = 0D0
            DO K = 1,KLIM+1
               CUMFPT2ABS = CUMFPT2ABS + (DERF((KLIST(K)+1-X0L)/SQT/2) - DERF((KLIST(K)-1+X0L)/SQT/2))
            ENDDO
            CUMFPT2ABS = CUMFPT2ABS*1/(1-X0L)
         ENDIF         
      END IF

      IF (CUMFPT2ABS.LT.-10*EPSILON(1D0).OR.CUMFPT2ABS.GT.1+10*EPSILON(1D0)) THEN
         PRINT*, 'ERROR IN CUMFPT2ABS: out of bounds', CUMFPT2ABS, &
              & USESINE,NLIM,KLIM, T,X0L, 10*EPSILON(1D0)
         stop 1
      ENDIF

      !PRINT*, 'TESTX1:', T, CUMFPT2ABS
    END FUNCTION CUMFPT2ABS
  END SUBROUTINE SAMPLEFPT2ABS
  
 
  
  DOUBLE PRECISION FUNCTION FSOLVE_INCR(F,FVAL,XL0,XU0,INTERPTYPE)
    !Search for roots of F(X) = FVAL
    ! Assuming F is a monotonically increasing function
    ! optionally, supply one or both end-points (XL,XU) for search interval
    ! otherwise will look over entire real line
    ! INTERPTYPE: type of interpolation to use when splitting an interval
    ! 1 = bisection
    ! 2 = linear
    ! 3 = logarithmic
    
    USE KEYS, ONLY : SOLVETOLF, SOLVETOLX,SOLVEITER
    
    IMPLICIT NONE
    DOUBLE PRECISION :: XL,XU,FL,FU,FM, XM
    DOUBLE PRECISION, INTENT(IN) :: FVAL
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: XL0, XU0
    INTEGER, INTENT(IN), OPTIONAL :: INTERPTYPE
    DOUBLE PRECISION :: TOLF, TOLX
    INTEGER :: MAXCT, NC, USEINTERP
    
    INTERFACE
       DOUBLE PRECISION FUNCTION F(X)
         DOUBLE PRECISION, INTENT(IN) :: X
       END FUNCTION F
    END INTERFACE

    ! Set type of interpolation
    IF (PRESENT(INTERPTYPE)) THEN
       USEINTERP = INTERPTYPE
    ELSE
       USEINTERP = 1
    END IF
    
    ! solve to the specified tolerance (along x axis)
    TOLF = SOLVETOLF
    TOLX = SOLVETOLX
    ! maximum iterations allowed
    MAXCT = SOLVEITER
    
    ! set an initial lower bound where function is negative
    IF (PRESENT(XL0)) THEN
       XL = XL0
       FL = F(XL)-FVAL
       IF (FL.EQ.0) THEN
          FSOLVE_INCR = XL
          RETURN
       ELSEIF (FL.GT.0) THEN
          PRINT*, 'ERROR IN FZERO: function is positive on lower bound'
          STOP 1
       ENDIF
    ELSE
       XL = 0D0 ! guess a lower bound
       NC = 0
       FL = F(XL)-FVAL
       DO WHILE (FL.GT.0)
          XL = -(10D0)**NC
          NC = NC+1
          IF (NC.GT.MAXCT) THEN
              PRINT*, 'ERROR IN FZERO: failed to find lower bound'
              STOP 1
          ENDIF
       END DO
    ENDIF

    ! set an initial upper bound where function is positive
    IF (PRESENT(XU0)) THEN
       XU = XU0
       FU = F(XU)-FVAL
       IF (FU.EQ.0) THEN
          FSOLVE_INCR = XL
          RETURN
       ELSEIF (FU.LT.0) THEN
          PRINT*, 'ERROR IN FZERO: function is negative on upper bound'
          STOP 1
       ENDIF
    ELSE
       XU = XL ! guess an upper bound
       NC = 0
       FU = F(XU)-FVAL
       DO WHILE (FU.LT.0)
          XU = XL + (10D0)**NC
          FU = F(XU)-FVAL
          NC = NC+1
          IF (NC.GT.MAXCT) THEN
              PRINT*, 'ERROR IN FZERO: failed to find upper bound', NC, MAXCT
              STOP 1
          ENDIF
       END DO
    ENDIF


    IF (USEINTERP.EQ.1) THEN
       ! bisect interval
       XM = (XL+XU)/2
    ELSEIF (USEINTERP.EQ.2) THEN
       ! linear interpolation    
       XM = XL + (XU-XL)*(-FL)/(FU-FL)
    ELSEIF(USEINTERP.EQ.3) THEN
       ! logarithmic interpolation
       IF (ABS(XL).LT.EPSILON(1D0)) THEN
          XM = (XL+XU)/2
       ELSE
          XM = XL*(XU/XL)**(-FL/(FU-FL))
       ENDIF
    ELSE
       PRINT*, 'ERROR IN FSOLVE_INCR: invalid USEINTERP', USEINTERP
       stop 1
    ENDIF
    FM = F(XM)-FVAL
    
    NC = 0 ! number of bisection cycles
    
    ! Bisection within interval to find solution
    DO WHILE(ABS(FM)>TOLF.OR.MIN(XU-XM, XM-XL)>TOLX)
       NC = NC+1
       IF (NC.GT.MAXCT) THEN
          PRINT*, 'ERROR IN FSOLVE_INCR: too many iterations', NC, MAXCT
          STOP 1
       ENDIF
       
       IF (FM.GT.0) THEN
          XU = XM
          FU = F(XU)-FVAL
       ELSE
          XL = XM
          FL = F(XL)-FVAL
       ENDIF

       
       IF (USEINTERP.EQ.1) THEN
          ! bisect interval
          XM = (XL+XU)/2
       ELSEIF (USEINTERP.EQ.2) THEN
          ! linear interpolation    
          XM = XL + (XU-XL)*(-FL)/(FU-FL)
       ELSEIF(USEINTERP.EQ.3) THEN
          ! logarithmic interpolation
          IF (ABS(XL).LT.EPSILON(1D0)) THEN
             XM = (XL+XU)/2
          ELSE
             XM = XL*(XU/XL)**(-FL/(FU-FL))
          ENDIF

       ENDIF

       FM = F(XM)-FVAL
      ! print*, 'testx3:', xl, xm, xu
    END DO
    FSOLVE_INCR = XM
  END FUNCTION FSOLVE_INCR

  SUBROUTINE BROWNDYNFPT2ABS(X0,NPART,NSTEP,DELT,WHICHLEAVE,LEAVETIME)
    ! use brownian dynamics simulations to sample
    ! FPT on a 1D interval with 2 absorbing boundaries
    ! Assumes length units normalized by L, time units by L^2/D
    ! input:
    ! X0: initial particle position
    ! NPART: number of particles
    ! NSTEPS: number of steps to run (will stop earlier if all absorb)
    ! DELT: timestep
    ! output:
    ! WHICHLEAVE: which side particle leaves on (0=left, 1=right, -1=never left)
    ! LEAVETIME: time of leaving (up to an error of delt)
    USE KEYS, ONLY : PRINTEVERY
    USE mt19937, ONLY : RNORM
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X0, DELT
    INTEGER, INTENT(IN) :: NPART,NSTEP
    INTEGER, INTENT(OUT) :: WHICHLEAVE(NPART)
    DOUBLE PRECISION, INTENT(OUT) :: LEAVETIME(NPART)
    DOUBLE PRECISION :: POS(NPART)
    INTEGER :: STEP, PC
    DOUBLE PRECISION :: SQT, CURTIME, DELX
    
    
    ! initialize particle positions
    POS = X0
    CURTIME = 0D0
    
    ! which side particle leaves on
    WHICHLEAVE = -1
    ! time of leaving
    LEAVETIME = HUGE(1D0)    
    
    ! run BD steps
    SQT = SQRT(2*DELT)
    DO STEP = 1,NSTEP       
       CURTIME = CURTIME + DELT

       IF (ALL(WHICHLEAVE.GE.0)) THEN
          ! all particles are finished
          PRINT*, 'all particles have been absorbed'
          EXIT
       END IF
       
       DO PC = 1,NPART ! propagate forward un-done particles
          IF (WHICHLEAVE(PC).GE.0) CYCLE ! particle not done yet

          ! sample step from a gaussian
          DELX = RNORM()*SQT         
          POS(PC) = POS(PC)+DELX

          IF (POS(PC).LT.0) THEN
             WHICHLEAVE(PC)=0
             LEAVETIME(PC) = CURTIME
          ELSEIF (POS(PC).GT.1D0) THEN
             WHICHLEAVE(PC) = 1
             LEAVETIME(PC) = CURTIME
          END IF
       ENDDO

       IF (MOD(STEP,PRINTEVERY).EQ.0) THEN
          PRINT*, STEP, CURTIME, COUNT(WHICHLEAVE.GE.0), POS(1)
       ENDIF
    ENDDO
    
  END SUBROUTINE BROWNDYNFPT2ABS
  
END MODULE DIFFUTIL
