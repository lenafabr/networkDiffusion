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
    INTEGER, PARAMETER :: MAXN = 500, MAXK = 100
    ! list of N values for sine series
    INTEGER :: NLIST(MAXN), NSET
    ! coefficients for sine series
    DOUBLE PRECISION :: COEFF(MAXN)
    ! precalculate a coefficient for figuring out how
    ! many N indices are needed for sin series convergence
    DOUBLE PRECISION :: COEFFCHECKN

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

    OPEN(UNIT=88,FILE='test.out')
    DO TC = 1,100
       T = 10**(-2+4*DBLE(TC)/100)
       H = CUMFPT2ABS(T)
       WRITE(88,*) TC, T, H
       PRINT*, TC, T, H
    ENDDO
    CLOSE(88)
    
  CONTAINS
    DOUBLE PRECISION FUNCTION CUMFPT2ABS(T)
      ! calculate the cumulative first passage time probability
      ! up to time t
      ! for a particle diffusing between two absorbing boundaries (in 1D)
      ! length units nondimensionalized by L, time units by L^2/D
      
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: T
      INTEGER :: NLIM, N
      DOUBLE PRECISION :: SQT
      LOGICAL :: USESINE = .TRUE.
      
      
      ! decide whether to use method of images or sine series      
      SQT = SQRT(T)     
      ! n index necessary for convergence
      NLIM = INT(COEFFCHECKN/SQT)     
     
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
               COEFF(NSET:NLIM) = (-1)**NLIST(NSET:NLIM)*&
                    & SIN(NLIST(NSET:NLIM)*PI*(1-X0L))/NLIST(NSET:NLIM)/(1-X0L)*2/PI
            ENDIF           
            NSET=NLIM
         ENDIF
         
         ! probability particle has hit before this time
         CUMFPT2ABS = 1 + SUM(COEFF(1:NLIM)*EXP(-(PI**2*T)*NLIST(1:NLIM)**2))  
      END IF

      IF (CUMFPT2ABS.LT.2*EPSILON(1D0).OR.CUMFPT2ABS.GT.1+2*EPSILON(1D0)) THEN
         PRINT*, 'ERROR IN CUMFPT2ABS: out of bounds', CUMFPT2ABS, USESINE,NLIM,T,X0L
         stop 1
      ENDIF

    END FUNCTION CUMFPT2ABS
  END SUBROUTINE SAMPLEFPT2ABS
  
 
  
  DOUBLE PRECISION FUNCTION FSOLVE_INCR(F,FVAL,XL0,XU0)
    !Search for roots of F(X) = FVAL
    ! Assuming F is a monotonically increasing function
    ! optionally, supply one or both end-points (XL,XU) for search interval
    ! otherwise will look over entire real line
    USE KEYS, ONLY : SOLVETOL, SOLVEITER
    
    IMPLICIT NONE
    DOUBLE PRECISION :: XL,XU,FL,FU,FM, XM
    DOUBLE PRECISION, INTENT(IN) :: FVAL
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: XL0, XU0
    DOUBLE PRECISION :: TOL
    INTEGER :: MAXCT, NC
    
    INTERFACE
       DOUBLE PRECISION FUNCTION F(X)
         DOUBLE PRECISION, INTENT(IN) :: X
       END FUNCTION F
    END INTERFACE

    ! solve to the specified tolerance (along x axis)
    TOL = SOLVETOL
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

    
    
    XM = XL + (XU-XL)*(-FL)/(FU-FL)    
    FM = F(XM)-FVAL
    
    NC = 0 ! number of bisection cycles
    
    ! Bisection within interval to find solution
    DO WHILE(ABS(FM)>TOL)
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

       XM = XL + (XU-XL)*(-FL)/(FU-FL)   
       FM = F(XM)-FVAL              
    END DO
    FSOLVE_INCR = XM
  END FUNCTION FSOLVE_INCR

END MODULE DIFFUTIL
