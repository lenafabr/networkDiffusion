MODULE DIFFUTIL
  ! utilities for sampling diffusive motion of particles

  IMPLICIT NONE
CONTAINS

  ! DOUBLE PRECISION FUNCTION CUMFPT2ABS(T,PARAM,NP)
  !   ! calculate the cumulative first passage time probability
  !   ! up to time t
  !   ! for a particle diffusing between two absorbing boundaries (in 1D)
  !   ! length units nondimensionalized by L, time units by L^2/D
  !   ! Parameters:
  !   ! PARAM(1)= x0 (initial position)
  !   ! PARAM(2) = flag (0,1, or -1)
  !   ! flag of 0 = first passage time to any exit
  !   ! flag of -1 = first passage time given left at left boundary
  !   ! flag of 1 = first passage time given left at right boundary
    
  !   IMPLICIT NONE
  !   DOUBLE PRECISION, INTENT(IN) :: T
  !   INTEGER, INTENT(IN) :: PARAM
  !   DOUBLE PRECISION, INTENT(IN) :: PARAM(NP)
  !   INTEGER, PARAMETER :: MAXN = 500, MAXK = 100
  !   LOGICAL :: USESINE = 1

  !   ! decide whether to use method of images or sine series      
  !   SQT = SQRT(T)
  !   ! n index necessary for convergence
  !   CHECKN = SQRT(-LOG(EPSILON(1D0))/T)/PI
  !   NLIM = INT(CHECKN)
    
  !   IF (USESINE) THEN
  !      ! use sine series summation (converges well for large t)

  !      IF (NLIM.GT.MAXN) THEN
  !         PRINT*, 'ERROR IN CUMFPT2ABS: NLIM out of bounds', NLIM, MAXN, T
  !         stop 1
  !      ELSE                 
  !         NLIST(1:NLIM) = (/(N, N = 1,NLIM)/) 
  !         IF (ROUND(PARAM(2)).EQ.1) THEN ! absorb to right boundary
  !            COEFF(1:NLIM) = (-1)**NLIST(1:NLIM)*SIN(NLIST(1:NLIM)*PI*X0)&
  !                 & /NLIST(1:NLIM)/X0*2/PI
  !         ELSEIF (ROUND(PARAM(2).EQ.-1) THEN ! absorb to left boundary 
  !            COEFF(1:NLIM) = (-1)**NLIST(1:NLIM)*&
  !                 & SIN(NLIST(1:NLIM)*PI*(1-X0))/NLIST(1:NLIM)/(1-X0)*2/PI    
  !         ELSE ! Overall absorbance
  !            print*, 'ERROR IN CUMFPT2ABS: overall absorbance (param(2)=0) not yet implemented'
  !            STOP 1
  !         ENDIF
  !      ENDIF

  !      ! probability particle has hit before this time
  !      CUMFPT2ABS = 1 + SUM(COEFF(1:NLIM)*EXP(-(PI**2*T*GAMMA/L**2)*NLIST(1:NLIM)**2))         
  !   END IF

  !   IF (CUMFPT2ABS.LT.0.OR.CUMFPT2ABS.GT.1) THEN
  !      PRINT*, 'ERROR IN CUMFPT2ABS: out of bounds', CUMFPT2ABS, USESINE,NLIM,T,PARAM(2)
  !      PRINT*, L0,L,LP,T*GAMMA,COEFF(1)
  !      stop 1
  !   ENDIF
    
  ! END FUNCTION CUMFPT2ABS
  
  DOUBLE PRECISION FUNCTION FSOLVE_INCR(F,FVAL,PARAM,XL0,XU0)
    !Search for roots of F(X) = FVAL
    ! Assuming F is a monotonically increasing function
    ! optionally, supply one or both end-points (XL,XU) for search interval
    ! otherwise will look over entire real line
    USE KEYS, ONLY : SOLVETOL, SOLVEITER
    
    IMPLICIT NONE
    DOUBLE PRECISION :: XL,XU,FL,FU,FM, XM
    DOUBLE PRECISION, INTENT(IN) :: FVAL,PARAM(:)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: XL0, XU0
    DOUBLE PRECISION :: TOL
    INTEGER :: MAXCT, NC, NP
    
    INTERFACE
       DOUBLE PRECISION FUNCTION F(X,P,NP)
         INTEGER, INTENT(IN) :: NP
         DOUBLE PRECISION, INTENT(IN) :: X
         DOUBLE PRECISION, INTENT(IN) :: P(NP)
       END FUNCTION F
    END INTERFACE

    NP = SIZE(PARAM)
    
    ! solve to the specified tolerance (along x axis)
    TOL = SOLVETOL
    ! maximum iterations allowed
    MAXCT = SOLVEITER
    
    ! set an initial lower bound where function is negative
    IF (PRESENT(XL0)) THEN
       XL = XL0
       FL = F(XL,PARAM,NP)-FVAL
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
       FL = F(XL,PARAM,NP)-FVAL
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
       FU = F(XU,PARAM,NP)-FVAL
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
       FU = F(XU,PARAM,NP)-FVAL
       DO WHILE (FU.LT.0)
          XU = XL + (10D0)**NC
          FU = F(XU,PARAM,NP)-FVAL
          NC = NC+1
          IF (NC.GT.MAXCT) THEN
              PRINT*, 'ERROR IN FZERO: failed to find upper bound', NC, MAXCT
              STOP 1
          ENDIF
       END DO
    ENDIF

    
    
    XM = XL + (XU-XL)*(-FL)/(FU-FL)    
    FM = F(XM,PARAM,NP)-FVAL
    
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
          FU = F(XU,PARAM,NP)-FVAL
       ELSE
          XL = XM
          FL = F(XL,PARAM,NP)-FVAL
       ENDIF

       XM = XL + (XU-XL)*(-FL)/(FU-FL)   
       FM = F(XM,PARAM,NP)-FVAL              
    END DO
    FSOLVE_INCR = XM
  END FUNCTION FSOLVE_INCR

END MODULE DIFFUTIL
