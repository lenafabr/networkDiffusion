MODULE DIFFUTIL
  ! utilities for sampling diffusive motion of particles

  IMPLICIT NONE
CONTAINS
  DOUBLE PRECISION FUNCTION FSOLVE_INCR(F,FVAL,PARAM,XL0,XU0)
    !Search for roots of F(X) = FVAL
    ! Assuming F is a monotonically increasing function
    ! optionally, supply one or both end-points (XL,XU) for search interval
    ! otherwise will look over entire real line
    USE KEYS, ONLY : SOLVETOL, SOLVEITER
    
    IMPLICIT NONE
    INTEGER :: NC
    DOUBLE PRECISION :: XL,XU,FL,FU,FM, XM
    DOUBLE PRECISION, INTENT(IN) :: FVAL,PARAM(2)
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
       DO WHILE (F(XL,PARAM,NP).GT.FVAL)
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
       DO WHILE (F(XU,PARAM,NP).LT.FVAL)
          XU = XL + (10D0)**NC
          NC = NC+1
          IF (NC.GT.MAXCT) THEN
              PRINT*, 'ERROR IN FZERO: failed to find upper bound', CT, MAXCT
              STOP 1
          ENDIF
       END DO
    ENDIF

    
     XM = (XL+XU)/2
     FSOLVE_INCR = XM
    
     NC = 0 ! number of bisection cycles
    
    !LOOP FOR BISECTION METHOD
    DO WHILE(ABS(F(XM,PARAM,NP)-FVAL)>TOL)
       FL = F(XL)
       FU = F(XU)
       FM = F(XM)

       !CHECK IF ROOT IS AT INTERVAL ENDPOINTS
       IF(ABS(FL)<TOL) THEN
          FZERO = XL
          EXIT
       ELSE IF(ABS(FU)<TOL) THEN
          FZERO = XU
          EXIT
       END IF

    !    IF(FL*FU>0.AND.NC.LT.10) THEN
    !       XL = XL/10
    !       XU = XU*10
    !       NC = NC+1
    !    ELSEIF(NC.GT.10) THEN			
    !       PRINT*, "FUNCTION DOES NOT CHANGE SIGN AT ENDPOINTS. CHECK FUNCTION"
    !       PRINT*, FL,FU,XL,XU
    !       STOP 1
    !    END IF

    !    !FIND ROOT CONTAINING INTERVAL
    !    IF (FL*FM<=0) THEN
    !       XU = XM
    !       XM = (XU+XL)/2
    !    ELSE IF(FU*FM<=0) THEN
    !       XL = XM
    !       XM = (XU+XL)/2
    !    END IF

    !    FZERO = XM	

    ! END DO
    
  END FUNCTION FSOLVE_INCR

END MODULE DIFFUTIL
