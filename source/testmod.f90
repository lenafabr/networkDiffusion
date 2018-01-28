MODULE TESTMOD
  ! testing subroutines that need to be in a module
  USE DIFFUTIL, only : fsolve_incr

CONTAINS
   SUBROUTINE TESTFZERO
    ! test function solver
    USE DIFFUTIL
    IMPLICIT NONE
    DOUBLE PRECISION :: XSOLVE, coeff(2)

    COEFF = (/2D0,-2D0/)
    
    XSOLVE = Fsolve_INCR(TESTFUNC,10D0,0d0)
    print*, 'xsolve:', xsolve

  CONTAINS    

  DOUBLE PRECISION FUNCTION TESTFUNC(X)
    ! test function for testing numerical solvers
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X
    
    TESTFUNC = coeff(1)*X**2 + coeff(2)

    PRINT*, 'TESTX2:', X, TESTFUNC
  END FUNCTION TESTFUNC
  END SUBROUTINE TESTFZERO
END MODULE TESTMOD
