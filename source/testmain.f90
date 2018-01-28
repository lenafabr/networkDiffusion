PROGRAM MAIN
  ! test various subroutines
  USE KEYS
  USE NETWORKUTIL
  USE GENUTIL
  USE TESTMOD
  USE DIFFUTIL
  IMPLICIT NONE
  TYPE(NETWORK), TARGET :: NET
  TYPE(NETWORK), POINTER :: NETP

  NETP=>NET

  CALL READKEY

  !CALL TESTFZERO
  !CALL TESTNETWORKFROMFILE
  CALL TESTSAMPLEFPT
  
  
CONTAINS
  SUBROUTINE TESTSAMPLEFPT
    ! test routine for sampling first passage time from a 2-abs domain
    IMPLICIT NONE
    DOUBLE PRECISION :: FPT
    INTEGER :: WHICHLEAVE
    INTEGER :: T
    
    ! DO T = 0,100
    !    PRINT*, 1.1**T, EXP(-1.1D0**T)
    ! ENDDO
   CALL SAMPLEFPT2ABS(0.5D0,1D0,1D0,FPT,WHICHLEAVE)
  END SUBROUTINE TESTSAMPLEFPT
  
  SUBROUTINE TESTNETWORKFROMFILE
    ! test ability to read in network from file
    IMPLICIT NONE
    INTEGER :: EI, NI
    
    CALL NETWORKFROMFILE(NETP,NETFILE)

    PRINT*, 'Edge connectivities and lengths:'
    DO EI = 1,NETP%NEDGE
       PRINT*, EI, NETP%EDGENODE(EI,:), NETP%EDGELEN(EI)
    END DO

    PRINT*, 'Edges connected to each node:'
    DO NI = 1,NETP%NNODE
       PRINT*, NI, NETP%NODEEDGE(NI,1:NETP%NODEDEG(NI))
    ENDDO
    
    PRINT*, 'Nodes connected to each node:'
    DO NI = 1,NETP%NNODE
       PRINT*, NI, NETP%NODENODE(NI,1:NETP%NODEDEG(NI))
    ENDDO
    
  END SUBROUTINE TESTNETWORKFROMFILE

 
END PROGRAM MAIN
