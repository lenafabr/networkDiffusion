PROGRAM MAIN
  ! test various subroutines
  USE KEYS, ONLY : ACTION

  ! read in parameter file
  
  CALL READKEY

  SELECT CASE(ACTION)
  CASE('GETFPT')
     ! calculate first passage time for all particles to each node
     CALL GETFPTDRIVER
  CASE DEFAULT
     PRINT*, 'Not a valid action: ', ACTION
  END SELECT
  
CONTAINS
  SUBROUTINE GETFPTDRIVER
    ! calculate first passage time for all particles to each node
    USE KEYS, ONLY : NETFILE, NPART, STARTNODE, NSTEP, DCOEFF, OUTFILE
    USE NETWORKUTIL, ONLY : NETWORKFROMFILE, NETWORK, CLEANUPNETWORK
    USE NETWORKPROPAGATE, ONLY : RUNPROPAGATE
    IMPLICIT NONE
    TYPE(NETWORK), TARGET :: NET
    TYPE(NETWORK), POINTER :: NETP
    DOUBLE PRECISION, ALLOCATABLE :: NODEFPT(:,:)
    INTEGER :: PC
    
    NETP=>NET

    CALL NETWORKFROMFILE(NETP,NETFILE)
    ALLOCATE(NODEFPT(NPART,NETP%NNODE))
    
    CALL RUNPROPAGATE(NETP,NPART,STARTNODE,NSTEP,DCOEFF,NODEFPT)

    ! for each particle output FPTs to file
    OPEN(UNIT=99,FILE=OUTFILE)
    DO PC = 1,NPART
       WRITE(99,*) PC, NODEFPT(PC,:)
    ENDDO
    CLOSE(99)
    
    DEALLOCATE(NODEFPT)    

    CALL CLEANUPNETWORK(NETP)
  END SUBROUTINE GETFPTDRIVER
END PROGRAM MAIN