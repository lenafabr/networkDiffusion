MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE  
  
  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE, NETFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS
  INTEGER :: SNAPSHOTEVERY

  ! ------------
  ! network geometry and setup
  ! ------------  
  INTEGER :: MAXBRANCH ! max number of branches (per node) that can be input in the network

  ! solving diffusion equations
  DOUBLE PRECISION :: SOLVETOL ! solution tolerance (in function)
  INTEGER :: SOLVEITER ! number of iterations allowed in function solution
  
END MODULE KEYS
