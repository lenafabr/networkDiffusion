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
  INTEGER :: PRINTEVERY
  
  ! ------------
  ! network geometry and setup
  ! ------------  
  INTEGER :: MAXBRANCH ! max number of branches (per node) that can be input in the network

  ! -----------------
  ! solving diffusion equations and propagating particles
  ! ----------------
  DOUBLE PRECISION :: SOLVETOLF, SOLVETOLX ! solution tolerance (in function and X value)
  INTEGER :: SOLVEITER ! number of iterations allowed in function solution
  INTEGER :: NPART ! number of particles to propagate
  INTEGER :: NSTEP ! number of steps to propagate for
  DOUBLE PRECISION :: DELT ! time-step for BD sims

  ! How close to a node does the particle have to come to be considered at the node
  ! for propagation through tubules
  DOUBLE PRECISION :: NODETOL 
END MODULE KEYS
