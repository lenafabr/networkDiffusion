MODULE NETWORKPROPAGATE
  ! subroutines for propagating particles within a network structure
  IMPLICIT NONE

CONTAINS
  SUBROUTINE RUNPROPAGATE(NETP,NPART,STARTNODE,MAXSTEPS,EDGEFPT)
    ! Propagate many particles over the network, starting at a specific node
    ! propagation is asynchronous, so different particles
    ! are at different times at each step
    ! --------
    ! Inputs:
    ! --------
    ! NETP: Pointer to network object
    ! (all structural and geometric info preset, edges at each node sorted by increasing length)
    ! NPART: number of particles to run
    ! STARTNODE: node on the network where they all start
    ! MAXSTEPS: maximum number of steps to run
    ! -------
    ! Outputs:
    ! -------
    ! EDGEFPT(P,E): for edge E, the first time particle P visits it    
    ! ----------
    ! WARNING:
    ! ---------
    ! when a particle propagates directly from one node to another
    ! then the FPT for that edge is marked as the arrival time
    ! this is not entirely right. Should probably put checkpoints in the middle
    ! of each edge, or else count FPTs to nodes
    ! SHOULD PROBABLY REWRITE THIS FOR FPTS TO NODES

    ! STOPPED HERE. THIS NEEDS VARIABLE DECLARATIONS, ETC
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NPART, STARTNODE, MAXSTEPS
    DOUBLE PRECISION, INTENT(OUT) :: EDGEFPT(NPART,NETP%NEDGE)
    INTEGER :: PLOC

    IF (.NOT.NETP%ARRAYSET) THEN
       PRINT*, 'ERROR IN STEPFROMNODE: network has not yet been set up!'
       STOP 1
    ELSEIF(NID.GT.NETP%NNODE) THEN
       PRINT*, 'ERROR IN STEPFROMNODE: network does not have enough nodes'
       STOP 1
    ENDIF

    EDGEFPT = -1D0

    ! start all particles on a single node
    ! ONNODE keeps track of which particles are on nodes vs edges
    ONNODE = .TRUE.
    ! PLOC is the node or edge index for each particle    
    PLOC = STARTNODE
    
    DO STEP = 1,MAXSTEPS
       PRINT*, 'STEP, # that have hit last edge:', STEP, COUNT(EDGEFPT(:,NETP%NEDGE)>0)
       DO PC = 1,NPART
          IF (DONE(PC)) CYCLE
          
          IF (ONNODE(PC)) THEN
             CALL STEPFROMNODE(NETP, PLOC(PC), DCOEFF,DELT, LEAVENODE, LEAVEEDGE, LEAVELEN)
             CURTIMES(PC) = CURTIMES(PC)+DELT
             
             IF (EDGEFPT(PC,LEAVEEDGE).LT.0) THEN
                ! mark down first passage time to this edge by this particle
                EDGEFPT(PC,LEAVEEDGE) = CURTIMES(PC)
             ENDIF
          ELSE ! propagate from edge to nearest node
             CALL SAMPLEFPT2ABS(L0,NETP%EDGELEN(PLOC(PC)),DCOEFF,DELT,WHICHLEAVE)
             CURTIMES(PC) = CURTIMES(PC)+DELT
             
             PLOC(PC)= NETP%EDGENODE(PLOC(PC),WHICHLEAVE+1)
             ONNODE(PC) = .TRUE.
          ENDIF

          ! Stop tracking the particle once it has hit every edge
          IF (ALL(EDGEFPT(PC,:).GT.0)) THEN
             DONE(P) = .TRUE.
          END IF
                          
       ENDDO
    ENDDO
  END SUBROUTINE RUNPROPAGATE
  
  ! SUBROUTINE STEPFROMEDGE(NETP,EID,L0,DCOEFF,DELT,LEAVEID)
  !   ! for a particle starting somewhere along an edge,
  !   ! propagate to one of the nearest nodes
  !   ! ---------
  !   ! Inputs
  !   ! --------
  !   ! NETP: Pointer to network object
  !   ! (all structural and geometric info preset, edges at each node sorted by increasing length)
  !   ! eID: edge ID where particle starts
  !   ! L0: initial position along the edge where the particle starts
  !   ! DCOEFF: diffusion coefficient
  !   ! --------
  !   ! Outputs
  !   ! -------
  !   ! DELT: propagation time
  !   ! LEAVEID: Node to which it propagates
  !   USE DIFFUTIL, ONLY : SAMPLEFPT2ABS
  !   USE NETWORKUTIL, ONLY : NETWORK
  !   USE MT19937, ONLY : GRND
  !   IMPLICIT NONE
  !   TYPE(NETWORK), POINTER :: NETP
  !   INTEGER, INTENT(IN) :: EID
  !   DOUBLE PRECISION, INTENT(IN) :: DCOEFF
  !   DOUBLE PRECISION, INTENT(OUT) :: DELT
  !   INTEGER, INTENT(OUT) :: LEAVEID

  !   IF (.NOT.NETP%ARRAYSET) THEN
  !      PRINT*, 'ERROR IN STEPFROMNODE: network has not yet been set up!'
  !      STOP 1
  !   ELSEIF(NID.GT.NETP%NNODE) THEN
  !      PRINT*, 'ERROR IN STEPFROMNODE: network does not have enough nodes'
  !      STOP 1
  !   ENDIF

  !   ! propagate through the interval
  !   CALL SAMPLEFPT2ABS(L0,NETP%EDGELEN(EID),DCOEFF,DELT,WHICHLEAVE)

  !   LEAVEID = NETP%EDGENODE(EID,WHICHLEAVE+1)
    
  ! END SUBROUTINE SETPFROMEDGE
  
  SUBROUTINE STEPFROMNODE(NETP, NID, DCOEFF,DELT, LEAVENODE, LEAVEEDGE,LEAVELEN)
    ! for a particle currently located at a node, get first passage time
    ! and propagate forward
    ! Treat the available branches as two half-intervals of different lengths
    ! If the half-interval length is shorter than the branch length propagate
    ! to the appropriate distance along the branch.
    ! If the half-interval length corresponds to the branch length, propagate to the next node
    ! If there are an odd number of branches, create a phantom branch
    ! of same length as the shortest branch
    ! A particle that exits on that branch is then placed at the appropriate
    ! distance on any one of the other branches
    ! ----
    ! Inputs
    ! -----
    ! NETP: pointer to network object
    ! (all structural and geometric info preset, edges at each node sorted by increasing length)
    ! NID: node ID where particle starts
    ! DCOEFF: diffusion coefficient
    ! ------
    ! Outputs:
    ! ------
    ! DELT: propagation time
    ! LEAVENODE: index of  node that the particle propagates to (if it reaches node)
    ! if particle does not reach node, then LEAVENODE=-1
    ! LEAVEEDGE: if particle propagates to node, then this is the edge traversed
    ! if particle does not reach node, this is the edge it ends up on
    ! LEAVELEN: for particle that propagates into an edge, length along the edge where it ends up
    USE DIFFUTIL, ONLY : SAMPLEFPT2ABS
    USE NETWORKUTIL, ONLY : NETWORK
    USE KEYS, ONLY : NODETOL
    USE MT19937, ONLY : GRND
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NID
    DOUBLE PRECISION, INTENT(IN) :: DCOEFF
    DOUBLE PRECISION, INTENT(OUT) :: DELT
    INTEGER, INTENT(OUT) ::  LEAVENODE, LEAVEEDGE
    DOUBLE PRECISION, INTENT(OUT) :: LEAVELEN

    INTEGER :: NE, NE1, PICKBRANCH, WHICHLEAVE
    DOUBLE PRECISION :: LEN1, LEN2, PICKLEN, X0, LTOT, tmp

    IF (.NOT.NETP%ARRAYSET) THEN
       PRINT*, 'ERROR IN STEPFROMNODE: network has not yet been set up!'
       STOP 1
    ELSEIF(NID.GT.NETP%NNODE) THEN
       PRINT*, 'ERROR IN STEPFROMNODE: network does not have enough nodes'
       STOP 1
    ENDIF

    LEAVELEN = 0D0
    
    ! degree of the node
    NE = NETP%NODEDEG(NID)
    
    ! split branches into 2 groups by length
    LEN1 = NETP%NODELEN(NID,1)    
    LEN2 = NETP%NODELEN(NID,NE/2+1)
    ! Number of branches in each group
    IF (MOD(NE,2).EQ.0) THEN
       NE1 = NE/2
    ELSE
       NE1 = (NE+1)/2
    ENDIF
    
    ! convert to a 1d interval
    X0 = LEN1
    LTOT = LEN1+LEN2

    ! propagate through the interval
    CALL SAMPLEFPT2ABS(X0,LTOT,DCOEFF,DELT,WHICHLEAVE)

    !print*, 'testx1:', X0, LTOT, NE1, DELT, WHICHLEAVE
    
    IF (WHICHLEAVE.EQ.0) THEN
       ! exited on shorter side
       PICKLEN = LEN1
       
       ! decide which branch it left on
       TMP = GRND()
       PICKBRANCH = CEILING(TMP*NE1)

       IF (MOD(NE,2).EQ.1.AND.PICKBRANCH.EQ.NE1) THEN
          ! propagated to phantom branch, distribute over all real branches
          TMP = GRND()
          PICKBRANCH = CEILING(TMP*NE)
       ENDIF

       !PRINT*, 'TESTX2:', PICKBRANCH
    ELSE
       ! exited on longer side
       PICKLEN = LEN2
       
       ! decide which branch it left on
       TMP = GRND()
       IF (MOD(NE,2).EQ.0) THEN          
          PICKBRANCH = CEILING(TMP*NE1)+NE1
       ELSE
          PICKBRANCH = CEILING(TMP*NE1)+NE1-1
       ENDIF       
    ENDIF
    
    ! Check if it has gone all the way to the end of the chosen edge
    IF (NETP%NODELEN(NID,PICKBRANCH) - PICKLEN.LT.NODETOL) THEN
       ! Propagated to next node
       LEAVENODE = NETP%NODENODE(NID,PICKBRANCH)
       LEAVEEDGE = NETP%NODEEDGE(NID,PICKBRANCH)
    ELSE
       ! Propagated to some point along the edge
       LEAVEEDGE = NETP%NODEEDGE(NID,PICKBRANCH)
       LEAVENODE = -1
       IF (NETP%EDGENODE(LEAVEEDGE,2).EQ.NID) THEN
          ! node is at end of edge, flip direction
          LEAVELEN = NETP%EDGELEN(LEAVEEDGE)-PICKLEN
       ELSE
          LEAVELEN = PICKLEN
       ENDIF
    ENDIF
       
    
  END SUBROUTINE STEPFROMNODE
END MODULE NETWORKPROPAGATE
