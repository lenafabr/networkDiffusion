MODULE NETWORKPROPAGATE
  ! subroutines for propagating particles within a network structure
  IMPLICIT NONE

CONTAINS
  SUBROUTINE RUNPROPAGATE(NETP,NPART,STARTNODE,MAXSTEPS,DCOEFF,NODEFPT,FINISHTIME)
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
    ! DCOEFF: diffusion coefficient
    ! -------
    ! Outputs:
    ! -------
    ! NODEFPT(P,N): for particle P, the first time it visits node N
    USE NETWORKUTIL, ONLY : NETWORK
    USE DIFFUTIL, ONLY : SAMPLEFPT2ABS
    USE KEYS, ONLY : PRINTEVERY, NODETOL
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NPART, STARTNODE, MAXSTEPS
    DOUBLE PRECISION, INTENT(IN) :: DCOEFF
    DOUBLE PRECISION, INTENT(OUT) :: NODEFPT(NPART,NETP%NNODE)
    ! index of node or edge where particle is located
    INTEGER :: PLOC(NPART), OLDNODE, OLDEDGE
    ! length along edge where particle is located
    DOUBLE PRECISION :: EDGEPOS(NPART), NEWPOS
    INTEGER :: STEP, PC, IND
    LOGICAL :: DONE(NPART), ONNODE(NPART)
    DOUBLE PRECISION :: CURTIMES(NPART)
    DOUBLE PRECISION :: DELT, LEAVELEN, X0, LENGTHTOT
    INTEGER :: LEAVENODE, LEAVEEDGE, WHICHLEAVE
    !DOUBLE PRECISION, POINTER :: EDGECOVER(NPART,NEDGE)
    INTEGER, DIMENSION(NPART,NETP%NEDGE) :: EDGECOVER
    DOUBLE PRECISION, DIMENSION(NPART,NETP%NEDGE) :: EDGELEFT, EDGERIGHT
    DOUBLE PRECISION, DIMENSION(NPART), INTENT(OUT) :: FINISHTIME
    

    IF (.NOT.NETP%ARRAYSET) THEN
       PRINT*, 'ERROR IN STEPFROMNODE: network has not yet been set up!'
       STOP 1
    ELSEIF(STARTNODE.GT.NETP%NNODE.OR.STARTNODE.LE.0) THEN
       PRINT*, 'ERROR IN RUNPROPAGATE: invalid startnode', STARTNODE
       STOP 1
    ENDIF

    ! none of the nodes visited yet    
    NODEFPT = -1D0    
    EDGEPOS = -1D0

    !set all edges to not have been covered yet
    EDGECOVER(:,:) = 0
    !set furthest progress from both directions along each edge to be zero
    EDGELEFT(:,:) = 0;
    DO PC=1,NPART
      EDGERIGHT(PC,:) = NETP%EDGELEN(:)
    END DO
    
    ! start all particles on a single node
    ! ONNODE keeps track of which particles are on nodes vs edges
    ONNODE = .TRUE.
    DONE = .FALSE.
    ! PLOC is the node or edge index for each particle    
    PLOC = STARTNODE
    NODEFPT(:,STARTNODE) = 0D0
    
    CURTIMES = 0D0
    
    DO STEP = 1,MAXSTEPS
       IF (MOD(STEP,PRINTEVERY).EQ.0) THEN
          !PRINT*, 'STEP, # that have hit all nodes, max time:', STEP, COUNT(DONE), CURTIMES(1), PLOC(1), ONNODE(1)
          !, MAXVAL(CURTIMES)
       ENDIF

       IF (ALL(DONE)) THEN
         FINISHTIME = CURTIMES
         EXIT
       ENDIF
       
       DO PC = 1,NPART
          IF (DONE(PC)) CYCLE
          
          IF (ONNODE(PC)) THEN
             CALL STEPFROMNODE(NETP, PLOC(PC), DCOEFF,DELT, LEAVENODE, LEAVEEDGE, LEAVELEN,&
                  & EDGECOVER, EDGELEFT, EDGERIGHT, NPART, PC)
             CURTIMES(PC) = CURTIMES(PC)+DELT

             IF (LEAVENODE.GT.0) THEN            
                IF (NODEFPT(PC,LEAVENODE).LT.-0.5) THEN
                   ! mark down first visit to this node by this particle
                   NODEFPT(PC,LEAVENODE) = CURTIMES(PC)
                ENDIF
                PLOC(PC) = LEAVENODE
                EDGECOVER(PC,LEAVEEDGE) = 1 !marking edge as completely visited
             ELSE
                OLDNODE = PLOC(PC)

                ONNODE(pc) = .FALSE. ! particle is along an edge now
                EDGEPOS(PC) = LEAVELEN
                PLOC(PC) = LEAVEEDGE

                !figuring out partial marking of edge visit
                IF(EDGECOVER(PC,LEAVEEDGE).EQ.0) THEN
                  IF (NETP%EDGENODE(LEAVEEDGE,1).EQ.OLDNODE) THEN
                    IF (LEAVELEN.GT.EDGELEFT(PC,LEAVEEDGE)) THEN
                      EDGELEFT(PC,LEAVEEDGE) = LEAVELEN
                    ENDIF
                  ELSE
                    IF(LEAVELEN.LT.EDGERIGHT(PC,LEAVEEDGE)) THEN
                      EDGERIGHT(PC,LEAVEEDGE) = LEAVELEN
                    ENDIF
                  ENDIF

                  IF (EDGERIGHT(PC,LEAVEEDGE) - EDGELEFT(PC,LEAVEEDGE).LT.NODETOL) THEN
                    EDGECOVER(PC,LEAVEEDGE) = 1
                  ENDIF
                ENDIF

             ENDIF
          ELSE IF ((EDGECOVER(PC,PLOC(PC)).EQ.1).OR.& ! propagate from edge to a neighboring node or to visit coverage on edge
                  &((EDGELEFT(PC,PLOC(PC)).LT.NODETOL).AND.&
                  &(NETP%EDGELEN(PLOC(PC)) - EDGERIGHT(PC,PLOC(PC)).LT.NODETOL))) THEN
             CALL SAMPLEFPT2ABS(EDGEPOS(PC),NETP%EDGELEN(PLOC(PC)),DCOEFF,DELT,WHICHLEAVE)
             CURTIMES(PC) = CURTIMES(PC)+DELT

             OLDEDGE = PLOC(PC)
             
             PLOC(PC)= NETP%EDGENODE(PLOC(PC),WHICHLEAVE+1)
             ONNODE(PC) = .TRUE.

             IF (NODEFPT(PC,PLOC(PC)).LT.-0.5) THEN
                ! Mark down first visit to this node by this particle
                NODEFPT(PC,PLOC(PC)) = CURTIMES(PC)
             ENDIF

             ! figuring out partial marking of edge visit
             IF(EDGECOVER(PC,OLDEDGE).EQ.0) THEN
               IF (NETP%EDGENODE(OLDEDGE,1).EQ.PLOC(PC)) THEN
                 IF(EDGEPOS(PC).GT.EDGELEFT(PC,OLDEDGE)) THEN
                   EDGELEFT(PC,OLDEDGE) = EDGEPOS(PC)
                 ENDIF
               ELSE
                 IF(EDGEPOS(PC).LT.EDGERIGHT(PC,OLDEDGE)) THEN
                   EDGERIGHT(PC,OLDEDGE) = EDGEPOS(PC)
                 ENDIF
               ENDIF

               IF (EDGERIGHT(PC,OLDEDGE) - EDGELEFT(PC,OLDEDGE).LT.NODETOL) THEN
                 EDGECOVER(PC,OLDEDGE) = 1
               ENDIF
             ENDIF
          ELSE !propagate from edge to neighboring node or to limit of visited part of edge
            IF(EDGELEFT(PC,PLOC(PC)) - EDGEPOS(PC).GT.-NODETOL) THEN
              X0 = EDGEPOS(PC)
              LENGTHTOT = EDGERIGHT(PC,PLOC(PC))
            ELSE
              X0 = EDGEPOS(PC) - EDGELEFT(PC,PLOC(PC))
              LENGTHTOT = NETP%EDGELEN(PLOC(PC)) - EDGELEFT(PC,PLOC(PC))
            ENDIF

            CALL SAMPLEFPT2ABS(X0,LENGTHTOT,DCOEFF,DELT,WHICHLEAVE)
            CURTIMES(PC) = CURTIMES(PC)+DELT

            IF(WHICHLEAVE.EQ.0) THEN
              NEWPOS = EDGEPOS(PC) - X0
            ELSE
              NEWPOS = EDGEPOS(PC) + LENGTHTOT - X0
            ENDIF

            IF((ABS(NEWPOS).LT.NODETOL).OR.(ABS(NEWPOS-NETP%EDGELEN(PLOC(PC))).LT.NODETOL)) THEN
              !its on a node
              OLDEDGE = PLOC(PC)
              PLOC(PC)= NETP%EDGENODE(PLOC(PC),WHICHLEAVE+1)
              ONNODE(PC) = .TRUE.
              IF (NODEFPT(PC,PLOC(PC)).LT.-0.5) THEN !make sure this is done correctly, its just pasted from somewhere else...
                ! mark down first visit to this node by this particle
                NODEFPT(PC,PLOC(PC)) = CURTIMES(PC)
              ENDIF

              IF((WHICHLEAVE.EQ.0).AND.(EDGELEFT(PC,OLDEDGE) - EDGEPOS(PC).GT.-NODETOL)) THEN
                !returned to same node (start node)
              ELSE IF ((WHICHLEAVE.EQ.1).AND.(EDGEPOS(PC) - EDGERIGHT(PC,OLDEDGE).GT.-NODETOL)) THEN
                !also returned to same node (end node)
              ELSE
                !went to other node, but edgeleft or edgeright was on this node
                EDGECOVER(PC,OLDEDGE) = 1
              ENDIF
            ELSE
              !its still on the edge, has covered entire edge
              EDGECOVER(PC,PLOC(PC)) = 1
              EDGEPOS(PC) = NEWPOS
            ENDIF
          ENDIF

          ! Stop tracking the particle once it has hit every node
          !IF (ALL(NODEFPT(PC,:).GT.-0.5)) THEN
           !  DONE(PC) = .TRUE.
          !END IF

          ! Stop tracking the particle once it has hit every edge
          IF (ALL(EDGECOVER(PC,:).EQ.1)) THEN
            DONE(PC) = .TRUE.
          ENDIF
                  
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
  
  SUBROUTINE STEPFROMNODE(NETP, NID, DCOEFF,DELT, LEAVENODE, LEAVEEDGE,LEAVELEN,&
             EDGECOVER, EDGELEFT, EDGERIGHT, NPART, PC)
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
    USE KEYS, ONLY : NODETOL, MAXBRANCH
    USE MT19937, ONLY : GRND
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NID
    DOUBLE PRECISION, INTENT(IN) :: DCOEFF
    DOUBLE PRECISION, INTENT(OUT) :: DELT
    INTEGER, INTENT(OUT) ::  LEAVENODE, LEAVEEDGE
    DOUBLE PRECISION, INTENT(OUT) :: LEAVELEN
    INTEGER, DIMENSION(NPART,NETP%NEDGE), INTENT(IN) :: EDGECOVER
    DOUBLE PRECISION, DIMENSION(NPART,NETP%NEDGE), INTENT(IN) :: EDGELEFT, EDGERIGHT
    INTEGER, INTENT(IN) :: NPART, PC

    INTEGER :: NE, NE1, PICKBRANCH, WHICHLEAVE, EI
    DOUBLE PRECISION :: LEN1, LEN2, PICKLEN, X0, LTOT, tmp
    DOUBLE PRECISION, DIMENSION(MAXBRANCH) :: NODELENVISIT
    INTEGER :: TMPARRAY(MAXBRANCH)

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

    !find lengths along edges to visited distance
    DO EI=1,NE
      IF(EDGECOVER(PC,NETP%NODEEDGE(NID,EI)).EQ.1) THEN
        NODELENVISIT(EI) = NETP%NODELEN(NID,EI)
      ELSE
        IF (NETP%EDGENODE(NETP%NODEEDGE(NID,EI),1).EQ.NID) THEN
          NODELENVISIT(EI) = EDGERIGHT(PC,NETP%NODEEDGE(NID,EI))
        ELSE
          NODELENVISIT(EI) = NETP%NODELEN(NID,EI) - EDGELEFT(PC,NETP%NODEEDGE(NID,EI))
        ENDIF
      ENDIF
    END DO

    !put lengths in order
    DO EI=1,NE
      TMPARRAY(EI) = EI
!NETP%NODEEDGE(NID,EI)
    END DO
    CALL SORT2INT(NE,NODELENVISIT,TMPARRAY)
    
    ! split branches into 2 groups by length
    !LEN1 = NETP%NODELEN(NID,1)    
    !LEN2 = NETP%NODELEN(NID,NE/2+1)
    LEN1 = NODELENVISIT(1)    
    LEN2 = NODELENVISIT(NE/2+1)
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

    PICKBRANCH = TMPARRAY(PICKBRANCH)
    
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
