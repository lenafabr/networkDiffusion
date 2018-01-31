MODULE NETWORKUTIL
  ! Utilities for dealing with connectivity and geometry of a network
  ! including defining, manipulating, input/output

  IMPLICIT NONE
  
  TYPE NETWORK
     ! define a network structure (connectivity and geometry)

     ! dimension of the space the network is embedded in
     INTEGER :: DIM 
     
     ! ----------------------
     ! information on network nodes
     ! ----------------------
     INTEGER :: NNODE ! number of nodes
     ! list of other node indices each node connects to
     INTEGER, POINTER :: NODENODE(:,:)
     ! degree (number of branches) of each node
     INTEGER, POINTER :: NODEDEG(:)
     ! list of branch indices each node connects to
     INTEGER, POINTER :: NODEEDGE(:,:)          
     ! lengths of all edges leaving the node
     DOUBLE PRECISION, POINTER :: NODELEN(:,:)
     ! spatial location of node
     DOUBLE PRECISION, POINTER :: NODEPOS(:,:)
     
     
     ! ------------------
     ! information on network branches
     ! ------------------
     INTEGER :: NEDGE ! Number of branches
     ! nodes at the start and end of each branch
     INTEGER, POINTER :: EDGENODE(:,:)
     ! spatial position of branch starting points
     ! branch direction and length
     DOUBLE PRECISION, POINTER :: EDGESTART(:,:), EDGEDIR(:,:), EDGELEN(:)

     
     ! Have arrays been set up already
     LOGICAL :: ARRAYSET = .FALSE.
  END type NETWORK

CONTAINS
  
  SUBROUTINE NETWORKFROMFILE(NETP,NETFILE)
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI
    USE KEYS, ONLY : MAXBRANCH
    USE GENUTIL, ONLY : NORMALIZE
    ! Set up (allocate) a network structure, reading in connectivity and
    ! geometry from an input file    
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    CHARACTER(LEN=*), INTENT(IN) :: NETFILE
    LOGICAL :: LDUM, CASESET
    LOGICAL :: FILEEND=.FALSE., DIMSET = .FALSE.
    CHARACTER*100 :: WORD
    INTEGER :: NITEMS, NNODE, NEDGE, NODE1, NODE2, DIM, EID, NID, NE
    DOUBLE PRECISION :: LEN, LENS(MAXBRANCH), LENSAVE(MAXBRANCH)
    INTEGER :: TMPARRAY(MAXBRANCH)
    
    INTEGER, PARAMETER :: NF = 55 ! input file unit number
   
    ! deallocate any previously set arrays
    IF (NETP%ARRAYSET) CALL CLEANUPNETWORK(NETP)

    ! go through file and count nodes and branches (total and from each node)
    PRINT*, 'Reading network structure file: ', NETFILE
    INQUIRE(FILE=NETFILE,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
       PRINT*, 'ERROR in NETWORKFROMFILE: network file ', TRIM(ADJUSTL(NETFILE)), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')

     ! Count number of nodes and edges, set spatial dimension
     DIMSET = .FALSE.
     NNODE = 0; NEDGE = 0
     DO 
        CALL READLINE(NF,FILEEND,NITEMS)
        
        IF (FILEEND.and.nitems.eq.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        IF (WORD.EQ.'NODE') THEN
           NNODE = NNODE+1
           ! set spatial dimension of network
           IF (DIMSET.AND.DIM.NE.NITEMS-2) THEN
              PRINT*, 'ERROR IN NETWORKFROMFILE: dimension of node positions is inconsistent', DIM, NITEMS
              STOP 1
           ENDIF           
           DIM = NITEMS-2           
        ELSEIF (WORD.EQ.'EDGE') THEN
           NEDGE = NEDGE + 1         
        ENDIF
     ENDDO
     CLOSE(NF)

     PRINT*, 'Number of nodes and edges found: ', NNODE, NEDGE
     
     ! allocate arrays
     CALL SETUPNETWORK(NETP,NNODE,NEDGE,DIM,MAXBRANCH)


     PRINT*, 'Reading in node positions and edge connectivity...'
     OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')
     ! Get node and edge information
     DO 
        CALL READLINE(NF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT
        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE
        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)
        ! Skip any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        IF (WORD.EQ.'NODE') THEN
           CALL READI(NID)
           IF (NID.LT.1.OR.NID.GT.NNODE) THEN
              PRINT*, 'ERROR IN NETWORKFROMFILE: invalid  node index', NNODE, NID
              STOP 1
           ENDIF
           CALL READF(NETP%NODEPOS(NID,1))
           CALL READF(NETP%NODEPOS(NID,2))           
        ELSEIF (WORD.EQ.'EDGE') THEN                              
           CALL READI(EID) ! edge id
           CALL READI(NODE1)
           CALL READI(NODE2)

           IF (NODE1.LT.1.OR.NODE2.LT.1.OR.EID.LT.1&
                & .OR.NODE1.GT.NNODE.OR.NODE2.GT.NNODE.OR.EID.GT.NEDGE) THEN
              PRINT*, 'ERROR IN NETWORK FROM FILE: &
                   & bad edge or node indices while reading edge.', &
                   ' EID, NODE1, NODE2, NNODE, NEDGE:', &
                   & EID, NODE1, NODE2, NNODE, NEDGE
              STOP 1
           ENDIF

           ! nodes connected to this edge           
           NETP%EDGENODE(EID,:) = (/NODE1,NODE2/)                
        ENDIF        
     END DO

     PRINT*, 'Setting up connectivity and edge lengths ...'

     DO EID = 1,NEDGE
        NODE1 = NETP%EDGENODE(EID,1)
        NODE2 = NETP%EDGENODE(EID,2)

        ! edge start, length, and (normalized) direction
        NETP%EDGESTART(EID,:) = NETP%NODEPOS(NODE1,:)
        NETP%EDGEDIR(EID,:) = NETP%NODEPOS(NODE2,:)-NETP%NODEPOS(NODE1,:)   
        CALL NORMALIZE(NETP%EDGEDIR(EID,:), LEN)          
        NETP%EDGELEN(EID) = LEN       

        ! increment degrees of the nodes
        NETP%NODEDEG(NODE1) = NETP%NODEDEG(NODE1)+1
        NETP%NODEDEG(NODE2) = NETP%NODEDEG(NODE2)+1

        IF ( MAX(NETP%NODEDEG(NODE1),NETP%NODEDEG(NODE2)).GT.MAXBRANCH) THEN
           PRINT*, 'ERROR IN NETWORKFROMFILE: node degree exceeds maximum.',&
                & NODE1, NODE2, MAXBRANCH, NETP%NODEDEG(NODE1), NETP%NODEDEG(NODE2)
           STOP 1
        ENDIF

        ! edges connected to each node
        NETP%NODEEDGE(NODE1,NETP%NODEDEG(NODE1)) = EID
        NETP%NODEEDGE(NODE2,NETP%NODEDEG(NODE2)) = EID

        ! nodes connected to each node
        NETP%NODENODE(NODE1,NETP%NODEDEG(NODE1)) = NODE2
        NETP%NODENODE(NODE2,NETP%NODEDEG(NODE2)) = NODE1              
     END DO

     ! for each node, sort edges by length
     DO NID = 1,NETP%NNODE
        NE = NETP%NODEDEG(NID)
        LENS(1:NE) = NETP%EDGELEN(NETP%NODEEDGE(NID,1:NE))
        LENSAVE = LENS
        TMPARRAY = NETP%NODEEDGE(NID,:)        
        CALL SORT2INT(NE,LENS,TMPARRAY)
        NETP%NODEEDGE(NID,1:NE) = TMPARRAY(1:NE)
        
        TMPARRAY = NETP%NODENODE(NID,:)
        CALL SORT2INT(NE,LENSAVE,TMPARRAY)
        NETP%NODENODE(NID,1:NE) = TMPARRAY(1:NE)

        ! length of all edges leaving the node
        NETP%NODELEN(NID,1:NE)= LENS(1:NE)
     ENDDO
     
  END SUBROUTINE NETWORKFROMFILE
  
  SUBROUTINE SETUPNETWORK(NETP,NNODE,NEDGE,DIM,MAXBRANCH)
    ! set up a network by allocating arrays
    ! NETP: Pointer to a network
    ! NNODE: number of nodes
    ! NEDGE: number of branches
    ! DIM: spatial dimensionality where network resides
    ! MAXBRANCH: maximum branches attached to each node
    
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP
    INTEGER, INTENT(IN) :: NNODE, NEDGE, DIM, MAXBRANCH

    NETP%NNODE = NNODE
    NETP%NEDGE = NEDGE
    NETP%DIM = DIM

    ! allocate node data
    ALLOCATE(NETP%NODENODE(NNODE,MAXBRANCH), NETP%NODEEDGE(NNODE,MAXBRANCH), &
         &  NETP%NODEPOS(NNODE,DIM), NETP%NODEDEG(NNODE),NETP%NODELEN(NNODE,MAXBRANCH))

    ! allocate branch data
    ALLOCATE(NETP%EDGENODE(NEDGE,2), NETP%EDGESTART(NEDGE,DIM),&
         & NETP%EDGEDIR(NEDGE,DIM), NETP%EDGELEN(NEDGE))

    NETP%ARRAYSET = .TRUE.
    NETP%NODELEN = HUGE(1D0)
    NETP%NODEDEG = 0
    NETP%NODEEDGE = 0; NETP%EDGENODE = 0; NETP%NODENODE = 0
  END SUBROUTINE SETUPNETWORK

  SUBROUTINE CLEANUPNETWORK(NETP)
    ! deallocate arrays for the network structure
    IMPLICIT NONE
    TYPE(NETWORK), POINTER :: NETP

    IF (NETP%ARRAYSET) THEN
    DEALLOCATE(NETP%NODENODE, NETP%NODEEDGE, NETP%NODEPOS, NETP%NODELEN)
    DEALLOCATE(NETP%EDGENODE, NETP%EDGESTART, NETP%EDGEDIR, NETP%EDGELEN)
 ENDIF
 
    NETP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPNETWORK
END MODULE NETWORKUTIL
