SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! input/output  
  OUTFILE = '*.out'
  DUMPSNAPSHOTS = .false. ! periodically dump chain snapshots
  SNAPSHOTEVERY = 1 ! how often to dump snapshots
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing

  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO 
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)
        CASE('OUTFILE')
           CALL READA(OUTFILE)
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('SNAPSHOTFILE')
           CALL READA (SNAPSHOTFILE)
        CASE('SNAPSHOTS')
           DUMPSNAPSHOTS = .TRUE.
           IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
           IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
           IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('VERBOSE')
           CALL READO(VERBOSE)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO


  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  
  
  

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument 
     ! and additionally the millisecond time 
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)  
  IF (DUMPSNAPSHOTS) THEN
     PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
  ENDIF    
  print*, '----------------------------------------------------'


END SUBROUTINE READKEY
