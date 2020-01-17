C*******************************************************************************
C Data and functions for TG packets
C*******************************************************************************

      MODULE PACKETS

      IMPLICIT NONE
      SAVE

C-----------------------------------
C Data types
C-----------------------------------

C Basic data structure for Lagrangian "packet" (aka particle, parcel, pulse ...)
C of air; packets are accessed as members of a PacketsInCell list
      TYPE :: Packet
         INTEGER :: ID_GD       ! unique packet identification number
         REAL(8) :: AGE         ! time since created (sec)
         REAL(8) :: X, Y, Z     ! coordinates of current location (m)
         REAL, POINTER :: Mix(:) ! dimension consistent with CGRID, i.e., NSPCSD
      END TYPE Packet

C A member of a doubly linked list of packets
      TYPE :: PacketListMember
         TYPE(Packet), POINTER :: P
         TYPE(PacketListMember), POINTER :: Prev, Next
      END TYPE PacketListMember

C Doubly linked list of packets
      TYPE :: PacketList
         INTEGER :: Count       ! number of packets in list
         TYPE(PacketListMember), POINTER :: Head, Tail
      END TYPE PacketList

C Variable resolution multipliers (horizontal only)
      TYPE :: VarResMult
         INTEGER :: X, Y
      END TYPE VarResMult

C Variable resolution DX and DY
      TYPE :: VarResDelta
         REAL(8) :: X, Y
      END TYPE VarResDelta

C-----------------------------------
C Data
C-----------------------------------

C Eulerian grid data
C
C GOTCHA: (XORIG_GD, YORIG_GD) in HGRD_DEFN.F are coordinates for the CROSS
C point (C=1, R=1), i.e., the CENTER of the cell (C=1, R=1), NOT for its
C lower, left corner, as Figure 12-5 on p. 12-13 in the document "Science
C Algorithms of the EPA Models-3 Community Multiscale Air Quality (CMAQ)
C Modeling System" indicates.
      REAL(8) :: XORIG_DOT, YORIG_DOT ! lower, left corner of cell (C=1, R=1)
      REAL(8) :: DX_GD, DY_GD   ! size of horizontal cells (m)

C Last packet id number assigned (incremented each time a packet is created)
      INTEGER :: LastID_GD

C For each cell of the Eulerian grid, the list of packets currently in the cell
      TYPE(PacketList), POINTER :: PacketsInCell(:,:,:)

C High resolution box:
C - cells are initialized with more packets
C - new packet introduction is more aggressive
      INTEGER :: HR_TOP = 2     ! top layer for high resolution box
      INTEGER :: HR_BEGCOL      ! beginning column for high resolution box
      INTEGER :: HR_ENDCOL      ! ending column for high resolution box
      INTEGER :: HR_BEGROW      ! beginning row for high resolution box
      INTEGER :: HR_ENDROW      ! ending row for high resolution box

C For each cell of the Eulerian grid, its horizontal resolution multipiers and
C deltas
      TYPE(VarResMult), POINTER :: ResMult(:,:,:)
      TYPE(VarResDelta), POINTER :: ResDel(:,:,:)

C Fill method
      CHARACTER(16) :: FILL_METHOD
      CHARACTER(16), PARAMETER :: NO_FILL = 'NO_FILL'
      CHARACTER(16), PARAMETER :: FILL_ALL = 'FILL_ALL'
      CHARACTER(16), PARAMETER :: SPARSE_FILL = 'SPARSE_FILL'

C Methods for setting the composition of a packet created by a fill method
      CHARACTER(16), PARAMETER :: CLOSEST_MIX = 'CLOSEST_MIX'
      CHARACTER(16), PARAMETER :: PU_SUN_MIX = 'PU_SUN_MIX'

C Pruning method
      CHARACTER(16) :: PRUNING_METHOD
      CHARACTER(16), PARAMETER :: NO_PRUNING = 'NO_PRUNING'
      CHARACTER(16), PARAMETER :: KEEP_CLOSEST = 'KEEP_CLOSEST'
      CHARACTER(16), PARAMETER :: KEEP_OLDEST = 'KEEP_OLDEST'
c$$$ TODO:      CHARACTER(16), PARAMETER :: KEEP_MAX_MIN_AVG = 'KEEP_MAX_MIN_AVG'
c$$$ TODO:      CHARACTER(16), PARAMETER :: PRUNE_SIMILAR = 'PRUNE_SIMILAR'

C Pruning parameters
      INTEGER :: PRUNING_FREQ = 5 ! pruning frequency (number of steps)
C High resolution cells
      INTEGER :: HR_KEEP_IN_CELL ! number of packets to keep after pruning
      INTEGER :: HR_KEEP_TOL    ! prune when NUM_IN_CELL > KEEP_IN_CELL + KEEP_TOL
C Normal resolution cells
      INTEGER :: NR_KEEP_IN_CELL = 2 ! number of packets to keep after pruning
      INTEGER :: NR_KEEP_TOL = 2 ! prune when NUM_IN_CELL > KEEP_IN_CELL + KEEP_TOL

C Sub-grid diffusion parameters
      REAL :: MAX_SGD_FAC = 0.1 ! maximum sub-grid diffusion factor allowed

C Output files and associated arrays
      CHARACTER(16), PARAMETER :: TG_PACKET_1 = 'TG_PACKET_1'
      INTEGER, POINTER :: NEW_PACKETS(:, :, :) ! new packets added to cell during output step

      CHARACTER(16), PARAMETER :: CLS_MIX_1 = 'CLS_MIX_1'
      REAL, POINTER :: CLS_MIX(:, :, :, :) ! packet closest to cell center

      CHARACTER(16), PARAMETER :: AVG_MIX_1 = 'AVG_MIX_1'
      REAL, POINTER :: AVG_MIX(:, :, :, :) ! average over packets in cell

      CHARACTER(16), PARAMETER :: MIN_MIX_1 = 'MIN_MIX_1'
      REAL, POINTER :: MIN_MIX(:, :, :, :) ! minimum over packets in cell

      CHARACTER(16), PARAMETER :: MAX_MIX_1 = 'MAX_MIX_1'
      REAL, POINTER :: MAX_MIX(:, :, :, :) ! maximum over packets in cell

      CHARACTER(16), PARAMETER :: OLD_MIX_1 = 'OLD_MIX_1'
      REAL, POINTER :: OLD_MIX(:, :, :, :) ! oldest packet in cell

      INTEGER, PARAMETER :: NMIX_FILES = 5
      CHARACTER(16) :: MIX_FILE(NMIX_FILES) =
     &     (/ CLS_MIX_1, AVG_MIX_1, MIN_MIX_1, MAX_MIX_1, OLD_MIX_1 /)

C Stuff
      REAL, POINTER :: AVG_LAY(:) ! average number of packets in a cell by layer
      INTEGER, POINTER :: MAX_LAY(:) ! maximum number of packets in a cell by layer
      INTEGER :: MAX_CELL       ! maximum packets in a cell
      INTEGER :: MAX_LOC(3)     ! cell with maximum packets
      INTEGER :: NUM_CELLS      ! number of interior cells (for averages)
      REAL :: SUM_AVG_CELL      ! sum over output step of average packets in a cell

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      CONTAINS

C*******************************************************************************
C FUNCTIONS and SUBROUTINES
C*******************************************************************************

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Initialize packet locations and compositions in the interior of the domain
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION PACKETS_INIT(CGRID, JDATE, JTIME, TSTEP)

      USE GRID_CONF             ! horizontal domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
	  USE M3UTILIO

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_CONST       ! constants
!      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
!      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
!      INCLUDE SUBST_COORD_ID    ! domain coordinate definitions (X3FACE_GD)
!      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

C Result:
      LOGICAL :: PACKETS_INIT

C Arguments:
      REAL, POINTER :: CGRID(:,:,:,:)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C External Functions:
!      INTEGER, EXTERNAL :: ENVINT
!      REAL, EXTERNAL    :: ENVREAL

C File Variables:

C Local Variables:
      INTEGER :: LOGDEV

      CHARACTER(16) :: PNAME = 'PACKETS_INIT'
      CHARACTER(96) :: XMSG = ' '
      CHARACTER(80) :: VARDESC  ! environment variable description

      INTEGER :: STATUS
      INTEGER :: ALLOCSTAT

      INTEGER :: HR_MULT = 2    ! resolution multiplier for high resolution box

      CHARACTER(16) :: HR_BCOL_ECOL = 'HR_BCOL_ECOL'
      CHARACTER(16) :: HR_BROW_EROW = 'HR_BROW_EROW'
      CHARACTER(16) :: VAL_LIST(2)
      INTEGER :: NVARS

      INTEGER :: C, R, L        ! column, row and layer indices
      INTEGER :: I, J           ! loop counters

      REAL(8) :: X, Y, Z        ! packet coordinates (m)
	  REAL(8) :: AVG_CELL

      REAL(8) :: XC, YC         ! coordinates of lower left corner of cell (m)

      REAL(8) :: HR_DX, HR_DY   ! for high resolution box (m)

      TYPE(PacketListMember), POINTER :: ITR
      TYPE(Packet), POINTER :: P

      INTERFACE
         SUBROUTINE GET_ENVLIST(ENV_VAR, NVARS, VAL_LIST)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: ENV_VAR
            INTEGER, INTENT(OUT) :: NVARS
            CHARACTER(16), INTENT(OUT) :: VAL_LIST(:)
         END SUBROUTINE GET_ENVLIST
      END INTERFACE

C-----------------------------------

      PACKETS_INIT = .TRUE.

      LOGDEV = INIT3()

C Get parameters for high resolution box
      VARDESC = 'Multiplier for high resolution box'
      HR_MULT = ENVINT('HR_MULT', VARDESC, HR_MULT, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

      CALL GET_ENVLIST(HR_BCOL_ECOL, NVARS, VAL_LIST)
      IF ( NVARS .NE. 2 ) THEN
         WRITE(LOGDEV, '(/5X, A, A, /5X, A, I3, A/)')
     &        'Error for environment variable ', TRIM(HR_BCOL_ECOL),
     &        'Using default: "0 ', NCOLS+1, '" (all columns)'
         HR_BEGCOL = 0
         HR_ENDCOL = NCOLS + 1
      ELSE
         READ(VAL_LIST(1), '(I4)') HR_BEGCOL
         READ(VAL_LIST(2), '(I4)') HR_ENDCOL
      END IF

      CALL GET_ENVLIST(HR_BROW_EROW, NVARS, VAL_LIST)
      IF ( NVARS .NE. 2 ) THEN
         WRITE(LOGDEV, '(/5X, A, A, /5X, A, I3, A/)')
     &        'Error for environment variable ', TRIM(HR_BROW_EROW),
     &        'Using default: "0 ', NROWS+1, '" (all rows)'
         HR_BEGROW = 0
         HR_ENDROW = NROWS + 1
      ELSE
         READ(VAL_LIST(1), '(I4)') HR_BEGROW
         READ(VAL_LIST(2), '(I4)') HR_ENDROW
      END IF

      VARDESC = 'Top layer in high resolution box'
      HR_TOP = ENVINT('HR_TOP', VARDESC, HR_TOP, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

C Get fill method
      VARDESC = 'Fill method for cells with too few packets'
      CALL ENVSTR('FILL_METHOD', VARDESC, FILL_ALL, FILL_METHOD, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF
      IF ( FILL_METHOD /= NO_FILL .AND.
     &     FILL_METHOD /= FILL_ALL .AND.
     &     FILL_METHOD /= SPARSE_FILL ) THEN
         WRITE(LOGDEV, '(5X, A, A)')
     &        'Invalid choice of FILL_METHOD: ', FILL_METHOD
         XMSG = 'FILL_METHOD must be NO_FILL, FILL_ALL or SPARSE_FILL'
         PACKETS_INIT = .FALSE.
      END IF

C Get pruning method
      VARDESC = 'Pruning method for cells with too many packets'
      CALL ENVSTR('PRUNING_METHOD', VARDESC, KEEP_CLOSEST,
     &  PRUNING_METHOD, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF
      IF ( PRUNING_METHOD /= NO_PRUNING .AND.
     &     PRUNING_METHOD /= KEEP_CLOSEST .AND.
     &     PRUNING_METHOD /= KEEP_OLDEST ) THEN
         WRITE(LOGDEV, '(5X, A, A)')
     &        'Invalid choice of PRUNING_METHOD: ', PRUNING_METHOD
         XMSG = 'PRUNING_METHOD must be ' //
     &        'NO_PRUNING, KEEP_CLOSEST or KEEP_OLDEST'
         PACKETS_INIT = .FALSE.
      END IF

C Get pruning parameters
      VARDESC = 'Pruning frequency (number of steps)'
      PRUNING_FREQ = ENVINT('PRUNING_FREQ', VARDESC, PRUNING_FREQ,
     &  STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

C High resolution cells
      VARDESC = 'Number of packets to keep after pruning' //
     &     ' high resolution cells'
      HR_KEEP_IN_CELL = ENVINT('HR_KEEP_IN_CELL', VARDESC, HR_MULT**2,
     &  STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

      VARDESC = 'Tolerance above HR_KEEP_IN_CELL before a cell is pruned'
      HR_KEEP_TOL = ENVINT('HR_KEEP_TOL', VARDESC, HR_KEEP_IN_CELL,
     &  STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

C Normal resolution cells
      VARDESC = 'Number of packets to keep after pruning' //
     &     ' normal resolution cells'
      NR_KEEP_IN_CELL = ENVINT('NR_KEEP_IN_CELL', VARDESC,
     &  NR_KEEP_IN_CELL, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

      VARDESC = 'Tolerance above NR_KEEP_IN_CELL before a cell is pruned'
      NR_KEEP_TOL = ENVINT('NR_KEEP_TOL', VARDESC, NR_KEEP_TOL, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

C Get sub-grid diffusion parameters
      VARDESC = 'Maximum sub-grid diffusion factor allowed'
      MAX_SGD_FAC = ENVREAL('MAX_SGD_FAC', VARDESC, MAX_SGD_FAC, STATUS)
      IF ( STATUS .NE. 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS .EQ. 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         PACKETS_INIT = .FALSE.
      END IF

      IF ( .NOT. PACKETS_INIT ) RETURN

C Set Eulerian grid data
      NUM_CELLS = NCOLS*NROWS*NLAYS

      IF ( GDTYP_GD == LATGRD3 ) THEN
         DX_GD = DG2M * XCELL_GD
     &        * COS( PI180*( YORIG_GD + YCELL_GD*REAL(GL_NROWS/2, 8) ) ) ! (m)
         DY_GD = DG2M * YCELL_GD ! (m)
      ELSE
         DX_GD = XCELL_GD       ! (m)
         DY_GD = YCELL_GD       ! (m)
      END IF

      XORIG_DOT = XORIG_GD - 0.5 * DX_GD
      YORIG_DOT = YORIG_GD - 0.5 * DY_GD

      LastID_GD = 0

C Allocate memory for data
      ALLOCATE(
     &     PacketsInCell(0:NCOLS+1, 0:NROWS+1, 1:NLAYS),
     &     ResMult(0:NCOLS+1, 0:NROWS+1, 1:NLAYS),
     &     ResDel(0:NCOLS+1, 0:NROWS+1, 1:NLAYS),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating PacketsInCell, ResMult, ResDel'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      ALLOCATE(
     &     AVG_LAY(1:NLAYS),
     &     MAX_LAY(1:NLAYS),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating AVG_LAY, MAX_LAY'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      ALLOCATE(
     &     NEW_PACKETS(1:NCOLS, 1:NROWS, 1:NLAYS),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating packet info arrays'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF
      NEW_PACKETS = 0

      ALLOCATE(
     &     CLS_MIX(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD),
     &     AVG_MIX(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD),
     &     MIN_MIX(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD),
     &     MAX_MIX(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD),
     &     OLD_MIX(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating gridded mixing ratios'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF
      CLS_MIX = 0.0
      AVG_MIX = 0.0
      MIN_MIX = 0.0
      MAX_MIX = 0.0
      OLD_MIX = 0.0

C Set variable resolution and pruning arrays
      HR_DX = DX_GD / REAL(HR_MULT, 8)
      HR_DY = DY_GD / REAL(HR_MULT, 8)

      DO L = 1, HR_TOP
         DO R = 0, NROWS + 1
            DO C = 0, NCOLS + 1
               IF ( C >= HR_BEGCOL .AND. C <= HR_ENDCOL .AND.
     &              R >= HR_BEGROW .AND. R <= HR_ENDROW ) THEN
                  ResMult(C, R, L)%X = HR_MULT
                  ResMult(C, R, L)%Y = HR_MULT
                  ResDel(C, R, L)%X = HR_DX
                  ResDel(C, R, L)%Y = HR_DY
               ELSE
                  ResMult(C, R, L)%X = 1
                  ResMult(C, R, L)%Y = 1
                  ResDel(C, R, L)%X = DX_GD
                  ResDel(C, R, L)%Y = DY_GD
               END IF
            END DO
         END DO
      END DO

      DO L = HR_TOP + 1, NLAYS
         DO R = 0, NROWS + 1
            DO C = 0, NCOLS + 1
               ResMult(C, R, L)%X = 1
               ResMult(C, R, L)%Y = 1
               ResDel(C, R, L)%X = DX_GD
               ResDel(C, R, L)%Y = DY_GD
            END DO
         END DO
      END DO

C Each cell is initialized with the number of packets determined by its
C horizontal resolution multipiers
      DO L = 1, NLAYS
         Z = 0.5 * (X3FACE_GD(L) + X3FACE_GD(L - 1))
         DO R = 0, NROWS + 1
            YC = YORIG_DOT + (R - 1) * DY_GD
            DO C = 0, NCOLS + 1
               XC = XORIG_DOT + (C - 1) * DX_GD

               CALL INITIALIZE_PACKET_LIST(PacketsInCell(C, R, L))

               Y = YC - 0.5 * ResDel(C, R, L)%Y
               DO J = 1, ResMult(C, R, L)%Y
                  Y = Y + ResDel(C, R, L)%Y
                  X = XC - 0.5 * ResDel(C, R, L)%X
                  DO I = 1, ResMult(C, R, L)%X
                     X = X + ResDel(C, R, L)%X
                     CALL NEW_PACKET(JDATE, JTIME, X, Y, Z, P)
                     CALL NEW_PACKET_LIST_MEMBER(JDATE, JTIME, P,
     &                    PacketsInCell(C, R, L))
                  END DO
               END DO

            END DO
         END DO
      END DO

C Initialize the composition of interior packets from CGRID
      DO L = 1, NLAYS
         DO R = 1, NROWS
            DO C = 1, NCOLS
               ITR => PacketsInCell(C, R, L)%Head
               DO
                  IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                  ITR%P%Mix = CGRID(C, R, L, :)
                  ITR => ITR%Next
               END DO
            END DO
         END DO
      END DO

C Test opening existing packet information output file for update
      IF ( .NOT. OPEN3(TG_PACKET_1, FSRDWR3, PNAME) ) THEN
         XMSG = 'Could not open ' // TRIM(TG_PACKET_1) // ' file for update - '
     &        // 'try to open new'
         CALL M3MESG(XMSG)

C Open the file new
         CALL OPEN_PACKET_INFO(JDATE, JTIME, TSTEP(1))
      END IF

      IF ( MYPE .EQ. 0 ) THEN
         IF ( .NOT. CLOSE3(TG_PACKET_1) ) THEN
            XMSG = 'Could not close ' // TRIM(TG_PACKET_1)
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF
      END IF

      AVG_CELL = 0.0
      DO L = 1, NLAYS
         DO R = 1, NROWS
            DO C = 1, NCOLS
               AVG_CELL = AVG_CELL + PacketsInCell(C, R, L)%Count
            END DO
         END DO
      END DO
      AVG_CELL = AVG_CELL / REAL(NUM_CELLS)
      SUM_AVG_CELL = SUM_AVG_CELL + AVG_CELL

      WRITE(LOGDEV, '(5X, A, F6.2/)') 'Avg in cell initial: ', AVG_CELL	  
	  
      END FUNCTION PACKETS_INIT

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Create a new packet and a new packet list member for the cell; the composition
C in the new packet depends on the context in which the packet is added, so Mix
C is set elsewhere
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE NEW_PACKET(JDATE, JTIME, X, Y, Z, P)

      USE CGRID_SPCS            ! CGRID species number and offsets
	  USE GRID_CONF
	  USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
!      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
	  INCLUDE SUBST_CONST
	  INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      REAL(8), INTENT(IN) :: X, Y, Z ! packet coordinates (m)
      TYPE(Packet), POINTER :: P

C Local Variables:
      CHARACTER(16) :: PNAME = 'NEW_PACKET'

      INTEGER :: ALLOCSTAT
      CHARACTER(96) :: XMSG = ' '

C-----------------------------------

C Allocate memory for new packet
      ALLOCATE(P, STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating new packet'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      ALLOCATE(P%Mix(1:NSPCSD), STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating Mix for new packet'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

C Set data for the new packet
      LastID_GD = LastID_GD + 1
      P%ID_GD = LastID_GD

      P%AGE = 0.0

      P%X = X
      P%Y = Y
      P%Z = Z

      END SUBROUTINE NEW_PACKET

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Delete a packet, given as a member of the PacketsInCell list for a cell; free
C associated memory
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DELETE_PACKET(JDATE, JTIME, PLM, C, R, L)

      USE GRID_CONF
	  USE UTILIO_DEFN
	  
	  IMPLICIT NONE

C Includes:
      INCLUDE SUBST_CONST
	  INCLUDE SUBST_FILES_ID
!	  INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      TYPE(PacketListMember), POINTER :: PLM
      INTEGER, INTENT(IN) :: C, R, L ! cell indices

C Local Variables:
      CHARACTER(16) :: PNAME = 'DELETE_PACKET'

      INTEGER :: ALLOCSTAT
      CHARACTER(96) :: XMSG = ' '

      TYPE(Packet), POINTER :: P

C-----------------------------------

      P => PLM%P

      DEALLOCATE(P%Mix, STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure deallocating Mix for packet'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      DEALLOCATE(P, STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure deallocating packet'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      CALL DELETE_PACKET_LIST_MEMBER(JDATE, JTIME, PLM,
     &  PacketsInCell(C, R, L))

      END SUBROUTINE DELETE_PACKET

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Create a new packet list member and add it to a packet list
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE NEW_PACKET_LIST_MEMBER(JDATE, JTIME, P, PACKET_LIST)

      USE GRID_CONF
      USE UTILIO_DEFN
      
      IMPLICIT NONE

C Includes:
!      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
      INCLUDE SUBST_CONST
      INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      TYPE(Packet), POINTER :: P
      TYPE(PacketList) :: PACKET_LIST

C Local Variables:
      CHARACTER(16) :: PNAME = 'NEW_PACKET_LIST_MEMBER'

      INTEGER :: ALLOCSTAT
      CHARACTER(96) :: XMSG = ' '

      TYPE(PacketListMember), POINTER :: PLM

C-----------------------------------

C Allocate memory for new packet list member
      ALLOCATE(PLM, STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating new packet list member'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

C Set data
      PLM%P => P

C Add to specified packet list
      CALL ADD_PACKET_LIST_MEMBER(PLM, PACKET_LIST)

      END SUBROUTINE NEW_PACKET_LIST_MEMBER

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Delete a packet list member from a packet list - free the memory for the PLM,
C but keep the packet
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DELETE_PACKET_LIST_MEMBER(JDATE, JTIME, PLM,
     &  PACKET_LIST)

      USE GRID_CONF
      USE UTILIO_DEFN
      
      IMPLICIT NONE

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      TYPE(PacketListMember), POINTER :: PLM
      TYPE(PacketList) :: PACKET_LIST

C Local Variables:
      CHARACTER(16) :: PNAME = 'DELETE_PACKET_LIST_MEMBER'

      INTEGER :: ALLOCSTAT
      CHARACTER(96) :: XMSG = ' '

C-----------------------------------

      CALL REMOVE_PACKET_LIST_MEMBER(PLM, PACKET_LIST)

      DEALLOCATE(PLM, STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure deallocating packet list member'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      END SUBROUTINE DELETE_PACKET_LIST_MEMBER

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Add an existing PLM to tail of a packet list
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ADD_PACKET_LIST_MEMBER(PLM, PACKET_LIST)

      USE GRID_CONF
      USE UTILIO_DEFN
           
      IMPLICIT NONE

C Arguments:
      TYPE(PacketListMember), POINTER :: PLM
      TYPE(PacketList) :: PACKET_LIST

C-----------------------------------
!	NULLIFY (PACKET_LIST%Head)
      IF ( .NOT. ASSOCIATED(PACKET_LIST%Head) ) THEN
         PACKET_LIST%Head => PLM
         NULLIFY(PLM%Prev)
      ELSE
         PACKET_LIST%Tail%Next => PLM
         PLM%Prev => PACKET_LIST%Tail
      END IF
      NULLIFY(PLM%Next)
      PACKET_LIST%Tail => PLM

      PACKET_LIST%Count = PACKET_LIST%Count + 1

      END SUBROUTINE ADD_PACKET_LIST_MEMBER

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Remove a packet list member from a packet list - keep both the PLM and its
C packet
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE REMOVE_PACKET_LIST_MEMBER(PLM, PACKET_LIST)

      IMPLICIT NONE

C Arguments:
      TYPE(PacketListMember), POINTER :: PLM
      TYPE(PacketList) :: PACKET_LIST

C-----------------------------------

      IF ( .NOT. ASSOCIATED(PLM%Prev) ) THEN ! PLM is PACKET_LIST%Head
         IF ( ASSOCIATED(PLM%Next) ) THEN
            PACKET_LIST%Head => PLM%Next
            NULLIFY(PACKET_LIST%Head%Prev)
         ELSE                   ! PLM is also PACKET_LIST%Tail
            NULLIFY(PACKET_LIST%Head, PACKET_LIST%Tail)
         END IF
      ELSE IF ( .NOT. ASSOCIATED(PLM%Next) ) THEN ! PLM is PACKET_LIST%Tail
         PACKET_LIST%Tail => PLM%Prev
         NULLIFY(PACKET_LIST%Tail%Next)
      ELSE
         PLM%Prev%Next => PLM%Next
         PLM%Next%Prev => PLM%Prev
      END IF

      PACKET_LIST%Count = PACKET_LIST%Count - 1

      END SUBROUTINE REMOVE_PACKET_LIST_MEMBER

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Initialize a packet list
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE INITIALIZE_PACKET_LIST(PACKET_LIST)

      IMPLICIT NONE

C Arguments:
      TYPE(PacketList) :: PACKET_LIST

C-----------------------------------

      PACKET_LIST%Count = 0
      NULLIFY(PACKET_LIST%Head, PACKET_LIST%Tail)

      END SUBROUTINE INITIALIZE_PACKET_LIST

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Extend a packet list with a list of additional packets by creating new members
C (do not change the list of additional packets or any of its members)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE EXTEND_PACKET_LIST(JDATE, JTIME, PACKET_LIST,
     &  ADDITIONAL_PACKETS)

      IMPLICIT NONE

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      TYPE(PacketList) :: PACKET_LIST, ADDITIONAL_PACKETS

C Local Variables:
      TYPE(PacketListMember), POINTER :: ITR

C-----------------------------------

      ITR => ADDITIONAL_PACKETS%Head
      DO
         IF ( .NOT. ASSOCIATED(ITR) ) EXIT
         CALL NEW_PACKET_LIST_MEMBER(JDATE, JTIME, ITR%P, PACKET_LIST)
         ITR => ITR%Next
      END DO

      END SUBROUTINE EXTEND_PACKET_LIST

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Delete a packet list
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DELETE_PACKET_LIST(JDATE, JTIME, PACKET_LIST)

      IMPLICIT NONE

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      TYPE(PacketList) :: PACKET_LIST

C-----------------------------------

      DO
         IF ( .NOT. ASSOCIATED(PACKET_LIST%Head) ) EXIT

         CALL DELETE_PACKET_LIST_MEMBER(JDATE, JTIME, PACKET_LIST%Head,
     &     PACKET_LIST)
      END DO

      END SUBROUTINE DELETE_PACKET_LIST

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Given a point in the x-y plane and a packet list, search for a packet in the
C list that is closer to the point (in the x-y plane) than the closest packet
C previously found
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GET_CLOSER_TO_POINT(X, Y, PACKET_LIST, D2_LEAST,
     &  CLOSEST)

      IMPLICIT NONE

C Arguments:
      REAL(8), INTENT(IN) :: X, Y ! coordinates of point (m)
      TYPE(PacketList) :: PACKET_LIST
      REAL(8), INTENT(IN OUT) :: D2_LEAST ! distance (squared) to point
      TYPE(Packet), POINTER :: CLOSEST

C Local Variables:
      REAL(8) :: D2             ! distance (squared) to center of cell

      TYPE(PacketListMember), POINTER :: ITR
      TYPE(Packet), POINTER :: P

C-----------------------------------

      ITR => PACKET_LIST%Head
      DO
         IF ( .NOT. ASSOCIATED(ITR) ) EXIT
         P => ITR%P

         D2 = (P%X - X)**2 + (P%Y - Y)**2
         IF ( D2 < D2_LEAST ) THEN
            D2_LEAST = D2
            CLOSEST => P
         END IF

         ITR => ITR%Next
      END DO

      END SUBROUTINE GET_CLOSER_TO_POINT

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Given a cell, find a packet in the same layer of the domain (including the
C boundary) that is closest to the center of the cell in the x-y plane;
C optionally, get list of nearby packets
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION CLOSEST_PACKET(JDATE, JTIME, C,R,L, NCX, NEARBY_PACKETS)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Result:
      TYPE(Packet), POINTER :: CLOSEST_PACKET

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
C      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: C, R, L ! column, row and layer indices
      INTEGER, INTENT(IN), OPTIONAL :: NCX ! max number of cells from given cell to search
      TYPE(PacketList), OPTIONAL :: NEARBY_PACKETS
	  TYPE(PacketListMember), POINTER :: PLM

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

      CHARACTER(16) :: PNAME = 'CLOSEST_PACKET'
      CHARACTER(96) :: XMSG = ' '

      REAL(8) :: XC, YC         ! coordinates of the center of the cell (m)
      REAL(8) :: D2_LEAST       ! distance (squared) to center of cell

      INTEGER :: CS, RS         ! column and row of cell to search
      INTEGER :: NC             ! number of cells from given cell to search
      INTEGER :: NCMAX          ! upper bound for NC

C-----------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
      END IF

      NULLIFY(CLOSEST_PACKET)

      IF ( PRESENT(NCX) ) THEN
         NCMAX = NCX
      ELSE
         NCMAX = 3
      END IF

      IF ( PRESENT(NEARBY_PACKETS) )
     &     CALL INITIALIZE_PACKET_LIST(NEARBY_PACKETS)

C GOTCHA: With the strategy below, if the packet chosen as "closest" is near
C a cell corner, then there may be a packet outside the searched cells that
C is actually closer

C If the cell contains exactly one packet and the list of nearby packets isn't
C needed, use that packet
      IF ( PacketsInCell(C, R, L)%Count == 1 .AND.
     &     .NOT. PRESENT(NEARBY_PACKETS) ) THEN
         CLOSEST_PACKET => PacketsInCell(C, R, L)%Head%P

C Otherwise
      ELSE
         D2_LEAST = 1.0E30

         XC = XORIG_GD + (C - 1) * XCELL_GD
         YC = YORIG_GD + (R - 1) * YCELL_GD

         IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
            CALL GET_CLOSER_TO_POINT(XC, YC, PacketsInCell(C, R, L),
     &           D2_LEAST, CLOSEST_PACKET)
            IF ( PRESENT(NEARBY_PACKETS) ) THEN
               CALL EXTEND_PACKET_LIST(JDATE, JTIME, NEARBY_PACKETS,
     &              PacketsInCell(C, R, L))
            ELSE IF ( ASSOCIATED(CLOSEST_PACKET) ) THEN
               RETURN
            END IF 
	 END IF

         DO NC = 1, NCMAX
            IF ( R >= NC ) THEN
               RS = R - NC      ! Southern edge
               DO CS = MAX(C-NC+1, 1), MIN(C+NC, NCOLS+1)
                  IF ( PacketsInCell(CS, RS, L)%Count > 0 ) THEN
                     CALL GET_CLOSER_TO_POINT(XC, YC,
     &                 PacketsInCell(CS, RS, L),
     &                    D2_LEAST, CLOSEST_PACKET)
                     IF ( PRESENT(NEARBY_PACKETS) ) THEN
                        CALL EXTEND_PACKET_LIST(JDATE, JTIME,
     &                    NEARBY_PACKETS,
     &                       PacketsInCell(CS, RS, L))
                     END IF
                  END IF
               END DO
            END IF

            IF ( C <= NCOLS - NC + 1 ) THEN
               CS = C + NC      ! Eastern edge
               DO RS = MAX(R-NC+1, 1), MIN(R+NC, NROWS+1)
                  IF ( PacketsInCell(CS, RS, L)%Count > 0 ) THEN
                     CALL GET_CLOSER_TO_POINT(XC, YC,
     &                 PacketsInCell(CS, RS, L),
     &                    D2_LEAST, CLOSEST_PACKET)
                     IF ( PRESENT(NEARBY_PACKETS) ) THEN
                        CALL EXTEND_PACKET_LIST(JDATE, JTIME,
     &                    NEARBY_PACKETS,
     &                       PacketsInCell(CS, RS, L))
                     END IF
                  END IF
               END DO
            END IF

            IF ( R <= NROWS - NC + 1 ) THEN
               RS = R + NC      ! Northern edge
               DO CS = MAX(C-NC, 0), MIN(C+NC-1, NCOLS)
                  IF ( PacketsInCell(CS, RS, L)%Count > 0 ) THEN
                     CALL GET_CLOSER_TO_POINT(XC, YC,
     &                 PacketsInCell(CS, RS, L),
     &                    D2_LEAST, CLOSEST_PACKET)
                     IF ( PRESENT(NEARBY_PACKETS) ) THEN
                        CALL EXTEND_PACKET_LIST(JDATE, JTIME,
     &                    NEARBY_PACKETS,
     &                       PacketsInCell(CS, RS, L))
                     END IF
                  END IF
               END DO
            END IF

            IF ( C >= NC ) THEN
               CS = C - NC      ! Western edge
               DO RS = MAX(R-NC, 0), MIN(R+NC-1, NROWS)
                  IF ( PacketsInCell(CS, RS, L)%Count > 0 ) THEN
                     CALL GET_CLOSER_TO_POINT(XC, YC,
     &                 PacketsInCell(CS, RS, L),
     &                    D2_LEAST, CLOSEST_PACKET)
                     IF ( PRESENT(NEARBY_PACKETS) ) THEN
                        CALL EXTEND_PACKET_LIST(JDATE, JTIME,
     &                    NEARBY_PACKETS,
     &                       PacketsInCell(CS, RS, L))
                     END IF
                  END IF
               END DO
            END IF

            IF ( ASSOCIATED(CLOSEST_PACKET) .AND.
     &           .NOT. PRESENT(NEARBY_PACKETS) ) RETURN
         END DO
      END IF

      IF ( .NOT. ASSOCIATED(CLOSEST_PACKET) ) THEN
         WRITE(LOGDEV, '(/5X, A, 3I4)') 'Cell (C, R, L): ', C, R, L
		 XMSG = 'Search for closest packet failed'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      END FUNCTION CLOSEST_PACKET

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Update one boundary cell:
C   - if the cell is empty, create new packet(s) in the cell
C   - update Mix for packets in the cell
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE UPDATE_BOUNDARY_CELL(JDATE, JTIME, C, R, L, K, BCON)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
C      INCLUDE SUBST_COORD_ID    ! domain coordinate definitions (X3FACE_GD)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: C, R, L ! column, row and layer indices
      INTEGER, INTENT(IN) :: K  ! horizontal cell index in BCON array
      REAL, INTENT(IN)    :: BCON(:, :, :) ! boundary conditions

C Local Variables:
      CHARACTER(16) :: PNAME = 'UPDATE_BOUNDARY_CELL'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: I, J           ! loop counters

      REAL(8) :: X, Y, Z        ! packet coordinates (m)

      REAL(8) :: XC, YC         ! coordinates of lower left corner of cell (m)
      REAL(8) :: XB, YB         ! base coordinates for new boundary packets (m)

      TYPE(PacketListMember), POINTER :: ITR
      TYPE(Packet), POINTER :: P

C-----------------------------------

      IF ( PacketsInCell(C, R, L)%Count == 0 ) THEN
         XC = XORIG_DOT + (C - 1) * DX_GD
         YC = YORIG_DOT + (R - 1) * DY_GD
         Z = 0.5 * (X3FACE_GD(L) + X3FACE_GD(L - 1))

         IF ( C == 0 ) THEN
            XB = XC - 0.95 * ResDel(C, R, L)%X
         ELSE IF ( C == NCOLS+1 ) THEN
            XB = XC - 0.05 * ResDel(C, R, L)%X
         ELSE
            XB = XC - 0.5 * ResDel(C, R, L)%X
         END IF

         IF ( R == 0 ) THEN
            YB = YC - 0.95 * ResDel(C, R, L)%Y
         ELSE IF ( R == NROWS+1 ) THEN
            YB = YC - 0.05 * ResDel(C, R, L)%Y
         ELSE
            YB = YC - 0.5 * ResDel(C, R, L)%Y
         END IF

         Y = YB
         DO J = 1, ResMult(C, R, L)%Y
            Y = Y + ResDel(C, R, L)%Y
            X = XB
            DO I = 1, ResMult(C, R, L)%X
               X = X + ResDel(C, R, L)%X
               CALL NEW_PACKET(JDATE, JTIME, X, Y, Z, P)
               CALL NEW_PACKET_LIST_MEMBER(JDATE, JTIME, P,
     &              PacketsInCell(C, R, L))
            END DO
         END DO
      END IF

C The TG version of RDBCON keeps composition in mixing ratio form and BCON is
C dimensioned for the full set of CGRID species (including, to avoid expanded
C code changes, the nasty slot for advected RHOJ)
      ITR => PacketsInCell(C, R, L)%Head
      DO
         IF ( .NOT. ASSOCIATED(ITR) ) EXIT
         ITR%P%Mix = BCON(K, L, :)
         ITR => ITR%Next
      END DO

      END SUBROUTINE UPDATE_BOUNDARY_CELL

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Given the (X, Y, Z) coordinates of a point, find the cell the point is in,
C i.e., (C, R, L), and its fractional location in the cell
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GET_CELL(X, Y, Z, C, R, L, FR_X, FR_Y, FR_Z)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
C      INCLUDE SUBST_COORD_ID    ! domain coordinate definitions (X3FACE_GD)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      REAL(8), INTENT(IN)     :: X, Y, Z ! coordinates of point (m)
      INTEGER, INTENT(OUT)    :: C, R ! column and row indices of cell
      INTEGER, INTENT(IN OUT) :: L ! layer index of cell
      REAL(8), INTENT(OUT), OPTIONAL :: FR_X, FR_Y, FR_Z ! fractional location in cell

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      REAL :: IDZ(1:NLAYS) ! 1/(cell depth)

      INTEGER :: LL             ! layer

      REAL(8) :: FR_C, FR_R, FR_L ! fractional indices of DOT cell

C-----------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         DO LL = 1, NLAYS
            IDZ(LL) = 1.0 / (X3FACE_GD(LL) - X3FACE_GD(LL-1))
         END DO
      END IF                    ! FIRSTIME

      FR_C = ((X - XORIG_DOT) / DX_GD) + 1
      C = FLOOR(FR_C)
      IF ( C == NCOLS + 1 .AND. FR_C == REAL(C, 8) ) C = NCOLS
      IF ( PRESENT(FR_X) ) FR_X = FR_C - C

      FR_R = ((Y - YORIG_DOT) / DY_GD) + 1
      R = FLOOR(FR_R)
      IF ( R == NROWS + 1 .AND. FR_R == REAL(R, 8) ) R = NROWS
      IF ( PRESENT(FR_Y) ) FR_Y = FR_R - R

C Note that the input layer is used to start the search for the current layer
      IF ( Z < X3FACE_GD(L) ) THEN
         DO LL = L, 1, -1
            IF ( Z >= X3FACE_GD(LL-1) ) EXIT
         END DO
      ELSE
         DO LL = L + 1, NLAYS
            IF ( Z < X3FACE_GD(LL) ) EXIT
         END DO
      END IF
      L = LL
      IF ( PRESENT(FR_Z) ) FR_Z = (Z - X3FACE_GD(L-1)) * IDZ(L)

      END SUBROUTINE GET_CELL

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Interpolate (in space) wind velocity to a point
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE INTERP_WIND(UWIND, VWIND, WWIND, LC, X, Y, Z, U, V, W)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Arguments:
      REAL, INTENT(IN)     :: UWIND(1:, 0:, 1:) ! components of wind velocity
      REAL, INTENT(IN)     :: VWIND(0:, 1:, 1:)
      REAL, INTENT(IN)     :: WWIND(:, :, :)
      INTEGER, INTENT(IN)  :: LC ! current layer of point
      REAL(8), INTENT(IN)  :: X, Y, Z ! coordinates of point (m)
      REAL(8), INTENT(OUT) :: U, V, W ! components of velocity at point

C Local Variables:
      INTEGER :: C_DOT, R_DOT   ! column and row indices of DOT point cell
      INTEGER :: C, R, L        ! column, row and layer indices for interpolation

      REAL(8) :: FR_X_DOT, FR_Y_DOT ! fractional location in DOT cell
      REAL(8) :: FR_X, FR_Y, FR_Z ! fractional location for interpolation

C-----------------------------------

      L = LC
      CALL GET_CELL(X, Y, Z, C_DOT, R_DOT, L, FR_X_DOT, FR_Y_DOT, FR_Z)

C Project the values returned by GET_CELL onto the closure of the interior of
C the domain so that interpolation won't require unavailable wind components
      IF ( C_DOT < 1 ) THEN
         C_DOT = 1
         FR_X_DOT = 0.0
      ELSE IF ( C_DOT > NCOLS ) THEN
         C_DOT = NCOLS
         FR_X_DOT = 1.0
      END IF

      IF ( R_DOT < 1 ) THEN
         R_DOT = 1
         FR_Y_DOT = 0.0
      ELSE IF ( R_DOT > NROWS ) THEN
         R_DOT = NROWS
         FR_Y_DOT = 1.0
      END IF

C Bilinear interpolation of UWIND from X direction flux points
      C = C_DOT
      FR_X = FR_X_DOT
      IF ( FR_Y_DOT >= 0.5 ) THEN
         R = R_DOT
         FR_Y = FR_Y_DOT - 0.5
      ELSE
         R = R_DOT - 1
         FR_Y = FR_Y_DOT + 0.5
      END IF
      U = (1.0 - FR_Y) *
     &     ((1.0 - FR_X) * UWIND(C, R, L) + FR_X * UWIND(C+1, R, L)) +
     &     FR_Y *
     &     ((1.0 - FR_X) * UWIND(C, R+1, L) + FR_X * UWIND(C+1, R+1, L))
	 

C Bilinear interpolation of VWIND from Y direction flux points
      IF ( FR_X_DOT >= 0.5 ) THEN
         C = C_DOT
         FR_X = FR_X_DOT - 0.5
      ELSE
         C = C_DOT - 1
         FR_X = FR_X_DOT + 0.5
      END IF
      R = R_DOT
      FR_Y = FR_Y_DOT
      V = (1.0 - FR_Y) *
     &     ((1.0 - FR_X) * VWIND(C, R, L) + FR_X * VWIND(C+1, R, L)) +
     &     FR_Y *
     &     ((1.0 - FR_X) * VWIND(C, R+1, L) + FR_X * VWIND(C+1, R+1, L))
	 

C Linear interpolation of WWIND from Z direction flux points
      C = C_DOT
      R = R_DOT
!      W = (1.0 - FR_Z) * WWIND(C, R, L) + FR_Z * WWIND(C, R, L+1)
		W = WWIND(C,R,L)
      END SUBROUTINE INTERP_WIND

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C For each packet (including boundary packets), update (U, V, W), (X, Y, Z),
C (C, R, L) and PacketsInCell; delete packets that leave the interior of the
C domain or the boundary through the outer edge of the boundary
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE UPDATE_POSITION(JDATE, JTIME, ASTEP)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
C      INCLUDE SUBST_COORD_ID    ! domain coordinate definitions (X3FACE_GD)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: ASTEP ! advection time step (HHMMSS)

C External Functions:
!      INTEGER, EXTERNAL :: TIME2SEC

C File Variables:
      REAL :: UWIND(1:NCOLS+1, 0:NROWS+1, 1:NLAYS) ! components of wind velocity
      REAL :: VWIND(0:NCOLS+1, 1:NROWS+1, 1:NLAYS)
      REAL :: WWIND(1:NCOLS, 1:NROWS, 1:NLAYS+1)
	  REAL :: HDIV(1:NCOLS+1, 0:NROWS+1, 1:NLAYS)
!	  REAL :: BCON(:,:,:) ! boundary conditions
C Local Variables:
      CHARACTER(16) :: PNAME = 'UPDATE_POSITION'
      CHARACTER(96) :: XMSG = ' '

      REAL(8) :: DT             ! advection time step (sec)

      INTEGER :: C, R, L, M        ! column, row and layer indices
      INTEGER :: CU, RU, LU     ! updated column, row and layer indices

      REAL(8) :: U, V, W        ! velocity components (m/sec)
      REAL(8) :: X, Y, Z        ! point coordinates (m)

      TYPE(PacketListMember), POINTER :: ITR, PLM
      TYPE(Packet), POINTER :: P
	  INTEGER, SAVE :: LOGDEV

      INTERFACE
         ! SUBROUTINE HCONTVEL(JDATE, JTIME, ASTEP, UWIND, VWIND)
            ! IMPLICIT NONE
            ! INTEGER, INTENT(IN) :: JDATE, JTIME, ASTEP
            ! REAL, INTENT(OUT)   :: UWIND(1:, 0:, 1:)
            ! REAL, INTENT(OUT)   :: VWIND(0:, 1:, 1:)
         ! END SUBROUTINE HCONTVEL
		 SUBROUTINE HVELOC( JDATE, JTIME, ASTEP, UWIND, VWIND, HDIV )
            INTEGER, INTENT( IN )  :: JDATE, JTIME, ASTEP
            REAL,    INTENT( IN ) :: UWIND( :,:,: ), VWIND( :,:,: )
            REAL,    INTENT( IN ) :: HDIV( :,:,: )
         END SUBROUTINE HVELOC
         SUBROUTINE VCONTVEL(JDATE, JTIME, ASTEP, WWIND)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: JDATE, JTIME, ASTEP
            REAL, INTENT(IN)   :: WWIND(:, :, :)
         END SUBROUTINE VCONTVEL
		 SUBROUTINE RDBCON ( JDATE, JTIME, ASTEP)
            INTEGER, INTENT( IN )  :: JDATE, JTIME, ASTEP
            ! REAL,    INTENT( OUT ) :: BCON2( :,: )
            ! LOGICAL, INTENT( INOUT ) :: LSTAT
         END SUBROUTINE RDBCON
      END INTERFACE

C-----------------------------------

C Update Mix for packets in the boundary cells
      CALL RDBCON(JDATE, JTIME, ASTEP)

C Get the 3D contravariant velocity components
!      CALL HCONTVEL(JDATE, JTIME, ASTEP, UWIND, VWIND)
	  CALL HVELOC( JDATE, JTIME, ASTEP, UWIND, VWIND, HDIV )

C Extrapolate UWIND to rows 0 and NROWS+1 (for use when interpolating the
C contravariant velocity to a location near the Southern or Northern edge
C of the domain)
C NOTE: This wouldn't be necessary if MCIP wrote UHAT_JD at dot points
      UWIND( :,0,: ) = 2.0 * UWIND( :,1,: ) - UWIND( :,2,: )
      UWIND( :,NROWS+1,: ) = 2.0 * UWIND( :,NROWS,: ) - UWIND( :,NROWS-1,: )	  

C Extrapolate VWIND to columns 0 and NCOLS+1 (for use when interpolating the
C contravariant velocity to a location near the Western or Eastern edge
C of the domain)
C NOTE: This wouldn't be necessary if MCIP wrote VHAT_JD at dot points
      VWIND( 0,:,: ) = 2.0 * VWIND( 1,:,: ) - VWIND( 2,:,: )
      VWIND( NCOLS+1,:,: ) = 2.0 * VWIND( NCOLS,:,: ) - VWIND( NCOLS-1,:,: )
  	  
C no need for vcontvel, in v5.2 hcontvel does all the work
      CALL VCONTVEL(JDATE, JTIME, ASTEP, WWIND)
	  M = 0.0
	DO L = 1, NLAYS + 1
	DO R = 1, NROWS
	DO C = 1, NCOLS	
	IF ( abs(WWIND(C,R,L)) .GE. 10.0 ) THEN
	M = 1.0
	WWIND(:,:,:) = 0.0
	PRINT *,"GOTCHA",C,R,L,WWIND(C,R,L)
	END IF	
	END DO
	END DO	
	END DO
	
	IF ( M .EQ. 1.0 ) THEN
	STOP
	END IF
!	  UWIND(:,:,:) = 1.0
!	  VWIND(:,:,:) = 1.0
!	  WWIND(:,:,:) = 1e-30 * WWIND(:,:,:)
!	  WWIND(:,:,:) = 0.0
!	  WWIND(0,:,:) = 0.0
	  
      DT = REAL(TIME2SEC(ASTEP), 8)

      DO L = 1, NLAYS
         DO R = 0, NROWS + 1
            DO C = 0, NCOLS + 1
               ITR => PacketsInCell(C, R, L)%Head
               DO
                  IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                  P => ITR%P

C Update age
                  P%AGE = P%AGE + DT

C Update velocity in boundary cells
C
C NOTES:
C 1) Corner velocities are initialized as zero and should remain zero as
C    long as only the component normal to the boundary is used
C 2) For outflow boundary cells, there is no point in a cycle that has packets
C    blowing out of the exterior boundary and then being replaced by new
C    packets, so boundary packets remain stationary until they move towards
C    the interior
                  U = 0.0
                  V = 0.0
                  W = 0.0
                  IF ( C < 1 ) THEN
                     IF ( 1 <= R .AND. R <= NROWS )
     &                    U = MAX(UWIND(1, R, L), 0.0)
                  ELSE IF ( C > NCOLS ) THEN
                     IF ( 1 <= R .AND. R <= NROWS )
     &                    U = MIN(UWIND(NCOLS + 1, R, L), 0.0)
                  ELSE IF ( R < 1 ) THEN
                     IF ( 1 <= C .AND. C <= NCOLS )
     &                    V = MAX(VWIND(C, 1, L), 0.0)
                  ELSE IF ( R > NROWS ) THEN
                     IF ( 1 <= C .AND. C <= NCOLS )
     &                    V = MIN(VWIND(C, NROWS + 1, L), 0.0)

C Update velocity in interior cells using a simple predictor-corrector method
C based on the assumption that, locally, wind is linear in space and constant
C in time
                  ELSE
                     CALL INTERP_WIND(UWIND, VWIND, WWIND,
     &                    L, P%X, P%Y, P%Z, U, V, W)

                     X = P%X + 0.5 * DT * U
                     Y = P%Y + 0.5 * DT * V
                     Z = P%Z + 0.5 * DT * W
					 
                     CALL INTERP_WIND(UWIND, VWIND, WWIND,
     &                    L, X, Y, Z, U, V, W)
                  END IF

C Update location
!apouyaei need change
                  P%X = P%X + DT * U
                  P%Y = P%Y + DT * V
                  P%Z = P%Z + DT * W
                  IF ( P%Z <= X3FACE_GD(0) .OR.
     &              P%Z >= X3FACE_GD(NLAYS) ) THEN
                     XMSG = 'Should not be possible to exit domain vertically'
					PRINT *, "out of domain",L,R,C,WWIND(C,R,L),W
!					P%Z = P%Z - DT * W
!					 CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
                  END IF

                  ITR => ITR%Next
               END DO
            END DO
         END DO
      END DO

C Easy to change cells or delete a packet if we loop through the cells and for
C each cell its PacketsInCell list
!		LU=0.0
!		RU=0.0
!		CU=0.0
      DO L = 1, NLAYS
         DO R = 0, NROWS + 1
            DO C = 0, NCOLS + 1
			LU=0.0
			RU=0.0
			CU=0.0
               ITR => PacketsInCell(C, R, L)%Head
               DO
                  IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                  P => ITR%P

                  PLM => ITR
                  ITR => ITR%Next

C Get current cell for the packet
                  CALL GET_CELL(P%X, P%Y, P%Z, CU, RU, LU)
				  

C Delete the packet if it is outside the boundary or it has left the interior
                  IF ( ( CU < 0 .OR. CU > NCOLS+1 .OR.
     &                 RU < 0 .OR. RU > NROWS+1 ) .OR.
     &                 ( ( C >= 1 .AND. C <= NCOLS .AND.
     &                 R >= 1 .AND. R <= NROWS ) .AND.
     &                 ( CU < 1 .OR. CU > NCOLS .OR.
     &                 RU < 1 .OR. RU > NROWS )).OR.
     &                 ( LU < 0 .OR. LU > NLAYS ) ) THEN
!					 PRINT *, "DELET_PACKET VERTICALLY", CU, RU, LU
                     CALL DELETE_PACKET(JDATE, JTIME, PLM, C, R, L)

C Otherwise, if the packet has changed cells, move it from the list of packets
C in its previous cell to the list for its new cell
                  ELSE IF ( CU /= C .OR. RU /= R .OR. LU /= L ) THEN
!                     PRINT *, "CHANGED CELLS", CU, RU, LU
					 CALL REMOVE_PACKET_LIST_MEMBER(PLM,
     &                 PacketsInCell(C, R, L))
                     CALL ADD_PACKET_LIST_MEMBER(PLM,
     &                 PacketsInCell(CU, RU, LU))
                  END IF
               END DO
            END DO
         END DO
      END DO

      END SUBROUTINE UPDATE_POSITION

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Fill empty interior cells by adding a packet at the center of the cell
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FILL_EMPTY_CELLS(JDATE, JTIME)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
C      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
C      INCLUDE SUBST_COORD_ID    ! domain coordinate definitions (X3FACE_GD)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

      CHARACTER(16) :: PNAME = 'FILL_EMPTY_CELLS'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: NUM_ADDED      ! number of interior packets added

      INTEGER :: C, R, L        ! column, row and layer indices
      INTEGER :: I, J           ! loop counters

      LOGICAL :: SKIP_CELL

      REAL(8) :: X, Y, Z        ! packet coordinates (m)

      TYPE(Packet), POINTER :: P, CLOSEST

C-----------------------------------

      IF ( FILL_METHOD == NO_FILL ) RETURN

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         LOGDEV = INIT3()
      END IF                    ! FIRSTIME

      NUM_ADDED = 0

C Fill empty cells
      SELECT CASE ( FILL_METHOD )

      CASE ( FILL_ALL )

C Never allow any interior cell to be empty
         DO L = 1, NLAYS
            Z = 0.5 * (X3FACE_GD(L) + X3FACE_GD(L - 1))
            DO R = 1, NROWS
               Y = YORIG_GD + (R - 1) * YCELL_GD
               DO C = 1, NCOLS
                  X = XORIG_GD + (C - 1) * XCELL_GD
                  IF ( PacketsInCell(C, R, L)%Count == 0 ) THEN
                     CALL NEW_PACKET(JDATE, JTIME, X, Y, Z, P)
                     CALL SET_MIX(JDATE, JTIME, P, C, R, L, CLOSEST_MIX)
                     CALL NEW_PACKET_LIST_MEMBER(JDATE, JTIME, P,
     &                    PacketsInCell(C, R, L))
!                     CALL SET_MIX(JDATE, JTIME, P, C, R, L, CLOSEST_MIX)
                     NUM_ADDED = NUM_ADDED + 1
                     NEW_PACKETS(C, R, L) = NEW_PACKETS(C, R, L) + 1
                  END IF
               END DO
            END DO
         END DO

      CASE ( SPARSE_FILL )

C Never allow interior high resolution cells to be empty
         DO L = 1, HR_TOP
            Z = 0.5 * (X3FACE_GD(L) + X3FACE_GD(L - 1))
            DO R = MAX(HR_BEGROW, 1), MIN(HR_ENDROW, NROWS)
               Y = YORIG_GD + (R - 1) * YCELL_GD
               DO C = MAX(HR_BEGCOL, 1), MIN(HR_ENDCOL, NCOLS)
                  X = XORIG_GD + (C - 1) * XCELL_GD
                  IF ( PacketsInCell(C, R, L)%Count == 0 ) THEN
                     CALL NEW_PACKET(JDATE, JTIME, X, Y, Z, P)
                     CALL SET_MIX(JDATE,JTIME,P,C,R,L,CLOSEST_MIX)
                     CALL NEW_PACKET_LIST_MEMBER(JDATE, JTIME, P,
     &                    PacketsInCell(C, R, L))
                     NUM_ADDED = NUM_ADDED + 1
                     NEW_PACKETS(C, R, L) = NEW_PACKETS(C, R, L) + 1
                  END IF
               END DO
            END DO
         END DO

C Add a packet at the center of every interior low resolution cell that is empty
C and is surrounded (in its layer) by empty cells
         DO L = 1, NLAYS
            Z = 0.5 * (X3FACE_GD(L) + X3FACE_GD(L - 1))
            DO R = 1, NROWS
               Y = YORIG_GD + (R - 1) * YCELL_GD
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count == 0 ) THEN
                     X = XORIG_GD + (C - 1) * XCELL_GD

         SKIP_CHECK: DO J = R - 1, R + 1
                        DO I = C - 1, C + 1
                           SKIP_CELL = (PacketsInCell(I, J, L)%Count>0)
                           IF ( SKIP_CELL ) EXIT SKIP_CHECK
                        END DO
                     END DO SKIP_CHECK
                     IF ( SKIP_CELL ) CYCLE

                     CALL NEW_PACKET(JDATE, JTIME, X, Y, Z, P)
                     CALL SET_MIX(JDATE, JTIME, P, C, R, L, CLOSEST_MIX)
                     CALL NEW_PACKET_LIST_MEMBER(JDATE, JTIME, P,
     &                    PacketsInCell(C, R, L))
                     NUM_ADDED = NUM_ADDED + 1
                     NEW_PACKETS(C, R, L) = NEW_PACKETS(C, R, L) + 1
                  END IF
               END DO
            END DO
         END DO

      CASE DEFAULT
         XMSG =  '"' // TRIM(FILL_METHOD) // '"' //
     &        ' is not a valid selection'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)

      END SELECT

      WRITE(LOGDEV, '(5X, A, I6)') 'Added to interior during step: ', NUM_ADDED

      END SUBROUTINE FILL_EMPTY_CELLS

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Prune cells with too many packets
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PRUNE_CELLS(JDATE, JTIME)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Pointer to a PacketListMember
      TYPE :: PacketListMemberPtr
         TYPE(PacketListMember), POINTER :: PLM
      END TYPE PacketListMemberPtr

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
C      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

      INTEGER, SAVE :: PRUNING_STEP

      CHARACTER(16) :: PNAME = 'PRUNE_CELLS'

      INTEGER :: ALLOCSTAT
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: KEEP_IN_CELL   ! number of packets to keep after pruning
      INTEGER :: KEEP_TOL       ! prune when NUM_IN_CELL > KEEP_IN_CELL + KEEP_TOL

      INTEGER :: NUM_PRUNED     ! number of interior packets pruned

      INTEGER :: C, R, L        ! column, row and layer indices

      TYPE(PacketListMember), POINTER :: ITR, PLM
      TYPE(Packet), POINTER :: CLOSEST

      REAL(8), POINTER, SAVE :: KEEP_VAL(:) ! best values of selection variable
      TYPE(PacketListMemberPtr), POINTER, SAVE :: KEEP_LIST(:) ! PacketsInCell to keep
      REAL(8) :: XC, YC         ! coordinates of the center of the cell (m)
      REAL(8) :: V              ! value of selection variable for a packet
      REAL(8) :: V_WORST        ! worst value of selection variable on KEEP_LIST
      INTEGER :: K, K_WORST     ! indices in KEEP_* arrays

C-----------------------------------

      IF ( PRUNING_METHOD == NO_PRUNING ) RETURN

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

         PRUNING_STEP = 0

         KEEP_IN_CELL = MAX(HR_KEEP_IN_CELL, NR_KEEP_IN_CELL)
         ALLOCATE(
     &        KEEP_VAL(1:KEEP_IN_CELL),
     &        KEEP_LIST(1:KEEP_IN_CELL),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT /= 0 ) THEN
            XMSG = 'Failure allocating arrays for KEEP_CLOSEST method'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF
      END IF                    ! FIRSTIME

      PRUNING_STEP = PRUNING_STEP + 1
      IF ( PRUNING_STEP == PRUNING_FREQ ) THEN
         PRUNING_STEP = 0
      ELSE
         RETURN
      END IF

      NUM_PRUNED = 0

C Prune cells with too many packets
      SELECT CASE ( PRUNING_METHOD )

      CASE ( KEEP_CLOSEST )

         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( L <= HR_TOP .AND.
     &                 HR_BEGROW <= R .AND. R <= HR_ENDROW .AND.
     &                 HR_BEGCOL <= C .AND. C <= HR_ENDCOL ) THEN
                     KEEP_IN_CELL = HR_KEEP_IN_CELL
                     KEEP_TOL = HR_KEEP_TOL
                  ELSE
                     KEEP_IN_CELL = NR_KEEP_IN_CELL
                     KEEP_TOL = NR_KEEP_TOL
                  END IF

                  IF ( PacketsInCell(C, R, L)%Count >
     &                 KEEP_IN_CELL + KEEP_TOL ) THEN
                     XC = XORIG_GD + (C - 1) * XCELL_GD
                     YC = YORIG_GD + (R - 1) * YCELL_GD

C Initialize KEEP_VAL and KEEP_LIST with the first KEEP_IN_CELL PLM's on the
C PacketsInCell list
                     ITR => PacketsInCell(C, R, L)%Head
                     DO K = 1, KEEP_IN_CELL
                        IF ( .NOT. ASSOCIATED(ITR) ) THEN
                           XMSG = 'Bad pointer to PacketListMember'
                           CALL M3EXIT(PNAME, JDATE,JTIME, XMSG, XSTAT2)
                        END IF

                        KEEP_VAL(K) = (ITR%P%X-XC)**2 +(ITR%P%Y - YC)**2
                        KEEP_LIST(K)%PLM => ITR
                        ITR => ITR%Next
                     END DO

                     K_WORST = 0
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT

C If necessary (flagged by K_WORST == 0), update V_WORST and K_WORST
                        IF ( K_WORST == 0 ) THEN
                           V_WORST = -1.0E30
                           DO K = 1, KEEP_IN_CELL
                              IF ( KEEP_VAL(K) > V_WORST ) THEN
                                 V_WORST = KEEP_VAL(K)
                                 K_WORST = K
                              END IF
                           END DO
                        END IF

C Decide whether to delete ITR%P or KEEP_LIST(K_WORST)%PLM%P
                        V = (ITR%P%X - XC)**2 + (ITR%P%Y - YC)**2

                        IF ( V < V_WORST ) THEN
                           PLM => KEEP_LIST(K_WORST)%PLM
                           KEEP_VAL(K_WORST) = V
                           KEEP_LIST(K_WORST)%PLM => ITR
                           K_WORST = 0
                        ELSE
                           PLM => ITR
                        END IF
                        ITR => ITR%Next

C Delete the selected packet
                        CALL DELETE_PACKET(JDATE, JTIME, PLM, C, R, L)
                        NUM_PRUNED = NUM_PRUNED + 1
                     END DO
                  END IF
               END DO
            END DO
         END DO

      CASE ( KEEP_OLDEST )

         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( L <= HR_TOP .AND.
     &                 HR_BEGROW <= R .AND. R <= HR_ENDROW .AND.
     &                 HR_BEGCOL <= C .AND. C <= HR_ENDCOL ) THEN
                     KEEP_IN_CELL = HR_KEEP_IN_CELL
                     KEEP_TOL = HR_KEEP_TOL
                  ELSE
                     KEEP_IN_CELL = NR_KEEP_IN_CELL
                     KEEP_TOL = NR_KEEP_TOL
                  END IF

                  IF ( PacketsInCell(C, R, L)%Count >
     &                 KEEP_IN_CELL + KEEP_TOL ) THEN

C Initialize KEEP_VAL and KEEP_LIST with the first KEEP_IN_CELL PLM's on the
C PacketsInCell list
                     ITR => PacketsInCell(C, R, L)%Head
                     DO K = 1, KEEP_IN_CELL
                        IF ( .NOT. ASSOCIATED(ITR) ) THEN
                           XMSG = 'Bad pointer to PacketListMember'
                           CALL M3EXIT(PNAME, JDATE, JTIME, XMSG,XSTAT2)
                        END IF

                        KEEP_VAL(K) = ITR%P%AGE
                        KEEP_LIST(K)%PLM => ITR
                        ITR => ITR%Next
                     END DO

                     K_WORST = 0
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT

C If necessary (flagged by K_WORST == 0), update V_WORST and K_WORST
                        IF ( K_WORST == 0 ) THEN
                           V_WORST = 1.0E30
                           DO K = 1, KEEP_IN_CELL
                              IF ( KEEP_VAL(K) < V_WORST ) THEN
                                 V_WORST = KEEP_VAL(K)
                                 K_WORST = K
                              END IF
                           END DO
                        END IF

C Decide whether to delete ITR%P or KEEP_LIST(K_WORST)%PLM%P
                        V = ITR%P%AGE

                        IF ( V > V_WORST ) THEN
                           PLM => KEEP_LIST(K_WORST)%PLM
                           KEEP_VAL(K_WORST) = V
                           KEEP_LIST(K_WORST)%PLM => ITR
                           K_WORST = 0
                        ELSE
                           PLM => ITR
                        END IF
                        ITR => ITR%Next

C Delete the selected packet
                        CALL DELETE_PACKET(JDATE, JTIME, PLM, C, R, L)
                        NUM_PRUNED = NUM_PRUNED + 1
                     END DO
                  END IF
               END DO
            END DO
         END DO

c$$$ TODO:         CASE ( KEEP_MAX_MIN_AVG )

c$$$ TODO:         CASE ( PRUNE_SIMILAR )

      CASE DEFAULT
         XMSG =  '"' // TRIM(PRUNING_METHOD) // '"' //
     &        ' is not a valid selection'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)

      END SELECT

      WRITE(LOGDEV, '(5X, A, I6)')
     &     'Number of interior packets pruned: ', NUM_PRUNED

      END SUBROUTINE PRUNE_CELLS

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Open packet information output file and write initial data for all cells
C (including boundary cells)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE OPEN_PACKET_INFO(JDATE, JTIME, TSTEP)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
C      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
C      INCLUDE SUBST_IOFDESC     ! file header data structure
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
C      INCLUDE SUBST_COORD_ID    ! coord. and domain definitions (req IOPARMS)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! starting date (YYYYDDD)
      INTEGER, INTENT(IN) :: JTIME ! starting time (HHMMSS)
      INTEGER, INTENT(IN) :: TSTEP ! output timestep (HHMMSS)

C Local Variables:
      INTEGER :: LOGDEV

      CHARACTER(16) :: PNAME = 'OPEN_PACKET_INFO'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(16) :: VNAME

      INTEGER :: COUNT(1:NCOLS, 1:NROWS, 1:NLAYS)

      INTEGER :: C, R, L        ! column, row and layer indices
      INTEGER :: K              ! loop counter

      REAL :: AVG_AGE(1:NCOLS, 1:NROWS, 1:NLAYS) ! average age of packets in cell (hr)
      REAL :: MAX_AGE(1:NCOLS, 1:NROWS, 1:NLAYS) ! age of oldest packet in cell (hr)

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

C Set output file characteristics based on COORD.EXT and open it
      FTYPE3D = GRDDED3
      SDATE3D = JDATE
      STIME3D = JTIME
      TSTEP3D = TSTEP
      NVARS3D = 4
      NCOLS3D = GL_NCOLS
      NROWS3D = GL_NROWS
      NLAYS3D = NLAYS
      NTHIK3D = 1
      GDTYP3D = GDTYP_GD
      P_ALP3D = P_ALP_GD
      P_BET3D = P_BET_GD
      P_GAM3D = P_GAM_GD
      XORIG3D = XORIG_GD
      YORIG3D = YORIG_GD
      XCENT3D = XCENT_GD
      YCENT3D = YCENT_GD
      XCELL3D = XCELL_GD
      YCELL3D = YCELL_GD
      VGTYP3D = VGTYP_GD
      VGTOP3D = VGTOP_GD
      DO L = 1, NLAYS3D + 1
         VGLVS3D(L) = VGLVS_GD(L)
      END DO
      GDNAM3D = GRID_NAME       ! from HGRD_DEFN

      FDESC3D(1) = 'TG packet information output file'
      DO K = 2, MXDESC3
         FDESC3D(K) = ' '
      END DO

      VTYPE3D(1) = M3INT
      VNAME3D(1) = 'COUNT'
      UNITS3D(1) = ' '
      VDESC3D(1) = 'Number of packets in cell'

      VTYPE3D(2) = M3INT
      VNAME3D(2) = 'NEW_PACKETS'
      UNITS3D(2) = ' '
      VDESC3D(2) = 'Number of packets added to cell during output step'

      VTYPE3D(3) = M3REAL
      VNAME3D(3) = 'AVG_AGE'
      UNITS3D(3) = 'hr'
      VDESC3D(3) = 'Average age of packets in cell'

      VTYPE3D(4) = M3REAL
      VNAME3D(4) = 'MAX_AGE'
      UNITS3D(4) = 'hr'
      VDESC3D(4) = 'Age of oldest packet in cell'

C Create header
      IF ( MYPE .EQ. 0 ) THEN   ! open new
         IF ( .NOT. OPEN3(TG_PACKET_1, FSNEW3, PNAME)) THEN
            XMSG = 'Could not open ' // TRIM(TG_PACKET_1) // ' file'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF
      END IF

C Write the initial output as step 0 on the file (assumed to be in correct
C output units)
      DO L = 1, NLAYS
         DO R = 1, NROWS
            DO C = 1, NCOLS
               COUNT(C, R, L) = PacketsInCell(C, R, L)%Count
            END DO
         END DO
      END DO

      AVG_AGE = 0.0
      MAX_AGE = 0.0

      VNAME = 'COUNT'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, COUNT) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      VNAME = 'NEW_PACKETS'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, NEW_PACKETS) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      VNAME = 'AVG_AGE'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, AVG_AGE) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      VNAME = 'MAX_AGE'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, MAX_AGE) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      WRITE(LOGDEV, '(/5X, 3(A, 1X), I8, ":", I6.6)')
     &     'Timestep written to', TRIM(TG_PACKET_1),
     &     'for date and time', JDATE, JTIME
      WRITE(LOGDEV, '(5X, A, 1X, I8, ":", I6.6/)')
     &     'from timestep on initial data files for date and time',
     &     JDATE, JTIME

      END SUBROUTINE OPEN_PACKET_INFO

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Write to packet information output file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE WRITE_PACKET_INFO(JDATE, JTIME, NUM_STEPS)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
C      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
C      INCLUDE SUBST_COORD_ID    ! coord. and domain definitions (req IOPARMS)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! starting date (YYYYDDD)
      INTEGER, INTENT(IN) :: JTIME ! starting time (HHMMSS)
      INTEGER, INTENT(IN) :: NUM_STEPS ! number of steps to average

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

      CHARACTER(16) :: PNAME = 'WRITE_PACKET_INFO'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(16) :: VNAME

      INTEGER :: COUNT(1:NCOLS, 1:NROWS, 1:NLAYS)

      INTEGER :: C, R, L        ! column, row and layer indices

      TYPE(PacketListMember), POINTER :: ITR
      REAL :: AVG_AGE(1:NCOLS, 1:NROWS, 1:NLAYS) ! average age of packets in cell (hr)
      REAL :: MAX_AGE(1:NCOLS, 1:NROWS, 1:NLAYS) ! age of oldest packet in cell (hr)

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

         IF ( .NOT. OPEN3(TG_PACKET_1, FSRDWR3, PNAME) ) THEN
            XMSG = 'Could not open ' // TRIM(TG_PACKET_1) // ' file for update'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF
      END IF                    ! FIRSTIME

      AVG_AGE = 0.0
      MAX_AGE = 0.0
      DO L = 1, NLAYS
         DO R = 1, NROWS
            DO C = 1, NCOLS
               COUNT(C, R, L) = PacketsInCell(C, R, L)%Count
               IF ( COUNT(C, R, L) > 0 ) THEN
                  ITR => PacketsInCell(C, R, L)%Head
                  DO
                     IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                     AVG_AGE(C, R, L) = AVG_AGE(C, R, L) + ITR%P%AGE
                     IF ( ITR%P%AGE > MAX_AGE(C, R, L) )
     &                    MAX_AGE(C, R, L) = ITR%P%AGE
                     ITR => ITR%Next
                  END DO
                  AVG_AGE(C, R, L) = AVG_AGE(C, R, L) / COUNT(C, R, L)
               END IF
            END DO
         END DO
         AVG_LAY(L) = REAL(SUM(COUNT(:,:,L))) / REAL(NCOLS*NROWS)
         MAX_LAY(L) = MAXVAL(COUNT(:,:,L))
      END DO

      AVG_AGE = AVG_AGE / 3600.0 ! convert seconds to hours
      MAX_AGE = MAX_AGE / 3600.0

      MAX_CELL = MAXVAL(MAX_LAY)
      MAX_LOC = MAXLOC(COUNT)

      VNAME = 'COUNT'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, COUNT) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      VNAME = 'NEW_PACKETS'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, NEW_PACKETS) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF
      NEW_PACKETS = 0

      VNAME = 'AVG_AGE'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, AVG_AGE) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      VNAME = 'MAX_AGE'
      IF ( .NOT. WRITE3(TG_PACKET_1, VNAME, JDATE, JTIME, MAX_AGE) ) THEN
         XMSG = 'Could not write ' // TRIM(VNAME) //
     &        ' to ' // TRIM(TG_PACKET_1)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      WRITE(LOGDEV, '(/5X, 3(A, 1X), I8, ":", I6.6)')
     &     'Timestep written to', TRIM(TG_PACKET_1),
     &     'for date and time', JDATE, JTIME

      WRITE(LOGDEV, '(/)')
      DO L = 1, NLAYS
         WRITE(LOGDEV, '(5X, A, I3, 1X, F6.2)') 'Avg in cell at layer: ', L, AVG_LAY(L)
      END DO
c$$$      WRITE(LOGDEV, '(/)')
c$$$      DO L = 1, NLAYS
c$$$         WRITE(LOGDEV, '(5X, A, I3, 1X, I5)') 'Max in cell at layer: ', L, MAX_LAY(L)
c$$$      END DO
c$$$      WRITE(LOGDEV, '(5X, A, I4)') 'Max in cell: ', MAX_CELL
c$$$      WRITE(LOGDEV, '(5X, A, 3I3)') 'Cell with max: ', MAX_LOC

      WRITE(LOGDEV, '(/5X, A, F6.2)') 'Avg in cell over output step: ',
     &     SUM_AVG_CELL / REAL(NUM_STEPS)
      WRITE(LOGDEV, '(5X, A, I8)') 'Num cells: ', NUM_CELLS
      WRITE(LOGDEV, '(5X, A, I8/)') 'Tot used: ', LastID_GD

      SUM_AVG_CELL = 0.0

      END SUBROUTINE WRITE_PACKET_INFO

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Update a gridded mixing ratio array
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE UPDATE_MIX_ARRAY(MIX_NAME, JDATE, JTIME)

      USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
C      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
C      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
C      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
      INCLUDE SUBST_CONST
	INCLUDE SUBST_FILES_ID

C Arguments:
      CHARACTER(16), INTENT(IN) :: MIX_NAME
      INTEGER, INTENT(IN)       :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN)       :: JTIME ! current model time, coded HHMMSS

C Local Variables:
      CHARACTER(16) :: PNAME = 'UPDATE_MIX_ARRAY'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: C, R, L        ! column, row and layer indices

      TYPE(PacketListMember), POINTER :: ITR
      TYPE(Packet), POINTER :: CLOSEST, OLDEST

C-----------------------------------

C Select the appropriate gridded mixing ratio array and update method
      SELECT CASE ( MIX_NAME )

      CASE ( CLS_MIX_1 )
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L)
                  CLS_MIX(C, R, L, :) = CLOSEST%Mix
               END DO
            END DO
         END DO

      CASE ( AVG_MIX_1 )
         AVG_MIX = 0.0
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                     ITR => PacketsInCell(C, R, L)%Head
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                        AVG_MIX(C, R, L,:) = AVG_MIX(C,R,L,:)+ITR%P%Mix
                        ITR => ITR%Next
                     END DO
                     AVG_MIX(C, R, L, :) =
     &                    AVG_MIX(C,R,L,:)/PacketsInCell(C, R, L)%Count
                  ELSE
                     CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L)
                     AVG_MIX(C, R, L, :) = CLOSEST%Mix
                  END IF
               END DO
            END DO
         END DO

      CASE ( MIN_MIX_1 )
         MIN_MIX = 1.0E30
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                     ITR => PacketsInCell(C, R, L)%Head
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                        MIN_MIX(C,R,L,:)=MIN(ITR%P%Mix,MIN_MIX(C,R,L,:))
                        ITR => ITR%Next
                     END DO
                  ELSE
                     CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L)
                     MIN_MIX(C, R, L, :) = CLOSEST%Mix
                  END IF
               END DO
            END DO
         END DO

      CASE ( MAX_MIX_1 )
         MAX_MIX = -1.0E30
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                     ITR => PacketsInCell(C, R, L)%Head
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                        MAX_MIX(C,R,L,:)=MAX(ITR%P%Mix,MAX_MIX(C,R,L,:))
                        ITR => ITR%Next
                     END DO
                  ELSE
                     CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L)
                     MAX_MIX(C, R, L, :) = CLOSEST%Mix
                  END IF
               END DO
            END DO
         END DO

      CASE ( OLD_MIX_1 )
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                     ITR => PacketsInCell(C, R, L)%Head
                     OLDEST => ITR%P
                     IF ( PacketsInCell(C, R, L)%Count > 1 ) THEN
                        ITR => ITR%Next
                        DO
                           IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                           IF ( ITR%P%AGE > OLDEST%AGE ) OLDEST => ITR%P
                           ITR => ITR%Next
                        END DO
                     END IF
                     OLD_MIX(C, R, L, :) = OLDEST%Mix
                  ELSE
                     CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L)
                     OLD_MIX(C, R, L, :) = CLOSEST%Mix
                  END IF
               END DO
            END DO
         END DO

      CASE DEFAULT
         XMSG =  '"' // MIX_NAME // '"' //
     &        ' is not a valid selection'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)

      END SELECT

      END SUBROUTINE UPDATE_MIX_ARRAY

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END MODULE PACKETS
