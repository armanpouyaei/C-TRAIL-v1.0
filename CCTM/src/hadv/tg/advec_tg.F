
C***********************************************************************
C   Portions of Models-3/CMAQ software were developed or based on      *
C   information from various groups: Federal Government employees,     *
C   contractors working on a United States Government contract, and    *
C   non-Federal sources (including research institutions).  These      *
C   research institutions have given the Government permission to      *
C   use, prepare derivative works, and distribute copies of their      *
C   work in Models-3/CMAQ to the public and to permit others to do     *
C   so.  EPA therefore grants similar permissions for use of the       *
C   Models-3/CMAQ software, but users are requested to provide copies  *
C   of derivative works to the Government without restrictions as to   *
C   use by others.  Users are responsible for acquiring their own      *
C   copies of commercial software associated with Models-3/CMAQ and    *
C   for complying with vendor requirements.  Software copyrights by    *
C   the MCNC Environmental Modeling Center are used with their         *
C   permissions subject to the above restrictions.                     *
C***********************************************************************

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ADVEC_TG(JDATE, JTIME, TSTEP)

C-----------------------------------------------------------------------
C Function:
C   3-D advection for the Trajectory Grid (TG) method.

C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.

C Revision History:
C   07/26/05 Peter Percell at UH-IMAQS: created for TG method
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal domain specifications
      USE PACKETS               ! data and functions for TG packets
	  USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
!      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!      INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters (NLAYS)
      INCLUDE SUBST_CONST      ! constants
      INCLUDE SUBST_FILES_ID   ! file name parameters

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP ! advection time step (HHMMSS)

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

      INTEGER :: C, R, L        ! column, row and layer indices

      REAL :: AVG_CELL          ! average number of packets in a cell

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

         SUM_AVG_CELL = 0.0
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

      WRITE(LOGDEV, '(5X, A, F6.2/)') 'Avg in cell before: ', AVG_CELL		  
	  
C For each packet (including boundary packets), update (U, V, W), (X, Y, Z),
C (C, R, L) and PacketsInCell; delete packets that leave the interior of the
C domain; prune cells with too many packets; fill empty cells
      CALL UPDATE_POSITION(JDATE, JTIME, TSTEP)

      WRITE(LOGDEV, '(/5X, A, 1X, I8, ":", I6.6)') 'At date and time', JDATE, JTIME

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

      WRITE(LOGDEV, '(5X, A, F6.2/)') 'Avg in cell after update: ', AVG_CELL	  

      CALL PRUNE_CELLS(JDATE, JTIME)
	  
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

      WRITE(LOGDEV, '(5X, A, F6.2/)') 'Avg in cell after prune: ', AVG_CELL	  

      CALL FILL_EMPTY_CELLS(JDATE, JTIME)

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

      WRITE(LOGDEV, '(5X, A, F6.2/)') 'Avg in cell after fill: ', AVG_CELL

      RETURN
      END
