
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/yoj/arc/CCTM/src/driver/yamo/wr_conc.F,v 1.6 2011/10/21 16:10:43 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE WR_CONC2 ( MIX_FILE_1, JDATE, JTIME, TSTEP )

C Revision History:
C   10/13/99 David Wong at LM
C      -- Called from driver, where CGRID is a pointer (subset) of PCGRID.
C         Necessary, to keep from referencing parts of PCGRID that don't
C         belong to CGRID.
C    1/31/2000 Jeff Young
C      -- f90 memory mgmt
C   Jeff - Dec 00 - move CGRID_MAP into f90 module
C   Jeff - Feb 01 - assumed shape arrays
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; DBUFF for WRITE3
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   30 May 05 J.Young: optional save derived vert. vel. to conc file
C   29 Aug 06 J.Young: in conjunction with state CGRID file
C   30 Jan 10 David Wong: modified the code to write one species at a time
C                         in case ALLVAR3 cannot deal with large CONC size
C   14 Sep 10 J.Young: fix ALLVAR3 bug in writing W_VEL
C   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
C   11 May 11 D.Wong: - let non I/O PEs open CTM_CONC_1 for read only
C   12 Aug 15 D.Wong: - Replaced MYPE with IO_PE_INCLUSIVE for parallel
C                       I/O implementation
C                     - Added code to allow non IO processor to open
C                       CTM_CONC_1
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE STD_CONC              ! standard CONC
      USE WVEL_DEFN             ! derived vertical velocity component
      USE UTILIO_DEFN
      USE PACKETS
      USE CGRID_SPCS

      IMPLICIT NONE

C Include Files:

      INCLUDE SUBST_FILES_ID    ! I/O definitions and declarations

!     REAL      :: SGRID( :,:,:,: )
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP        ! output timestep (HHMMSS)
      
C Arguments:
      CHARACTER(16),    INTENT(IN) ::   MIX_FILE_1

C Local variables:

      REAL, POINTER      :: MGRID( :,:,:,: )     ! alias for a gridded mixing ratio array

      CHARACTER( 16 ) :: PNAME = 'WR_CONC2'
      CHARACTER( 96 ) :: XMSG = ' '

      REAL, ALLOCATABLE, SAVE :: JACF( :,:,: ) ! full-layer Jacobian
      REAL, ALLOCATABLE, SAVE :: MSFX2( :,: )  ! map scale factor ** 2
      REAL, ALLOCATABLE :: DBUFF ( :,:,: )     ! input for W_VEL
      REAL, ALLOCATABLE :: WBUFF ( :,:,: )     ! output for W_VEL
      INTEGER      ALLOCSTAT

      INTEGER, SAVE :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2

      INTEGER      C, R, K, L, V   ! loop induction variables
      INTEGER       SPC, VAR        ! species loop counters

      CHARACTER( 16 ), ALLOCATABLE, SAVE :: VNAME( : )
      INTEGER, SAVE :: NVARS

C-----------------------------------------------------------------------
      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

!         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C open conc file for update

         DO K = 1, NMIX_FILES
            IF ( .NOT. OPEN3( MIX_FILE(K), FSRDWR3, PNAME ) ) THEN
               XMSG = 'Could not open ' // MIX_FILE(K) // ' file for update'
            END IF
         END DO

      END IF

C Update the gridded mixing ratio array corresponding to the file and make
C MGRID an alias for the array
      CALL UPDATE_MIX_ARRAY(MIX_FILE_1, JDATE, JTIME)

      SELECT CASE ( MIX_FILE_1 )
         CASE ( CLS_MIX_1 )
            MGRID => CLS_MIX

         CASE ( AVG_MIX_1 )
            MGRID => AVG_MIX

         CASE ( MIN_MIX_1 )
            MGRID => MIN_MIX

         CASE ( MAX_MIX_1 )
            MGRID => MAX_MIX

         CASE ( OLD_MIX_1 )
            MGRID => OLD_MIX

         CASE DEFAULT
            XMSG =  '"' // MIX_FILE_1 // '"' //
     &           ' is not a valid file selection'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END SELECT

      ALLOCATE ( DBUFF( MY_NCOLS,MY_NROWS,NLAYS ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'Failure allocating DBUFF'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

      DO SPC = 1, N_GC_CONC
         VAR = GC_STRT - 1 + GC_CONC_MAP( SPC )

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,VAR )

         IF ( .NOT. WRITE3( MIX_FILE_1, GC_CONC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' // TRIM( GC_CONC( SPC ) ) // ' to '
     &           // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      DO SPC = 1, N_AE_CONC
         VAR = AE_STRT - 1 + AE_CONC_MAP( SPC )

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,VAR )

         IF ( .NOT. WRITE3( MIX_FILE_1, AE_CONC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' // TRIM( AE_CONC( SPC ) ) // ' to '
     &           // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      DO SPC = 1, N_NR_SPC
         VAR = NR_STRT - 1 + SPC

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,VAR )

         IF ( .NOT. WRITE3( MIX_FILE_1, NR_SPC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' // TRIM( NR_SPC( SPC ) ) // ' to '
     &           // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      DO SPC = 1, N_TR_SPC
         VAR = TR_STRT - 1 + SPC

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,VAR )

         IF ( .NOT. WRITE3( MIX_FILE_1, TR_SPC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' // TRIM( TR_SPC( SPC ) ) // ' to '
     &           // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      DEALLOCATE ( DBUFF )

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &     'Timestep written to', TRIM( MIX_FILE_1 ),
     &     'for date and time', JDATE, JTIME

      RETURN
      END
