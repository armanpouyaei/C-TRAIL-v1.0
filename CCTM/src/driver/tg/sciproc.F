
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
C $Header: /project/yoj/arc/CCTM/src/driver/yamo/sciproc.F,v 1.5 2011/10/21 16:10:43 yoj Exp $ 

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCIPROC ( CGRID, JDATE, JTIME, TSTEP, ASTEP )

C-----------------------------------------------------------------------
C Function:
C    Controls all of the physical and chemical processes for a grid
C    Operator splitting symmetric around chemistry
 
C Preconditions:
C    Dates and times represented YYYYDDD:HHMMSS.
C    No "skipped" dates and times.  All boundary input variables (layered or
C    non-layered) have the same perimeter structure with a thickness of NTHIK
 
C Subroutines and functions called:
C    All physical and chemical subroutines, 
C    DECOUPLE, COUPLE
 
C Revision History:
C    Oct. 24, 1995 by M. Talat Odman and Clint L. Ingram at NCSC: created
C    13 Dec 97 - Jeff - uncouple diffusion processes
C    27 Jun 98 - Jeff - sync step = chem step
C     7 Jul 01 - Shawn - mv cloud processing before chem
C       Jan 02 - Jeff - dyn alloc; remove PCGRID argument to ping
C    23 Jun 03 J.Young: for layer dependent advection tstep
C    18 Aug 03 J. Pleim - move vdiff before advection
C    29 Aug 03 J.Young: eliminate symmetric processing option and unused
C                       SUBST_GRID_ID string
C    30 May 05 J.Young: mass-conserving advection (yamo)
C     7 Jul 07 J.Young: No cksummer call after HADV - delay possible truncation of
C                       HADV negative conc's until after ZADV
C    25 Jul 08 S.Roselle: removed PING from modeled processes
C    16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
C    11 May 11 D.Wong: incorporated twoway model implementation
C    25 Mar 16 R.Mathur: incorporated pot. vorticity strat. O3 scaling
C-----------------------------------------------------------------------

      USE UTILIO_DEFN
	  USE GRID_CONF
	  USE PACKETS
	  USE CGRID_SPCS

      IMPLICIT NONE   

C Include files:

!      INCLUDE SUBST_PACTL_ID    ! PA control parameters

C Arguments:
      REAL, POINTER :: CGRID( :,:,:,: )

      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 3 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
                                ! TSTEP(3) = twoway model time step w.r.t. wrf time
                                !            step and wrf/cmaq call frequency

!     INTEGER      ASTEP( NLAYS )  ! layer advection time step
      INTEGER      ASTEP( : )      ! layer advection time step
	  REAL		AVG_MAX

C Parameters:

C External Functions:

C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'SCIPROC'

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 120 ) :: XMSG = ' '
      CHARACTER( 36 ) :: NMSG = 'After NEXTIME: returned JDATE, JTIME'

      INTEGER, SAVE :: LOGDEV

      INTEGER      SDATE        ! current science process date, coded YYYYDDD
      INTEGER      STIME        ! current science process time, coded HHMMSS
 
      INTEGER      ALLOCSTAT
      INTEGER      STATUS          ! ENV... status
      CHARACTER( 80 ) :: VARDESC   ! environment variable description
      CHARACTER( 80 ) :: MSG = ' '
      CHARACTER( 16 ) :: CTM_CKSUM = 'CTM_CKSUM'     ! env var for cksum on
      LOGICAL, SAVE   :: CKSUM     ! flag for cksum on, default = [T]

      INTERFACE
         SUBROUTINE CKSUMMER ( PNAME, CGRID, JDATE, JTIME )
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
         END SUBROUTINE CKSUMMER
         SUBROUTINE PA_UPDATE ( PNAME, CGRID, JDATE, JTIME, TSTEP )
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE PA_UPDATE
         SUBROUTINE HADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
            INTEGER, INTENT( IN )     :: ASTEP( : )
         END SUBROUTINE HADV
         SUBROUTINE ZADV ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE ZADV
         SUBROUTINE COUPLE ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE COUPLE
         SUBROUTINE DECOUPLE ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE DECOUPLE
         SUBROUTINE HDIFF ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE HDIFF
         SUBROUTINE VDIFF ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE VDIFF
         SUBROUTINE CLDPROC ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE CLDPROC
         SUBROUTINE CHEM ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE CHEM
#ifdef potvorto3
         SUBROUTINE PVO3 ( CGRID, JDATE, JTIME )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
         END SUBROUTINE PVO3
#endif
         SUBROUTINE AERO ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE AERO
      END INTERFACE
  
C-----------------------------------------------------------------------

		CGRID => AVG_MIX

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

         CKSUM = .TRUE.         ! default
         VARDESC = 'Cksum on flag'
         CKSUM = ENVYN( CTM_CKSUM, VARDESC, CKSUM, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            MSG = 'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            MSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) MSG, JTIME
         END IF

      END IF       ! if firstime

C Synchronization Step Physical Processes Sequence

C Do 3-D advection for the TG method - ADVSTEP has been modified so that
C TSTEP(2) satisfies the Courant condition for all layers

      CALL UPDATE_MIX_ARRAY(AVG_MIX_1, JDATE, JTIME) ! for aerosols in RDDEPV

	  CALL VDIFF ( CGRID, JDATE, JTIME, TSTEP )
		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'VERTICAL DIFF MAX', AVG_MAX	  
	  
!      IF ( CKSUM ) CALL CKSUMMER ( 'VDIFF', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'VDIF', CGRID, JDATE, JTIME, TSTEP )

  
C couple CGRID for advection and horizontal diffusion
      CALL ADVEC_TG(JDATE, JTIME, TSTEP(2))
		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'ADVECCC MAX', AVG_MAX	
		
!      CALL COUPLE_WRF( AVG_MIX, JDATE, JTIME, TSTEP )
!      IF ( CKSUM ) CALL CKSUMMER ( 'COUPLE', CGRID, JDATE, JTIME )

!      CALL HADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
!     IF ( CKSUM ) CALL CKSUMMER ( 'HADV', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'HADV', CGRID, JDATE, JTIME, TSTEP )

!      CALL ZADV ( CGRID, JDATE, JTIME, TSTEP )
!      IF ( CKSUM ) CALL CKSUMMER ( 'ADV', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'ZADV', CGRID, JDATE, JTIME, TSTEP )

	  CALL HDIFF ( CGRID, JDATE, JTIME, TSTEP )
		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'HORIZ DIFF MAX', AVG_MAX	
!      IF ( CKSUM ) CALL CKSUMMER ( 'HDIFF', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'HDIF', CGRID, JDATE, JTIME, TSTEP )

C decouple CGRID for cloud and chemistry

!      SDATE = JDATE
!      STIME = JTIME
!      CALL NEXTIME ( SDATE, STIME, TSTEP( 2 ) )

!      CALL DECOUPLE_WRF ( AVG_MIX, SDATE, STIME, TSTEP )
!      IF ( CKSUM ) CALL CKSUMMER ( 'DECOUPLE', CGRID, JDATE, JTIME )

!     CALL CLDPROC ( CGRID, JDATE, JTIME, TSTEP )
!      IF ( CKSUM ) CALL CKSUMMER ( 'CLDPROC', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'CLDS', CGRID, JDATE, JTIME, TSTEP )

!      CALL CHEM ( CGRID, JDATE, JTIME, TSTEP )
!      IF ( CKSUM ) CALL CKSUMMER ( 'CHEM', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'CHEM', CGRID, JDATE, JTIME, TSTEP ) 

#ifdef potvorto3
!      CALL PVO3 ( CGRID, JDATE, JTIME )
!      IF ( CKSUM ) CALL CKSUMMER ( 'PVO3', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'PVO3', CGRID, JDATE, JTIME, TSTEP )
#endif

!      CALL AERO ( CGRID, JDATE, JTIME, TSTEP )
!      IF ( CKSUM ) CALL CKSUMMER ( 'AERO', CGRID, JDATE, JTIME )
!      IF ( LIPR ) CALL PA_UPDATE ( 'AERO', CGRID, JDATE, JTIME, TSTEP )

      CALL NEXTIME ( JDATE, JTIME, TSTEP( 2 ) )
      WRITE( LOGDEV,'(/ 5X, A, I8, I7.6)' ) NMSG, JDATE, JTIME

      RETURN
      END
