
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

!C RCS file, release, date & time of last delta, author, state, [and
!locker]
!C $Header: /project/yoj/arc/CCTM/src/emis/emis/LTNG_DEFN.F,v 1.6
!2012/01/27 19:57:47 sjr Exp $

!C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE LTNG_DEFN

!C----------------------------------------------------------------------
!C Function: production of NO from lightning

!C Revision History:
!C     1 Jul 2010: Jeff Young
!C     1 Jan 2011: Rob Pinder: added support for online calculation of 
!C                 lightning NO from convective precip and NLDN data
!C    10 Mar 2011: Jeff Young: code revisions
!C    11 Apr 2011: S.Roselle: replaced I/O API include files with UTILIO_DEFN
!C    11 May 2011 D.Wong: incorporated twoway model implementation
!C     6 May 2013 D.Wong: replaced species name RC with RCA in the twoway model
!C    10 Jun 2013 D.Wong: modified the code to work in the twoway model mode
!C    24 Sep 2013 D.Wong: modified the code to allow finer met data time step
!C                        rather than fixed with one hour data
!C    15 Aug 2016 D.Wong: Replaced MYPE with IO_PE_INCLUSIVE for parallel I/O
!C                        implementation
!C    16 Aug 2016 D.Kang: Updating the Lightning NO calculation with NLDN hourly
!C                        data directly and using NLDN-RC log linear relationship,
!C                        completely changing the data input process   
!C     1 Feb 2017 D.Kang: Modify the parameter (linear-log linear regression)scheme to
!C                        set log-linear when RC is greater than the
!C                        exp(log-intercept), and move the LTNOx
!C                        production rate to environmental variable to be set in run script
!C    18 Apr 2017 D. Kang: Make this routine compatible with the two-way coupled model
!C    17 Jul 2017 D. Wong: In subroutine GET_LTNG: 
!C                         * implemented a bug fix for determining 
!C                           LASTTIC base on accumulation of time step, TSTEP(2)
!C                         * called SUBHFILE only once to improve code efficiency
!C    02 Aug 2017 D. Wong: In subroutine GET_LTNG: 
!C                         * removed block of code that computes CDATE and CTIME
!C                         * replaced CDATE and CTIME with JDATE and JTIME, respectively 
!C                           when calling INTERPX
!C----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL,              SAVE :: LTNG_NO          ! flag for lightning NO emissions
      INTEGER,              SAVE :: LTNG_MAP         ! map to GC NO
      REAL,    ALLOCATABLE, SAVE :: VDEMIS_LT( :,:,: )   ! lightning emis

      PUBLIC LTNG_NO, LTNG_MAP, VDEMIS_LT, LTNG_INIT, GET_LTNG, LTSPC
      PRIVATE

!C lightning emis species name
      CHARACTER( 16 ), PARAMETER :: LTSPC = 'NO'

!C lightning log linear regression parameters with RC 
      REAL            :: SQUAREKM
      REAL, SAVE      :: SCL_FACTOR
      REAL, SAVE      :: MOLSNCG
      REAL, SAVE      :: MOLSNIC

      LOGICAL,              SAVE :: NLDNSTRIKE     ! flag to use NLDN STRIKE directly 
      LOGICAL,              SAVE :: LTNGDIAG   ! flag to turn on lightning NO diagnostics
      INTEGER,              SAVE :: LTLYRS     ! no. of emissions layers
      REAL,    ALLOCATABLE, SAVE :: VDEMIS_DIAG( :,:,: ) ! lightning NO diagnostic
      REAL,    ALLOCATABLE, SAVE :: COLUMN_DIAG( :,: ) ! column total NO

!C file name for off line lightning NO input; set to "InLine" for inline production
      CHARACTER( 16 ), SAVE :: LNAME
!C file name for inline lightning parameters, (only used if LNAME = "InLine")
      CHARACTER( 16 ), SAVE :: ANAME
      CHARACTER( 16 ), SAVE :: ONAME
!C file name for met file for inline lightning
      CHARACTER( 16 ), SAVE :: MNAME
      CHARACTER( 16 ), SAVE :: RNAME
      CHARACTER( 16 ), SAVE :: RC_NAME       ! RC name: old is RC and CCnew is RCA

!C allocate these if LNAME = 'InLine'
      REAL,    ALLOCATABLE, SAVE :: LTNG_PRSFC  ( :,: ) ! surface pressure
      REAL,    ALLOCATABLE, SAVE :: LTNG_RC     ( :,: ) ! convective rainfall
!C allocate these if LPARAM
      REAL,    ALLOCATABLE, SAVE :: NLDN_STRIKE ( :,: ) ! Hourly NLDN strike data
      REAL,    ALLOCATABLE, SAVE :: OCEAN_MASK  ( :,: ) ! reduce offshore strikes
      REAL,    ALLOCATABLE, SAVE :: SLOPE       ( :,: ) ! reduce offshore strikes
      REAL,    ALLOCATABLE, SAVE :: INTERCEPT   ( :,: ) ! reduce offshore strikes
      REAL,    ALLOCATABLE, SAVE :: SLOPE_lg    ( :,: ) ! reduce offshore strikes
      REAL,    ALLOCATABLE, SAVE :: INTERCEPT_lg( :,: ) ! reduce offshore strikes
      REAL,    ALLOCATABLE, SAVE :: ICCG        ( :,: ) ! intercloud strikes per cloud to ground strike
!C Vertical coord values
      REAL,    ALLOCATABLE, SAVE :: VGLVSLT( : )

!C    scenario time/date needed for diagnostic output
      INTEGER, SAVE    :: SDATE     ! scenario start date
      INTEGER, SAVE    :: STIME     ! scenario start time

      INTEGER, SAVE    :: NSTEPS    ! run duration
      INTEGER, SAVE    :: NTICS = 0 ! no. of substeps within an output tstep
      INTEGER, SAVE    :: MTICS = 0 ! temporary sub for NTICS before it set to 0 
      INTEGER, SAVE    :: LDATE     ! test date to update emissions diag avg
      INTEGER, PRIVATE :: LTIME     ! test time to update emissions diag avg
      INTEGER, PRIVATE :: FTIME     ! Start time to Write to Diagnostic files
      INTEGER, SAVE    :: LOGDEV


      INTEGER, SAVE    :: GXOFF, GYOFF   ! global origin offset from file
      INTEGER, SAVE    :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      INTEGER :: LIGHTNING_FILE_TSTEP, 
     &           LIGHTNING_FILE_SDATE,
     &           LIGHTNING_FILE_STIME

      CONTAINS

!C======================================================================
!C Initialize lightning routines

         FUNCTION LTNG_INIT ( JDATE, JTIME, TSTEP, EMLAYS ) RESULT ( SUCCESS )

         USE GRID_CONF ! horizontal & vertical domain specifications
         USE CGRID_SPCS         ! CGRID mechanism species
         USE UTILIO_DEFN

         IMPLICIT NONE

!C Includes:
         INCLUDE SUBST_CONST     ! constants
         INCLUDE SUBST_FILES_ID  ! file name parameters

!C Arguments:
C        INTEGER, INTENT( IN ) :: JDATE, JTIME, TSTEP
         INTEGER, INTENT( IN ) :: EMLAYS     ! number of emissions layers
         INTEGER :: JDATE, JTIME, TSTEP

         LOGICAL SUCCESS

!C External Functions:
         INTEGER, EXTERNAL :: SETUP_LOGDEV
         LOGICAL, EXTERNAL :: CHKGRID
!        CHARACTER( 16 )   :: iccgFile = 'ICCG' 
!        INTEGER :: ICCG_UNIT    ! scenario start time
         INTEGER :: COL, ROW     ! local col and row variable 
!C Local Variables:
         CHARACTER( 16 ), SAVE :: CTM_LTNG_NO = 'CTM_LTNG_NO'  ! env var Cflag
!C value to start log linear regression  
         CHARACTER( 16 )       :: PNAME = 'LTNG_INIT'
         CHARACTER( 80 )       :: VARDESC   ! env variable description
         CHARACTER( 120 )      :: XMSG = ' '
         REAL                  :: LAT, LON, X, Y, VAL
         INTEGER               :: C, R
         CHARACTER( 80 )       :: DUMP

         LOGICAL LTNGPARAM           ! env var to use lightning NO parameters file
         INTEGER LSPCS               ! no. of lightning species
         REAL    CONVEM              ! conversion for emissions rates to Kg/s
         INTEGER               :: I, J, K, L, V, STATUS
         INTEGER               :: CJDATE
         CHARACTER( 7 )        :: SJDATE

         LOGICAL OK
         INTEGER               :: SPC

!-----------------------------------------------------------------------

         SUCCESS = .TRUE.

         LOGDEV = SETUP_LOGDEV()

!C Lightning NO production
!C Set LTNG_NO to Y or T to turn on lightning NO production
         VARDESC = 'Flag for lightning emissions'
         LTNG_NO = ENVYN( CTM_LTNG_NO, VARDESC, .FALSE., STATUS )
         IF ( .NOT. LTNG_NO ) RETURN

!C Read in scenario time/date
         SDATE = ENVINT( 'CTM_STDATE', 'Scenario Start (YYYYJJJ)', 0, STATUS )
         STIME = ENVINT( 'CTM_STTIME', 'Scenario Start (HHMMSS)', 0, STATUS )
         NSTEPS = ENVINT( 'CTM_RUNLEN', 'Scenario Duration (HHMMSS)', 24, STATUS )
         NSTEPS = NSTEPS / TSTEP           ! initscen guarantees divisibility

         WRITE( SJDATE, '(I7)' ) JDATE
         READ( SJDATE, '(4X, I3)' ) CJDATE

!C Check grid definition (initialize if first call)
!C Set up domain based on met file
!C Open met file
#ifdef twoway
            IF ( IO_PE_INCLUSIVE ) THEN
               MNAME = PROMPTMFILE(
     &                 'Enter name for gridded met input file',
     &                 FSRDWR3, CTM_CONC_1, PNAME )
            ELSE
               MNAME = PROMPTMFILE(
     &                 'Enter name for gridded met input file',
     &                 FSREAD3, CTM_CONC_1, PNAME )
            END IF
#else
            MNAME = PROMPTMFILE(
     &              'Enter name for gridded met input file',
     &              FSREAD3, 'MET_CRO_3D', PNAME )
#endif

!C Read description of domain (there may be a better way to get this info)
            IF ( .NOT. DESC3( MNAME ) ) THEN
               XMSG = 'Could not get description of file "'
     &              // TRIM( MNAME ) // '"'
               CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF

!C Check grid definition (initialize if first call)
            OK = CHKGRID( MNAME )

!C Get domain window info for met_cro_3d file
            CALL SUBHFILE ( MNAME, GXOFF, GYOFF,
     &              STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         VARDESC = 'Flag for lightning NO production rate'
         MOLSNCG = ENVREAL( 'MOLSNCG', VARDESC, 350.0, STATUS )
         MOLSNIC = ENVREAL( 'MOLSNIC', VARDESC, 350.0, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Env. variable set, but empty: Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         END IF

         VARDESC = 'lightning NO file name or InLine flag'
         CALL ENVSTR( 'LTNGNO', VARDESC, 'LTNGNO', LNAME, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Env. variable set, but empty: Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         END IF

!C Is lightning NO production inline, or from a file?
         IF ( LNAME .EQ. "InLine" ) THEN  ! inline lightning NO production
            XMSG = 'Using in-line lightning NO production'
            CALL M3MSG2( XMSG )

!C Is there a file specifying hourly NLDN strikes?
            VARDESC = 'Use NLDN STRIKE'
            NLDNSTRIKE = ENVYN( 'USE_NLDN', VARDESC, .FALSE., STATUS )
            IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
            IF ( STATUS .EQ. 1 ) THEN
               XMSG = 'Environment variable improperly formatted'
               CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            ELSE IF ( STATUS .EQ. -1 ) THEN
               XMSG = 'Env. variable set, but empty: Using default:'
               WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
            ELSE IF ( STATUS .EQ. -2 ) THEN
               XMSG = 'Environment variable not set ... Using default:'
               WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
            END IF

!C Open lightning parameters file?
               ONAME = PROMPTMFILE(
     &                    'Enter name for lightning parameters file',
     &                    FSREAD3, 'LTNGPARMS_FILE', PNAME )
               RNAME = PROMPTMFILE(
     &                    'Enter name for gridded met input file',
     &                    FSREAD3, 'MET_CRO_2D', PNAME )
!C Read description of normalized emissions file
            IF(NLDNSTRIKE) THEN
                  XMSG = 'Using hourly NLDN Strike data'
                  CALL M3MSG2( XMSG )
                  ANAME = PROMPTMFILE(
     &                    'Enter name for NLDN STRIKE Parameters',
     &                    FSREAD3, 'NLDN_STRIKES', PNAME )

               IF ( .NOT. DESC3( ANAME ) ) THEN
                  XMSG = 'Could not get description of file "'
     &                 // TRIM( ANAME ) // '"'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
                  SUCCESS = .FALSE.; RETURN
               END IF

               LIGHTNING_FILE_TSTEP = TIME2SEC ( TSTEP3D )
               LIGHTNING_FILE_SDATE = SDATE3D
               LIGHTNING_FILE_STIME = STIME3D

            ELSE
                 XMSG = 'Using derived parameters for KF scheme '
                 CALL M3MSG2( XMSG )
            END IF  !END if NLDNSTRIKE

!C Set Lambert projection
            IF ( .NOT. SETLAM( SNGL( P_ALP3D ),
     &            SNGL( P_BET3D ),
     &            SNGL( P_GAM3D ),
     &            SNGL( XCENT3D ),
     &            SNGL( YCENT3D ) ) ) THEN
                  CALL M3EXIT( 'CELL_LL', 0, 0,
     &            'Lambert projection setup error', 2 )
            END IF

            SQUAREKM = REAL( ( XCELL3D * YCELL3D * 1.0D-6 ), 4 )
!C Set up vertical layers
            ALLOCATE( VGLVSLT( 0:EMLAYS ), STAT = STATUS )
            CALL CHECKMEM( STATUS, 'VGLVSLT', PNAME )
            SPC = INDEX1( 'RCA', NVARS3D, VNAME3D )
            IF ( SPC .GT. 0 ) THEN
               RC_NAME = 'RCA'
            ELSE
               RC_NAME = 'RC'
            END IF

            ALLOCATE( ICCG        ( NCOLS,NROWS ),
     &                   OCEAN_MASK  ( NCOLS,NROWS ),
     &                   SLOPE       ( NCOLS,NROWS ),
     &                   INTERCEPT   ( NCOLS,NROWS ),
     &                   SLOPE_lg    ( NCOLS,NROWS ),
     &                   INTERCEPT_lg( NCOLS,NROWS ),
     &                   LTNG_PRSFC  ( NCOLS,NROWS ),
     &                   NLDN_STRIKE ( NCOLS,NROWS ),
     &                   LTNG_RC     ( NCOLS,NROWS ),STAT = STATUS )
               IF ( STATUS .NE. 0 ) THEN
                  XMSG = 'ICCG, OCEAN_MASK, LTNG_Parameters'
     &                 // '  memory allocation failed'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
                  SUCCESS = .FALSE.; RETURN
               END IF
!C read in the ocean mask and parameter file if needed, if NLDN strike
!C is used, only the OCEAN MASK values are needed, the rest will be ignored          
               IF ( .NOT. XTRACT3( ONAME, "OCNMASK", 1, 1,  
     &                             STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                             JDATE, JTIME, OCEAN_MASK ) ) THEN
                  XMSG = 'Could not read OCNMASK from' // TRIM( ONAME ) // '"'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

               IF ( .NOT. XTRACT3( ONAME, "SLOPE", 1, 1,
     &                             STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3, 
     &                             JDATE, JTIME, SLOPE ) ) THEN
                  XMSG = 'Could not read SLOPE from' // TRIM( ONAME ) // '"'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               IF ( .NOT. XTRACT3( ONAME, "INTERCEPT", 1, 1, 
     &                             STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                             JDATE, JTIME, INTERCEPT ) ) THEN
                  XMSG = 'Could not read INTERCEPT from' // TRIM( ONAME) // '"'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               IF ( .NOT. XTRACT3( ONAME, "SLOPE_lg", 1, 1,
     &                             STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                             JDATE, JTIME, SLOPE_lg ) ) THEN
                  XMSG = 'Could not read SLOPE_lg from' // TRIM( ONAME ) // '"'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               IF ( .NOT. XTRACT3( ONAME, "INTERCEPT_lg", 1, 1,
     &                             STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                             JDATE, JTIME, INTERCEPT_lg ) ) THEN
                  XMSG = 'Could not read INTERCEPT_lg from' // TRIM( ONAME ) // '"'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               IF ( CJDATE .GT. 90 .AND. CJDATE .LT. 274 ) THEN
                  IF ( .NOT. XTRACT3( ONAME, "ICCG_SUM", 1, 1,
     &                                STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                             JDATE, JTIME, ICCG ) ) THEN
                     XMSG = 'Could not read ICCG from' // TRIM( ONAME ) // '"'
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                  END IF
               ELSE
                  IF ( .NOT. XTRACT3( ONAME, "ICCG_WIN", 1, 1,
     &                                STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                             JDATE, JTIME, ICCG ) ) THEN
                     XMSG = 'Could not read ICCG from' // TRIM( ONAME ) // '"'
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                  END IF
               END IF
!            END IF  !LPARAM
!C Store local layer information
            DO L = EMLAYS, 0, -1
                  VGLVSLT( L ) = VGLVS3D( L+1 )
                  WRITE( LOGDEV,'(5X, A, I3, A, F11.7)' ) 'VGLVSLT(', L, ' ):', VGLVSLT( L )
            END DO

!C Lightning Diagnostic file?
            LTNGDIAG = ENVYN( 'LTNGDIAG',
     &              'Write lightning NO production diagnostic file?',
     &              .FALSE., STATUS )

            LDATE = SDATE; LTIME = STIME
            IF ( LTNGDIAG ) THEN
!C Build description for, and open lightning diagnostic file
!C (all but variables-table and horizontal domain in description is borrowed from MNAME)
               FTIME = STIME
               CALL NEXTIME( LDATE, FTIME, TSTEP )
               SDATE3D = LDATE
               STIME3D = FTIME
               TSTEP3D = TSTEP

               FTYPE3D = GRDDED3
               NCOLS3D = GL_NCOLS
               NROWS3D = GL_NROWS
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
               GDNAM3D = GRID_NAME  ! from HGRD_DEFN

               NLAYS3D = EMLAYS
               DO L = 1, NLAYS3D + 1
                  VGLVS3D( L ) = VGLVS_GD( L )
               END DO

               NVARS3D = 1
               VNAME3D( 1 ) = LTSPC
               VDESC3D( 1 ) = 'hourly average NO produced from lightning'
               VTYPE3D( 1 ) = M3REAL
               UNITS3D( 1 ) = 'mol/s'

               FDESC3D = ' '   ! array assignment
               FDESC3D( 1 ) = 'Gridded lightning NO production from CMAQ'
               FDESC3D( 2 ) = '/from/ ' // PNAME
               FDESC3D( 3 ) = '/Version/ CMAQ'

!C Open output file (mol/s)
               IF ( IO_PE_INCLUSIVE ) THEN
                  IF ( .NOT. OPEN3( CTM_LTNGDIAG_1, FSUNKN3, PNAME ) ) THEN
                     XMSG = 'Could not open ' // TRIM( CTM_LTNGDIAG_1 )
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                  END IF
                  NLAYS3D = 1
                  VDESC3D( 1 ) = 'Column NO produced from lightning'
                  IF ( .NOT. OPEN3( CTM_LTNGDIAG_2, FSUNKN3, PNAME ) ) THEN
                     XMSG = 'Could not open ' // TRIM( CTM_LTNGDIAG_2 )
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                  END IF

               END IF

               ALLOCATE( VDEMIS_DIAG( NCOLS,NROWS,EMLAYS ), STAT = STATUS )
               IF ( STATUS .NE. 0 ) THEN
                  XMSG = 'VDEMIS_DIAG memory allocation failed'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
                  SUCCESS = .FALSE.; RETURN
               END IF
               VDEMIS_DIAG = 0.0   ! array assignment
               ALLOCATE( COLUMN_DIAG( NCOLS,NROWS), STAT = STATUS )
               IF ( STATUS .NE. 0 ) THEN
                  XMSG = 'COLUMN_DIAG memory allocation failed'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
                  SUCCESS = .FALSE.; RETURN
               END IF
               COLUMN_DIAG = 0.0   ! array assignment
            END IF   ! LTNGDIAG
         ELSE   ! lightning emissions off line

!C Lightning NO production from an input file
            CALL M3MSG2( 'Using lightning NO production from a file' )

!C Open lightning emissions input file
            LNAME = PROMPTMFILE(
     &              'Enter name for lightning emissions input file',
     &              FSREAD3, 'LTNGNO', PNAME )

!C Get description of emissions file
            IF ( .NOT. DESC3( LNAME ) ) THEN
               XMSG = 'Could not get description of file "'
     &              // TRIM( LNAME ) // '"'
               CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF

!C Check grid definition (intialize, if first call)
            OK = CHKGRID( LNAME )


         END IF   ! IF ( LNAME .EQ. "InLine" ) inline or offline lightning NO production

         LTLYRS = EMLAYS

!C Lightning to gas-phase species map

         LTNG_MAP = INDEX1( LTSPC, N_GC_EMIS, GC_EMIS )
         IF ( LTNG_MAP .EQ. 0 ) THEN
            XMSG = TRIM( LTSPC ) // ' not found in GC_EMIS table'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE( VDEMIS_LT( NCOLS,NROWS,LTLYRS ), STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = 'VDEMIS_LT memory allocation failed'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         END IF

         RETURN

         END FUNCTION LTNG_INIT

!C======================================================================
!C Get NO produced from lightning in VDEMIS_LT

         SUBROUTINE GET_LTNG ( JDATE, JTIME, TSTEP )

         USE GRID_CONF             ! horizontal & vertical domain specifications
         USE UTILIO_DEFN

         IMPLICIT NONE

!     need to change this to set this path in the Makefile
!        include "/usr/local/intel/mkl/include/mkl_vml.f77"
         INCLUDE SUBST_CONST     ! constants
         INCLUDE SUBST_FILES_ID  ! file name parameters

         INTEGER JDATE, JTIME, TSTEP( 3 )
         INTEGER STATUS

         REAL,   PARAMETER :: CONVPA = 1.0E-2  ! convert Pa to hPa
!        REAL,   PARAMETER :: WK = 8.0         ! shape parameter for weibull distribution
!        REAL,   PARAMETER :: WLAMBDA = 700.0  ! scale parameter for weibull distribution
         REAL,   PARAMETER :: WMU = 350.0      ! mean
         REAL,   PARAMETER :: WSIGMA = 200.0   ! standard deviation
         REAL,   PARAMETER :: W2MU = 600.0     ! mean
         REAL,   PARAMETER :: W2SIGMA = 50.0   ! standard deviation
         REAL,   PARAMETER :: SQRT2 = 1.414213562731

         INTEGER COL, ROW, LAY ! iterator variables

         REAL    PCALC    ! pressure level for NO vertical distribution (hPa)
         REAL    BOTTOM   ! pressure at bottom of grid cell (hPa)
         REAL    TOP      ! pressure at top of grid cell (hPa)        
         REAL    BOTTOM_FRAC, TOP_FRAC ! their difference is the fraction of lightning NO in this grid cell
         REAL    BOTTOM_FRAC2, TOP_FRAC2
         REAL    SUM_FRAC ! stores the sum of vertical fractions to re-normalize the column
         REAL    WEIGHT ! used to normalize emissions to total amount
         REAL    inErfB, inErfT !  nputs to error funciton calculation
         REAL    outErfB, outErfT ! outputs from error funciton calculation
         REAL :: LTEMIS( LTLYRS )
         REAL    XCELLR, YCELLR   ! cell spacing ratio to 36Km
         REAL    FLASH_FAC        ! lightning flashes factor

         LOGICAL, SAVE :: LASTTIC   ! true: last sync step this output tstep
         REAL          :: DIVFAC   ! averaging factor for diagnostic file

         CHARACTER( 16 )  :: MNAME
         CHARACTER( 16 )  :: PNAME = 'GET_LTNG'
         CHARACTER( 120 ) :: XMSG = ' '

         INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2

         INTEGER               :: CJDATE
         INTEGER               :: DTSTEP
         CHARACTER( 7 )        :: SJDATE

         LOGICAL, SAVE         :: R_READY = .TRUE.
         REAL,    ALLOCATABLE  :: COLUMN_LTNG_NO( :,: ) ! column total NO

         INTEGER, SAVE :: TOT_TSTEP
         LOGICAL :: FIRSTIME = .TRUE.

C statement function for ERF approximation
         REAL :: ERF                ! ERF approx. statement function
         REAL :: X                  ! dummy argument for ERF
         ERF( X ) = SIGN( 1.0, X ) * SQRT( 1.0 - EXP( -4.0 * X * X / PI ) )
!-----------------------------------------------------------------------

         IF ( LNAME .EQ. "InLine" ) THEN
!C case of inline lightning NO production
!C initialize output array
            VDEMIS_LT = 0.0

            WRITE( LOGDEV,'(5X, A)' ) 'InLine LNOx production'
            DTSTEP = TIME2SEC( TSTEP( 1 ) ) / TIME2SEC( TSTEP( 2 ) ) 
            !C Open me filet
!C Get domain window info for met_cro_2d file
            IF (FIRSTIME) THEN
               CALL SUBHFILE ( RNAME, GXOFF, GYOFF,
     &                         STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2 )
               TOT_TSTEP = 0
               FIRSTIME = .FALSE.
            END IF

            WRITE( SJDATE, '(I7)') JDATE
            READ( SJDATE, '(4X, I3)') CJDATE

!C read in the hourly lightning strike file
#ifdef twoway
            IF(.NOT. NLDNSTRIKE ) THEN
              IF ( .NOT. INTERPX( RNAME, RC_NAME, PNAME,
     &                       STRTCOLMC2,ENDCOLMC2,STRTROWMC2,ENDROWMC2, 1,1,
     &                       JDATE, JTIME, LTNG_RC) ) THEN
                     XMSG = 'Could not read LTNG_RC from ' // TRIM(RNAME)
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
              END IF
!C Interpret the timestep RC value into hourly value
              LTNG_RC = LTNG_RC*DTSTEP
            END IF
            IF ( .NOT. INTERPX( RNAME, 'PRSFC', PNAME,
     &                       STRTCOLMC2,ENDCOLMC2,STRTROWMC2,ENDROWMC2, 1,1, 
     &                       JDATE, JTIME, LTNG_PRSFC ) ) THEN
                     XMSG = 'TwoWay: Could not extract PRSFC from ' // TRIM(RNAME)
                     WRITE( LOGDEV,'(5X,A,I7.7,1X,I6.6)')'For ',
     &               LDATE, LTIME
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF
#else
            IF( R_READY ) THEN
               IF(.NOT. NLDNSTRIKE ) THEN
                 IF ( .NOT. XTRACT3( RNAME, RC_NAME, 1, 1,
     &                       STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3,
     &                       LDATE, LTIME, LTNG_RC) ) THEN
                     XMSG = 'Could not extract LTNG_RC from ' // TRIM(RNAME)
                     CALL M3EXIT( PNAME, LDATE, LTIME, XMSG, XSTAT2 )
                 END IF
               END IF
               IF ( .NOT. XTRACT3( RNAME, 'PRSFC', 1, 1,
     &                       STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3, 
     &                       LDATE, LTIME, LTNG_PRSFC ) ) THEN
                     XMSG = 'Could not extract PRSFC from ' // TRIM(RNAME)
                     WRITE( LOGDEV,'(5X,A,I7.7,1X,I6.6)')'For ',
     &               LDATE, LTIME
                     CALL M3EXIT( PNAME, LDATE, LTIME, XMSG, XSTAT2 )
               END IF
               IF(.NOT. NLDNSTRIKE ) THEN
                  R_READY = .FALSE.
               END IF
            END IF
#endif

            IF( R_READY ) THEN
               IF ( NLDNSTRIKE ) THEN
                  IF ( .NOT. XTRACT3( ANAME, "NLDNstrk", 1, 1,
     &                          STRTROWMC3,ENDROWMC3,STRTCOLMC3,ENDCOLMC3, 
     &                             LDATE, LTIME, NLDN_STRIKE ) ) THEN
                        XMSG = 'Could not read NLDN_STRIKE from' // TRIM( ANAME )
                        CALL M3EXIT( PNAME, LDATE, LTIME, XMSG, XSTAT2 )
                  END IF
                     R_READY = .FALSE.
               END IF  ! NLDNSTRIKE
            END IF     !R_READY
            ALLOCATE( COLUMN_LTNG_NO( NCOLS,NROWS ), STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = 'COLUMN_LTNG_NO memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF
            COLUMN_LTNG_NO( :,: ) = 0.0
!C Iterate over each grid cell and distribute lightning NO vertically

               DO ROW = 1, MY_NROWS
                  DO COL = 1, MY_NCOLS
                     IF ( NLDNSTRIKE ) THEN
                        COLUMN_LTNG_NO( COL,ROW ) =
     &                     ( NLDN_STRIKE( COL,ROW )
     &                       * SQUAREKM     ! the NLD_STRIKE in unit of km-2
     &                       * OCEAN_MASK( COL,ROW ) ! reduce offshore strikes
     &                       * ( MOLSNCG + MOLSNIC * ICCG( COL,ROW ) ) )
     &                       / ( 60.0 * 60.0 ) ! get time units right:emisions are in the unit of moles/s
                     ELSE
                        IF ( LTNG_RC( COL,ROW ) .GT. 0 ) THEN
                           SCL_FACTOR = EXP(INTERCEPT_lg(COL, ROW))
                           IF ( LTNG_RC( COL,ROW ) .GT. SCL_FACTOR .AND. OCEAN_MASK( COL,ROW ) .GT. 0.2 )THEN
                              COLUMN_LTNG_NO( COL,ROW ) =
     &                         ( EXP( SLOPE_LG( COL,ROW ) * LOG( LTNG_RC( COL,ROW ) )
     &                           + INTERCEPT_LG( COL,ROW) )
     &                           * SQUAREKM   ! the relation is built on the unit of flash/km2*hr
     &                           * OCEAN_MASK( COL,ROW ) ! reduce offshore strikes
     &                           * ( MOLSNCG  ! moles N per flash intercloud strikes per cloud-to-ground strike
     &                             + ( MOLSNIC * ICCG( COL,ROW ) ) ) )
     &                         / ( 60.0 * 60.0 )         ! get time units right
                              IF ( COLUMN_LTNG_NO( COL,ROW ) .LT. 0 ) COLUMN_LTNG_NO( COL,ROW ) = 0.0
                           ELSE
                              COLUMN_LTNG_NO( COL,ROW ) =
     &                         ( ( SLOPE( COL,ROW ) * LTNG_RC( COL,ROW ) + INTERCEPT( COL,ROW ) )
     &                           * SQUAREKM   ! the relation is built on the unit of flash/km2*hr
     &                           * OCEAN_MASK( COL,ROW ) ! reduce offshore strikes
     &                           * ( MOLSNCG  ! moles N per flash intercloud strikes per cloud-to-ground strike
     &                             + ( MOLSNIC * ICCG( COL,ROW ) ) ) )
     &                         / ( 60.0 * 60.0 )         ! get time units right 
                              IF ( COLUMN_LTNG_NO( COL,ROW ) .LT. 0 ) COLUMN_LTNG_NO( COL,ROW ) = 0.0
                           END IF
                        ELSE
                           COLUMN_LTNG_NO( COL,ROW ) = 0.0
                        END IF
                     END IF   ! NLDNSTRIKE
                  END DO   ! COL
               END DO   ! ROW

            VDEMIS_LT = 0.0   ! array assignment

            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS

!C check to see if there are lightning strikes for this grid cell
!C only calculate lightning for cloud top greater than 6500 meters

                  IF ( COLUMN_LTNG_NO( COL,ROW ) .LE. 0.0 ) CYCLE
                  SUM_FRAC = 0.0
                  LTEMIS = 0.0   ! array assignment

                  DO LAY = 1, LTLYRS

!C Get pressures: Use SIGMA values and surface pres.
!p=sigma*(psfc-ptop)+ptop
                     BOTTOM = ( VGLVSLT( LAY-1 )
     &                      * ( LTNG_PRSFC( COL,ROW ) - VGTOP_GD )
     &                      + VGTOP_GD ) * CONVPA
!                           write( logdev,* ) "bottom: ", bottom
                     TOP    = ( VGLVSLT( LAY )
     &                      * ( LTNG_PRSFC( COL,ROW ) - VGTOP_GD )
     &                      + VGTOP_GD ) * CONVPA

!C Find the bottom and top of each layer, and calculate the fraction 
!C of the column emissions for that layer
!C Use normal distribution, mean = wmu, standard deviation = wsigma
                     inErfB      = ( BOTTOM - WMU ) / ( WSIGMA * SQRT2 )
                     inErfT      = ( TOP - WMU ) / ( WSIGMA * SQRT2 )
                     outErfB     = ERF( inErfB )
                     outErfT     = ERF( inErfT )
                     BOTTOM_FRAC = 0.5 * ( 1.0 + outErfB )
                     TOP_FRAC    = 0.5 * ( 1.0 + outErfT )

!C Find the bottom and top of each layer, and calculate the fraction
!C of the column emissions for that layer
!C use normal distribution, mean = wmu, standard deviation = wsigma
                     inErfB       = ( BOTTOM - W2MU ) / ( W2SIGMA * SQRT2 )
                     inErfT       = ( TOP - W2MU ) / ( W2SIGMA * SQRT2 )
                     outErfB      = ERF( inErfB )
                     outErfT      = ERF( inErfT )
                     BOTTOM_FRAC2 = 0.5 * ( 1.0 + outErfB )
                     TOP_FRAC2    = 0.5 * ( 1.0 + outErfT )

!C Add weighted contribution to this level
                     WEIGHT = ( BOTTOM_FRAC - TOP_FRAC )
     &                      + ( BOTTOM_FRAC2 - TOP_FRAC2 ) * 0.2

                     LTEMIS( LAY ) = WEIGHT * COLUMN_LTNG_NO( COL,ROW )

!C Sum weights in order to normalize to 1
                     SUM_FRAC = SUM_FRAC + WEIGHT

!C If emissions are less than 0, generate an error message in the log
                     IF ( LTEMIS( LAY ) .LT. 0.0 ) THEN
                        WRITE( LOGDEV,* ) LTEMIS( LAY ),
     &                                    COLUMN_LTNG_NO( COL,ROW ),
     &                                    BOTTOM_FRAC, TOP_FRAC
                        XMSG = '*** Ltng NO emis is less than zero'
                        CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                     END IF

                  END DO   ! end layers loop 

                  DO LAY = 1, LTLYRS
!C Re-normalize, in some cases area under the error function is not 1
                     VDEMIS_LT( COL,ROW,LAY ) = LTEMIS( LAY ) / SUM_FRAC
                  END DO        ! layers renormalized

               END DO   ! columns
            END DO   ! rows

!C Determine the time to read and/or write the hourly files
            NTICS = NTICS + 1
            TOT_TSTEP = TOT_TSTEP + TIME2SEC( TSTEP( 2 ) ) 
            LASTTIC = TOT_TSTEP .GE. LIGHTNING_FILE_TSTEP
            IF ( LASTTIC ) THEN
               MTICS = NTICS
               NTICS = 0
               TOT_TSTEP = 0
               CALL NEXTIME( LDATE, LTIME, TSTEP( 1 ) )
               R_READY = .TRUE.
            END IF

!C Write lightning NO to the diagnostic file
            IF ( LTNGDIAG ) THEN
               VDEMIS_DIAG = VDEMIS_DIAG + VDEMIS_LT   ! array assignment
               COLUMN_DIAG = COLUMN_DIAG + COLUMN_LTNG_NO
               IF ( LASTTIC ) THEN   ! time to write out
                  DIVFAC = 1.0 / REAL( MTICS, 4 )
                  VDEMIS_DIAG = VDEMIS_DIAG * DIVFAC   ! array assignment
                  COLUMN_DIAG = COLUMN_DIAG*DIVFAC
                  IF ( .NOT. WRITE3( CTM_LTNGDIAG_1, LTSPC, LDATE, LTIME, VDEMIS_DIAG ) )  THEN
                     XMSG = 'Could not write to ' // TRIM( CTM_LTNGDIAG_1 )
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                  ELSE
                     WRITE( LOGDEV,94040 )
     &                    'Timestep written to', TRIM( CTM_LTNGDIAG_1 ),
     &                    'for date and time', LDATE, LTIME
                  END IF
                  IF ( .NOT. WRITE3( CTM_LTNGDIAG_2, LTSPC, LDATE, LTIME, COLUMN_DIAG ) )  THEN
                     XMSG = 'Could not write to ' // TRIM( CTM_LTNGDIAG_2 )
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                  ELSE
                     WRITE( LOGDEV,94040 )
     &                    'Timestep written to', TRIM( CTM_LTNGDIAG_2 ),
     &                    'for date and time', LDATE, LTIME
                  END IF

                  VDEMIS_DIAG = 0.0   ! array assignment
                  COLUMN_DIAG = 0.0
               END IF ! LASTTIC

            END IF  ! diagnostics turned on

         ELSE  ! LTNGO is not "InLine", but instead specifies a file

!C Read in lightning NO production from an input file
            VDEMIS_LT = 0.0   ! array assignment

            IF ( .NOT. INTERPX( LNAME, LTSPC, PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,LTLYRS,
     &                          JDATE, JTIME, VDEMIS_LT ) ) THEN
               XMSG = 'Could not read NO' //
     &                ' from file "' // TRIM( LNAME ) // '"'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

         END IF  ! end lightning NO production inline or from a file

         WRITE( LOGDEV,94042 ) 'VDEMIS_LT sum:', JDATE, JTIME, SUM( VDEMIS_LT )

         DEALLOCATE( COLUMN_LTNG_NO )


         RETURN

C------------------  Format  Statements   ------------------------------

94040 FORMAT( /5X, 3( A, :, 1X ), I8, ":", I6.6 )
94042 FORMAT( /5X, A, 1X, I8, ":", I6.6, 1X, 1PE13.5 )

         END SUBROUTINE GET_LTNG

      END MODULE LTNG_DEFN
