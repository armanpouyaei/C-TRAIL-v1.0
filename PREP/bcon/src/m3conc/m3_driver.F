
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
C $Header: /project/yoj/arc/BCON/src/m3conc/m3_driver.F,v 1.2 2011/10/21 16:52:35 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE INBNDY( LOGUNIT, SDATE, STIME, RUNLEN, NML_BC_NAME )   

C***********************************************************************
 
C  Function: Opens the Models-3 CTM conc file(s) and checks for 
C            data consistency, and calls the routine to 
C            compute the BCs and write the output BC file  
              
C  Preconditions: None
  
C  Key Subroutines/Functions Called: CK_CTM_FLS
C                                    GC_SPC_MAP
C                                    NGC_SPC_MAP
C                                    LST_SPC_MAP
C                                    BCOUT  
 
C  Revision History:
C    Prototype created by Jerry Gipson, January, 1998          
C    Modified 4/22/99 by JG to provide default treatment for missing RUNLEN
C    Modified by JG May, 1999 to treat PinG concs
C    02/25/02 Steve Howard (Jeff Young) - dynamic allocation
C    01/05/05 J.Young: vert dyn alloc - Use VGRD_DEFN
C    13 Jul 11 J.Young: Replaced I/O API include files with M3UTILIO and
C                       Namelist for species definitions
C    23 May 12 J.Young: Replaced BC_PARMS include file with an F90 module
                    
C***********************************************************************

      USE HGRD_DEFN   ! Module to store and load the horizontal grid variables
      USE VGRD_DEFN   ! vertical layer specifications
      USE M3UTILIO    ! IOAPI module
      USE CGRID_SPCS
      USE BC_PARMS    ! BCON parameters

      IMPLICIT NONE     

C Include Files: None

C Arguments:
      INTEGER LOGUNIT            ! Unit number for output log
      INTEGER SDATE              ! Date for BC output
      INTEGER STIME              ! Time for BC output
      INTEGER RUNLEN             ! No. of hours for BC output
      CHARACTER( 16 ) :: NML_BC_NAME( * )  ! Names of species in *.NMLs

C Parameters: None

C External Functions: None

C Local Variables:
      CHARACTER(  16 ) :: ENV_DFLT       ! Environment variable default value
      CHARACTER(  80 ) :: ENV_DESC       ! Environment variable description
      CHARACTER(  16 ) :: FL_NAME        ! Input CTM file name
      CHARACTER(  80 ) :: MSG            ! Log message
      CHARACTER(  16 ) :: PNAME = 'INBNDY' ! Procedure name
      CHARACTER( 256 ) :: RET_VAL        ! Returned value of environment variable
      CHARACTER(  80 ) :: VARDESC        ! environment variable description

      CHARACTER(  16 ) :: CTM_FL_NAME( MXCTMS )           ! Names of CTM files
      CHARACTER(  16 ) :: INFL_SP_NAME( MXCTMS*MXVARS3 ) ! Names of CTM species
!     CHARACTER(  16 ) :: SPNAME_OUT( NSPCSD )             ! Output file species names
      CHARACTER(  16 ), ALLOCATABLE :: SPNAME_OUT( : )

      INTEGER JDATE        ! Current date
      INTEGER JTIME        ! Current time
      INTEGER N            ! File loop index
      INTEGER N_GC_OUT     ! No. of GC species
      INTEGER N_CTM_FLS    ! Number of input CTM files
      INTEGER NCOLS_IN     ! No. of columns in input conc file
      INTEGER NLAYS_IN     ! No. of layers in input conc file
      INTEGER NROWS_IN     ! No. of rows in input conc file
      INTEGER NSPCS_IN     ! CTM files species counter
      INTEGER NSPCS_OUT    ! Number of IC species on output file(s)
      INTEGER NSPCS        ! Number of NML species
      INTEGER NSTINC       ! No. of steps on file before start of output
      INTEGER STATUS       ! Status code
      INTEGER STRTINC      ! No. of sec. on file before start of output
      INTEGER STEPSECS     ! Seconds per time step
      INTEGER TOTSECS      ! Run duration on file seconds
      INTEGER NSTEPS       ! Run duration, as number of output time steps
      INTEGER TSTEP        ! Time step on CTM_CONC_1 file
      INTEGER V            ! Variable loop index
      INTEGER ALST         ! Allocate status

!     INTEGER :: SPMAP( NSPCSD )  ! Map to input file species from BC species
      INTEGER, ALLOCATABLE :: SPMAP( : )   ! Map to input file species from BC species

      LOGICAL LCHKD                        ! Flag for met file checked

C***********************************************************************

      WRITE( LOGUNIT, 92000 )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Determine the number of input CTM conc files that need to be read
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      N_CTM_FLS = 0

      WRITE( LOGUNIT, 92020 )

      DO N = 1, MXCTMS

         WRITE( FL_NAME, '( ''CTM_CONC_'', I1 )' ) N
         WRITE( ENV_DESC, '( ''CTM Concentration file no. '', I1 )' ) N 
         ENV_DFLT = ' '        
         CALL ENVSTR( FL_NAME, ENV_DESC, ENV_DFLT, RET_VAL, STATUS)

         IF ( STATUS .EQ. 0 ) THEN
            N_CTM_FLS = N_CTM_FLS + 1
            CTM_FL_NAME( N_CTM_FLS ) = FL_NAME
         END IF
     
      END DO

      IF ( N_CTM_FLS .EQ. 0 ) THEN
         MSG = 'No CTM CONC files found'
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open the input Models-3 CTM concentration file(s) and put the species 
c  names on each file in one contiguous array
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NSPCS_IN = 0

      DO N = 1, N_CTM_FLS

         IF ( .NOT. OPEN3( CTM_FL_NAME( N ), FSREAD3, PNAME ) ) THEN
            MSG = 'Could not open ' // CTM_FL_NAME( N ) // ' file'
            CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT2 )
         END IF

         IF ( .NOT. DESC3( CTM_FL_NAME( N ) ) ) THEN
            MSG = 'Could not read DESC of  ' // CTM_FL_NAME( N ) 
     &            // ' file'
            CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT2 )
         END IF

         DO V = 1, NVARS3D
            NSPCS_IN = NSPCS_IN + 1
            INFL_SP_NAME( NSPCS_IN ) = VNAME3D( V )
         END DO

      END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check multiple input CTM files for consistency
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( N_CTM_FLS .GT. 1 ) CALL CK_CTM_FLS( LOGUNIT, N_CTM_FLS,
     &                                         CTM_FL_NAME )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Reset start and runlen times if necessary
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( .NOT. DESC3( CTM_FL_NAME( 1 ) ) ) THEN
         MSG = 'Could not read DESC of  ' // CTM_FL_NAME( N ) 
     &         // ' file'
         CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT2 )
      END IF

      NCOLS_IN = NCOLS3D
      NROWS_IN = NROWS3D
      NLAYS_IN = NLAYS3D
      TSTEP    = TSTEP3D

      IF ( SDATE .EQ. 0 ) THEN
         MSG = 'Input SDATE equals zero; SDATE & STIME ' //
     &         'reset to CTM file start date & time.'
         CALL M3MESG( MSG )
         SDATE = SDATE3D
         STIME = STIME3D
      END IF

      IF ( .NOT. CHECK3( CTM_FL_NAME( 1 ), 'ALL', SDATE, STIME ) ) THEN
         MSG = 'Requested starting time is not on the ' // 
     &          CTM_FL_NAME( N ) // ' file'
         CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT2 )
      END IF
      
c Find number of steps on file from starting point
      STEPSECS = TIME2SEC( TSTEP3D )
      STRTINC = SECSDIFF( SDATE3D, STIME3D, SDATE, STIME )
      NSTINC = STRTINC / STEPSECS 

c Check and reset RUNLEN if necessary
      IF ( RUNLEN .NE. 0 ) THEN

         IF ( RUNLEN .LT. 1000000 ) THEN
            TOTSECS  = TIME2SEC( RUNLEN )
         ELSE                          ! HH > 99
            RUNLEN = RUNLEN - 1000000
            TOTSECS  = TIME2SEC( RUNLEN )
            TOTSECS  = TOTSECS + 360000
         END IF

         IF ( MOD( TOTSECS, STEPSECS ) .EQ. 0 ) THEN
            NSTEPS = TOTSECS / STEPSECS + 1
         ELSE
            MSG = 'Output time step ' // HHMMSS( TSTEP3D  ) //
     &            ' does not divide duration ' // HHMMSS( RUNLEN )
            CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT2 )
         END IF

         IF ( NSTEPS .GT. MXREC3D - NSTINC ) THEN
            MSG = 'Input RUNLEN exceeds time steps on input file. '
            WRITE( LOGUNIT, '( /5X, A )' ) MSG
            MSG = 'Resetting RUNLEN to correspond to CTM file ending ' //
     &            'date & time.' 
            WRITE( LOGUNIT, '( 5X, A )' ) MSG
            NSTEPS = MXREC3D - NSTINC
            RUNLEN = SEC2TIME( ( NSTEPS - 1 ) * STEPSECS )
         END IF

      ELSE

         MSG = 'Input RUNLEN not set or equal to zero. '
         WRITE( LOGUNIT, '( /5X, A )' ) MSG
         MSG = 'Resetting RUNLEN to correspond to CTM file ending ' //
     &         'date & time.'       
         WRITE( LOGUNIT, '( 5X, A )' ) MSG
         NSTEPS = MXREC3D - NSTINC

      END IF

c Check to make sure end date is on file     
      JDATE = SDATE
      JTIME = STIME
      CALL NEXTIME( JDATE, JTIME, RUNLEN )

      IF ( .NOT. CHECK3( CTM_FL_NAME( 1 ), 'ALL', JDATE, JTIME ) ) THEN
         MSG = 'Requested ending time is not on the ' // 
     &         CTM_FL_NAME( 1 ) // ' file'
         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT2 )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Check MET_CRO_3D and MET_BDY_3 if necessary
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LCHKD = .FALSE.

      IF ( VGTYP_GD .NE. VGTYP3D .OR. VGTOP_GD .NE. VGTOP3D ) THEN 
         LCHKD = .TRUE.
         CALL CK_MET_FL( LOGUNIT, SDATE, STIME, TSTEP, NSTEPS, CTM_FL_NAME )
         CALL CK_BDY_FL( LOGUNIT, SDATE, STIME, TSTEP, NSTEPS )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Build a species map and report 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NSPCS = N_GC_SPC + N_AE_SPC + N_NR_SPC
      ALLOCATE ( SPMAP( NSPCS ), SPNAME_OUT( NSPCS ), STAT = ALST )
      IF ( ALST .NE. 0 ) THEN
         MSG = 'Error allocating SPMAP or SPNAME_OUT'
         CALL M3EXIT ( PNAME, 0, 0, MSG, XSTAT2 )
      END IF

      CALL GC_SPC_MAP( LOGUNIT, NSPCS_IN, N_GC_OUT,
     &                 SPMAP, INFL_SP_NAME, SPNAME_OUT )

      CALL NGC_SPC_MAP( LOGUNIT, NSPCS_IN, N_GC_OUT, NSPCS_OUT,
     &                  SPMAP, INFL_SP_NAME, SPNAME_OUT )

      CALL LST_SPC_MAP( LOGUNIT, NSPCS_IN, NSPCS_OUT,
     &                  SPMAP, NML_BC_NAME, INFL_SP_NAME, SPNAME_OUT )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute and outputs BCs from the input CTM conc file(s)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL BCOUT( LOGUNIT,
     &            N_CTM_FLS,
     &            SDATE, STIME, NSTEPS,
     &            NCOLS_IN, NROWS_IN, NLAYS_IN,
     &            NSPCS_IN, NSPCS_OUT, SPMAP,   
     &            CTM_FL_NAME,
     &            INFL_SP_NAME,
     &            SPNAME_OUT )

      RETURN

C************************* FORMAT Statements ***************************

92000 FORMAT( // 1X, 79( '#' ) 
     &         / 1X, '#  Input section '
     &         / 1X, 79( '#' )
     &        // 5X, 'Boundary Concentrations from a CMAQ CTM ',
     &               'concentration file.' )

92020 FORMAT( // 5X, 'Input Models3 CTM file names: ' )

      END
