
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
C $Header: /project/yoj/arc/BCON/src/profile/prof_driver.F,v 1.2 2011/10/21 16:52:36 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE INBNDY( LOGUNIT, SDATE, STIME, NHRS, NML_BC_NAME )   

C***********************************************************************
 
C  Function: Opens the input profile data file, and calls the routine to 
C             compute the BCs and write the output BC file  
              
C  Preconditions: None
  
C  Key Subroutines/Functions Called: GC_DPC_MAP
C                                    NGC_SPC_MAP
C                                    LST_SPC_MAP
C                                    BCOUT
 
C  Revision History:
C    Prototype created by Jerry Gipson, January, 1998
C    Modified April 1998 by JG to remove ppb option for input profile file
C    13 Jul 11 J.Young: Replaced I/O API include files with M3UTILIO
C                       and Namelist for species definitions
C    23 May 12 J.Young: Replaced BC_PARMS include file with an F90 module
 
C***********************************************************************

      USE M3UTILIO    ! IOAPI module
      USE CGRID_SPCS
      USE BC_PARMS    ! BCON parameters

      IMPLICIT NONE     

C Include Files: None

C Arguments: 
      INTEGER LOGUNIT                 ! Unit number for output log
      INTEGER SDATE                   ! Date for BC Output
      INTEGER STIME                   ! Time for BC output
      INTEGER NHRS                    ! No. of hours for BC output
      CHARACTER( 16 ) :: NML_BC_NAME( * )  ! Names of species in *.NMLs

C Parameters: None

C External Functions: None
      
C Local Variables:
      CHARACTER(  16 ) :: ENV_DFLT = ' ' ! Environment variable default value
      CHARACTER(  80 ) :: ENV_DESC       ! Environment variable description
      CHARACTER( 256 ) :: LINEIN         ! Input line
      CHARACTER(  80 ) :: MSG            ! Log message
      CHARACTER(  16 ) :: PNAME = 'INBNDY'          ! Program name
      CHARACTER( 256 ) :: PROF_FL_NAME   ! Input profile file name
      CHARACTER(  16 ) :: BC_PROFILE = 'BC_PROFILE' ! Logical name of input profile file
      CHARACTER(  16 ) :: SPEC_NAME      ! Species name

      CHARACTER(  16 ) :: PROF_SP_NAME( MXCTMS*MXVARS3 ) ! Species names on profile file
!     CHARACTER(  16 ) :: SPNAME_OUT( NSPCSD )       ! Output file species names
      CHARACTER(  16 ), ALLOCATABLE :: SPNAME_OUT( : )            ! Output file species names

      INTEGER INDX         ! Array indices for species
      INTEGER NEDGE        ! Loop index for no. of edges (4)  
      INTEGER N, SPC       ! Loop indices
      INTEGER N_GC_OUT     ! No. of GC species for conversion 
      INTEGER NLAYS_IN     ! No. of layers in input conc file
      INTEGER NSPCS_IN     ! CTM files species counter
      INTEGER NSPCS_OUT    ! Number of BC species on output file(s)
      INTEGER PFILE        ! Unit number of profile file    
      INTEGER STATUS       ! Status code
      INTEGER ALST         ! Allocate status

!     INTEGER :: SPMAP( NSPCSD )  ! Map to input file species from BC species
      INTEGER, ALLOCATABLE :: SPMAP( : )  ! Map to input file species from BC species

      LOGICAL LERROR                     ! Error flag

C***********************************************************************

      WRITE( LOGUNIT, 92000 )

      NHRS = 0
      SDATE = 0
      STIME = 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Open the input profile file, read the number of species and their 
c  names, and list the file on the output log
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ENV_DESC = 'Input BC Profile file'        
      CALL ENVSTR( BC_PROFILE, ENV_DESC, ENV_DFLT, PROF_FL_NAME, STATUS )

      IF ( STATUS .NE. 0 ) THEN
         MSG = 'No BC profile file assigned to BC_PROFILE'
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
      END IF

      PFILE = JUNIT()

      OPEN( UNIT = PFILE, FILE = PROF_FL_NAME , ERR = 999 )

      WRITE( LOGUNIT, 92020 ) PROF_FL_NAME( 1:100 )

C Report the 3-line header
      DO N = 1, 3
         READ( PFILE, 94000 ) LINEIN
         WRITE( LOGUNIT, 92040 ) TRIM( LINEIN )
      END DO

      READ( PFILE, * ) NLAYS_IN, NSPCS_IN
!     READ( PFILE, * ) NSPCS_IN, NLAYS_IN
      BACKSPACE PFILE 
      READ( PFILE, 94000 ) LINEIN
      WRITE( LOGUNIT, 92040 ) TRIM( LINEIN )

      READ( PFILE, 94000 ) LINEIN
      WRITE( LOGUNIT, 92040 ) TRIM( LINEIN )

      LERROR = .FALSE.
      IF ( NLAYS_IN .GT. MXLAYS3 + 1 ) THEN
         MSG = 'ERROR: No of layers on profile file exceeds MXLAYS + 1' 
         WRITE( LOGUNIT, 92080  ) MSG
         LERROR = .TRUE.
      END IF

      IF ( NSPCS_IN .GT. MX_INFL_SP ) THEN
         MSG = 'ERROR: No of layers on profile file exceeds  MXINFL_SP' 
         WRITE( LOGUNIT, 92080  ) MSG
         LERROR = .TRUE.
      END IF

      IF ( LERROR ) THEN
         MSG = 'Dimension problem in BC profile file'
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
      END IF

      LERROR = .FALSE.
      DO NEDGE = 1, 4   ! file order: North, East, South. West

         READ( PFILE, 94000 ) LINEIN
         WRITE( LOGUNIT, 92040 ) TRIM( LINEIN )

         DO SPC = 1, NSPCS_IN
            READ( PFILE, * ) SPEC_NAME
            BACKSPACE PFILE
            READ( PFILE, 94000 ) LINEIN
            IF ( NEDGE .EQ. 1 ) THEN
               PROF_SP_NAME( SPC ) = SPEC_NAME
               WRITE( LOGUNIT, 92040 ) TRIM( LINEIN )
            ELSE
               INDX = INDEX1( SPEC_NAME, NSPCS_IN, PROF_SP_NAME( SPC ) )
               IF ( INDX .EQ. 0 ) THEN
                  LERROR = .TRUE.
                  WRITE( LOGUNIT, 92050 ) TRIM( LINEIN )
               ELSE
                  WRITE( LOGUNIT, 92040 ) TRIM( LINEIN )
               END IF
            END IF
         END DO
                  
      END DO

      IF ( LERROR ) THEN
         MSG = 'Inconsistent species names in BC profile file'
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
      END IF
      
      WRITE( LOGUNIT, 92060 ) 

      CLOSE( PFILE )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Build a species map and report
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ALLOCATE ( SPMAP( NSPCSD ), SPNAME_OUT( NSPCSD ), STAT = ALST )
      IF ( ALST .NE. 0 ) THEN
         MSG = 'Error allocating SPMAP or SPNAME_OUT'
         CALL M3EXIT ( PNAME, 0, 0, MSG, XSTAT2 )
      END IF

      SPMAP = 0         ! array assignment
      SPNAME_OUT = ' '  ! array assignment

      CALL GC_SPC_MAP( LOGUNIT, NSPCS_IN, N_GC_OUT,
     &                 SPMAP, PROF_SP_NAME, SPNAME_OUT )

      CALL NGC_SPC_MAP( LOGUNIT, NSPCS_IN, N_GC_OUT, NSPCS_OUT,
     &                  SPMAP, PROF_SP_NAME, SPNAME_OUT )

      CALL LST_SPC_MAP( LOGUNIT, NSPCS_IN, NSPCS_OUT,
     &                  SPMAP, NML_BC_NAME,
     &                  PROF_SP_NAME, SPNAME_OUT )
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Call the routine that computes and outputs BCs from the input
c  CTM conc file(s)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL BCOUT( LOGUNIT,
     &            SDATE, STIME, NHRS,
     &            NSPCS_IN, NSPCS_OUT, SPMAP,
     &            PROF_FL_NAME,
     &            PROF_SP_NAME,
     &            SPNAME_OUT )

      RETURN

999   CONTINUE

      MSG = 'Could not open file ' // BC_PROFILE
      CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )

C************************* FORMAT STATEMENTS ***************************

92000 FORMAT( // 1X, 79( '#' ) 
     &         / 1X, '#  Input file section '
     &         / 1X, 79( '#' )
     &        // 5X, 'Boundary Concentrations from an input BC profile ',
     &               'file.'
     &         / 5X, 'BCs do not vary with time.' )

92020 FORMAT(  / 5X, 'Profile file name: ', A 
     &         / 5X, 'File contents: '/ )

92040 FORMAT(  7X, A )
92050 FORMAT(  2X, '---> ', A )
92060 FORMAT(  / 5X, 'End of file contents.')
92080 FORMAT( // 1X, A )
92200 FORMAT( // 5X, 'Profile concentration units assumed to be ppm' ) 
92220 FORMAT( // 5X, 'Profile concentration units assumed to be ppb' ) 
94000 FORMAT( A )

      END
