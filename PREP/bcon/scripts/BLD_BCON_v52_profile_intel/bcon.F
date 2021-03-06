
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
C $Header: /project/yoj/arc/BCON/src/common/bcon.F,v 1.2 2011/10/21 16:52:32 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      PROGRAM BCON   

C***********************************************************************
 
C Function: To produce a Models-3 CTM Boundary Concentration file.
C           Valid input options include: 
C           1) a Models-3 CTM output file,
C           2) a user specified ASCII file of vertical profiles,
C           3) a RADM model output conc file,
C           4) subroutines that produce catalogued tracer species concs.
              
C Preconditions: None
  
C Key Subroutines/Functions Called: INBNDY
 
C Revision History:
C    Created by Jerry Gipson, January, 1998; based on the original Models-3
C                                            prototype.
C    Modified 4/22/99 by JG to provide default treatment for missing RUNLEN,
C                     SDATE, and STIME
C    02/25/02 Steve Howard (Jeff Young) - dynamic allocation
C    01/05/05 J.Young: vert dyn alloc - Use VGRD_DEFN
C    06 Jul 11 J.Young: Replaced I/O API include files with M3UTILIO
C    13 Jul 11 J.Young: Convert for Namelist redesign
C    23 May 12 J.Young: Replaced BC_PARMS include file with an F90 module
      
C***********************************************************************

      USE HGRD_DEFN     ! Module to store and load the horizontal grid variables
      USE VGRD_DEFN     ! vertical layer specifications
      USE M3UTILIO      ! IOAPI module
      USE CGRID_SPCS    ! CGRID mechanism species
      USE BC_PARMS      ! BCON parameters

      IMPLICIT NONE   

C Include Files: None

C Arguments: None

C Parameters: None

C External Procedures: None

C Local Variables:
!     CHARACTER( 16 ) :: NML_BC_NAME( NSPCSD )  ! Names of species in *.NMLs
      CHARACTER( 16 ), ALLOCATABLE :: NML_BC_NAME( : )  ! Names of species in *.NMLs
      CHARACTER( 16 ) :: PNAME = 'BCON'         ! Procedure name
      CHARACTER( 80 ) :: MSG                    ! Log message 

      INTEGER    LOGUNIT     ! Unit number for output log
      INTEGER    ISP         ! Species counter
      INTEGER    N           ! Loop index
      INTEGER    RUNLEN      ! No. of hours for BC output
      INTEGER :: SDATE = 0   ! Start date 
      INTEGER :: STIME = 0   ! Start time 
      INTEGER    STATUS      ! Status code
      INTEGER, SAVE :: NPROCS = 1 ! dummy for non-parallel

C***********************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Initialize IO/API and get environment variables
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      LOGUNIT = INIT3( )

      IF ( NSPCSD .LE. 0 ) THEN
         MSG = 'Number of species in *.NML files must be greater than zero'
         CALL M3EXIT ( PNAME, SDATE, STIME, MSG, XSTAT2 )
      END IF

      SDATE = ENVINT( 'SDATE', 'START DATE', 0, STATUS )
      IF ( STATUS .NE. 0 ) SDATE = 0

      STIME = ENVINT( 'STIME', 'START TIME', 0, STATUS )
      IF ( STATUS .NE. 0 ) STIME = 0

      RUNLEN = ENVINT( 'RUNLEN', 'RUN LENGTH', 0, STATUS )
      IF ( STATUS .NE. 0 ) RUNLEN = 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Load horizontal grid variables calling GRID_DEF->HGRD_INIT function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( .NOT. HGRD_INIT() ) then
         MSG = '*** Failure defining horizontal domain'
         CALL M3EXIT ( PNAME, SDATE, STIME, MSG, XSTAT2 )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c define vertical layer structure
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( .NOT. VGRD_INIT() ) THEN
         MSG = '*** Failure defining vertical structure'
         CALL M3EXIT ( PNAME, SDATE, STIME, MSG, XSTAT2 )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Set CMAQ CGRID mechanism species
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( .NOT. CGRID_SPCS_INIT() ) THEN
         MSG = 'Error in CGRID_SPCS:CGRID_SPCS_INIT'
         CALL M3EXIT ( PNAME, SDATE, STIME, MSG, XSTAT2 )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Build one array containing the names of the species in *.NML files 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ALLOCATE ( NML_BC_NAME( NSPCSD ), STAT = STATUS )
      IF ( STATUS .NE. 0 ) THEN
         MSG = 'Error allocating NML_BC_NAME'
         CALL M3EXIT ( PNAME, 0, 0, MSG, XSTAT1 )
      END IF

      ISP = 0
      DO N = 1, N_GC_SPC
         ISP = ISP + 1
         NML_BC_NAME( ISP ) = GC_SPC( N )
      END DO

      DO N = 1, N_AE_SPC
         ISP = ISP + 1
         NML_BC_NAME( ISP ) = AE_SPC( N )
      END DO

      DO N = 1, N_NR_SPC
         ISP = ISP + 1
         NML_BC_NAME( ISP ) = NR_SPC( N )
      END DO

      DO N = 1, N_TR_SPC
         ISP = ISP + 1
         NML_BC_NAME( ISP ) = TR_SPC( N )
      END DO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Generate the BC output file and quit
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL INBNDY( LOGUNIT, SDATE, STIME, RUNLEN, NML_BC_NAME )

      IF ( SHUT3() ) THEN
         WRITE ( *,92999 )
     &  '>>---->  Program  BCON completed successfully  <----<< '
      ELSE
         WRITE ( *,91000 )
     &      'FATAL ERROR shutting down Models-3 I/O'
      END IF

      STOP

C************************* FORMAT STATEMENTS ***************************

92999 FORMAT ( //5X, A, // )

91000 FORMAT ( //5X , '*** ERROR ABORT in program BCON ***',
     &            /5X , A , // )        !  generic error message format

      END
