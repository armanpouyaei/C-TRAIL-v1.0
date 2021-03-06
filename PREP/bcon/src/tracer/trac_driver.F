
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
C $Header: /project/yoj/arc/BCON/src/tracer/trac_driver.F,v 1.2 2011/10/21 16:52:37 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%


      SUBROUTINE INBNDY( LOGUNIT, SDATE, STIME, RUNLEN )   

C***********************************************************************
 
C Function: Call the routine to compute catalogued special tracer species
C           BCs and write the Models-3 output file
              
C Preconditions: None
  
C Key Subroutines/Functions Called:
 
C Revision History:
C    Prototype created by Jerry Gipson, January, 1998          
C    02/25/02 Steve Howard (Jeff Young) - dynamic allocation
C    19 Jul 11 J.Young: Convert for Namelist redesign + M3UTILIO
C    23 May 12 J.Young: Replaced BC_PARMS include file with an F90 module
                    
C***********************************************************************

      USE HGRD_DEFN    ! Module to store and load the horizontal grid variables
      USE VGRD_DEFN    ! vertical layer specifications
      USE M3UTILIO     ! IOAPI
      USE CGRID_SPCS   ! CGRID mechanism species
      USE BC_PARMS     ! BCON parameters

      IMPLICIT NONE     

C Include Files:

C Arguments: 
      INTEGER LOGUNIT                 ! Unit number for output log
      INTEGER SDATE                   ! Date for IC Output
      INTEGER STIME                   ! Time for IC output
      INTEGER RUNLEN                  ! No. of hours for BC output

C Parameters: None

C External Functions: 
      INTEGER FINDEX    ! Finds the index of a number in a list 
      
C Local Variables:
      CHARACTER( 16 ) :: BLANK  = ' '     ! String of blanks
      CHARACTER( 80 ) :: MSG              ! Log message
      CHARACTER( 16 ) :: PNAME = 'INBNDY' ! Procedure name
      CHARACTER( 16 ) :: SUROUT           ! Surrogate output to log
      CHARACTER( 16 ) :: VNAME            ! Output variable name

      CHARACTER( 16 ) :: SPNAME_OUT( N_TR_SPC )    ! Output tracer names

      CHARACTER( 16 ) :: TRNML_SPC_NF( N_TR_SPC )  ! TR_SPC.NML species not found
      CHARACTER( 16 ) :: TRNML_SUR_NF( N_TR_ICBC ) ! TR_SPC.NML surr not found
      CHARACTER( 16 ) :: TR_FNAME( MXCTMS )        ! Output file name
                            
      INTEGER FLN          ! Output file number
      INTEGER ICBC_IND     ! Index of surrogate in TR_ICBC
      INTEGER N            ! Loop index
      INTEGER NSP          ! Counter
      INTEGER NSPCS_OUT    ! Number of IC species on output file(s)
      INTEGER TRIND        ! Index for catalogued tracer
      INTEGER TSTEP        ! Time step for output file
      INTEGER ALST         ! Allocate status

      INTEGER :: TRMAP( N_TR_SPC ) ! Map to get catalogued tracer for NML species
      LOGICAL :: LUSED( N_TR_SPC ) ! Flag for output for TR_SPC species

!     REAL    :: TCOUT( NBNDY,NLAYS )   ! Output tracer concs
      REAL, ALLOCATABLE :: TCOUT( :,: )   ! Output tracer concs

C***********************************************************************

      WRITE( LOGUNIT, 92000 )
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Build a species map and report
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NSPCS_OUT = 0
      DO N = 1, N_TR_SPC

         LUSED( N ) = .TRUE.
         TRMAP( N ) = 0

C get the index of the NML TR_SPC name in the catalogued list
         TRIND = INDEX1( TR_SPC( N ), NTRACERS, TRAC_NAME )
         TRMAP( N )  = TRIND
         SUROUT = BLANK

         IF ( TRIND .LE. 0 ) THEN   ! look for a surrogate
C get the index of the NML TR_ICBC surrogate name in the catalogued list
            ICBC_IND  = FINDEX( N, N_TR_ICBC, TR_ICBC_MAP )

            IF ( ICBC_IND .NE. 0 ) THEN
               TRIND = INDEX1( TR_ICBC( ICBC_IND ), NTRACERS, TRAC_NAME )
               TRMAP( N ) = TRIND
               SUROUT = TR_ICBC( ICBC_IND )
            END IF
         END IF
         
C list the ones found and tag the ones not found
         IF ( TRMAP( N ) .GT. 0 .AND. TRMAP( N ) .LE. NTRACERS ) THEN
            WRITE( LOGUNIT, 92020 ) TR_SPC( N ), SUROUT, 
     &                            TRAC_NAME( TRMAP( N ) ), 
     &                            TRAC_DESC( TRMAP( N ) )
            NSPCS_OUT = NSPCS_OUT + 1
            SPNAME_OUT( NSPCS_OUT ) = TRAC_NAME( TRMAP( N ) )
         ELSE
            LUSED( N ) = .FALSE.
         END IF

      END DO

      IF ( NSPCS_OUT .EQ. 0 ) THEN
          MSG = 'ERROR: No NML species could be mapped to a'
     &           // ' catalogued tracer species.'
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
      END IF 

C list *.NML species with no BC output
      NSP = 0
      DO N = 1, N_TR_SPC
         IF ( .NOT. LUSED( N ) ) THEN
            NSP = NSP + 1
            TRNML_SPC_NF( NSP ) = TR_SPC( N )
         END IF
      END DO

      IF ( NSP .GT. 0 ) THEN
         WRITE( LOGUNIT, 92040 )
         WRITE( LOGUNIT, 92060 ) ( TRNML_SPC_NF( N ), N = 1, NSP )
      END IF
         
C list surrogate species not found on the 
      NSP = 0
      DO N = 1, N_TR_ICBC
         TRIND = INDEX1( TR_ICBC( N ), NTRACERS, TRAC_NAME )
         IF ( TRIND .EQ. 0 ) THEN
            NSP = NSP + 1
            TRNML_SUR_NF( NSP ) = TR_ICBC( N )
         END IF
      END DO

      IF ( NSP .GT. 0 ) THEN
         WRITE( LOGUNIT, 92080 )
         WRITE( LOGUNIT, 92060 ) ( TRNML_SUR_NF( N ), N = 1, NSP )
      END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  open the Models3 IC output file(s)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      TSTEP = 0
      CALL OPN_BC_FILE( LOGUNIT, SDATE, STIME, TSTEP, NSPCS_OUT, 
     &                  SPNAME_OUT, TR_FNAME, 1 )  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute the BCs and then write the ooutput file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ALLOCATE ( TCOUT( NBNDY,NLAYS ), STAT = ALST )
      IF ( ALST .NE. 0 ) THEN
         MSG = 'Error allocating TCOUT'
         CALL M3EXIT ( PNAME, 0, 0, MSG, XSTAT2 )
      END IF

      DO N = 1, NSPCS_OUT

         FLN = ( N - 1 ) / MXVARS3 + 1
         VNAME = SPNAME_OUT( N )

         CALL TRAC_BC( VNAME, TCOUT )

         IF ( .NOT. WRITE3( TR_FNAME( FLN ), VNAME, SDATE, STIME,
     &                      TCOUT ) ) THEN
            MSG =  'Could not WRITE species ' //  VNAME // 
     &             'to file ' // TR_FNAME( FLN ) 
            CALL M3EXIT( PNAME, SDATE, STIME, MSG, XSTAT1 )
         END IF

      END DO

      RETURN

C************************* FORMAT STATEMENTS ***************************

92000 FORMAT( // 1X, 79( '#' ) 
     &         / 1X, '#  TRACER species section '
     &         / 1X, 79( '#' )
     &        // 5X, 'The tracer species module has been invoked; ',
     &               'BCs will be generated for tracer species ONLY.'
     &        // 5X, 'BCs set for the following tracer species:'
     &        //10X, 'NML Species       ICBC Surrogate    ',
     &               'Catalogued Tracer' )

92020 FORMAT(   10X, A16, 2X, A16, 2X, A16,': ', A )

92040 FORMAT( // 5X, 'WARNING: The following *.NML species could not be '
     &               'mapped to a catalogued tracer.'
     &         / 5X, '         BCs were not generated for them. ' )

92060 FORMAT(   (15X, 4( A, 2X ) ) )


92080 FORMAT( // 5X, 'WARNING: The following *.NML surrogate species ',
     &               'are not catalogued tracers.'
     &         / 5X, '         BCs were generated only when the *.NML '
     &               'species could be mapped to a catalogued tracer.')


94000 FORMAT( A )

      END
