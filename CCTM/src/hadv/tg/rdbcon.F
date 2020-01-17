
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
C $Header: /project/yoj/arc/CCTM/src/hadv/yamo/rdbcon.F,v 1.7 2012/01/19 14:51:13 yoj Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE RDBCON ( JDATE, JTIME, TSTEP)

C-----------------------------------------------------------------------
C Function:
C   Read boundary concentrations data for advection and couple with
C   SqRDMT, Where SqRDMT = Sq. Root [det ( metric tensor )]
C                        = Vertical Jacobian / (map scale factor)**2
C   Load Air Density X SqRDMT = RHOJ into last BCON slot for advection

C Preconditions:

C Subroutines and Functions Called:
C   INTERPX, INTERPB, M3EXIT, TRIMLEN, ADVBC_MAP, TIME2SEC, SEC2TIME, NEXTIME

C Revision History:
C   Jeff - Aug 1997 Based on beta version, keep in ppmV units (not aerosols)
C   Jeff - Dec 97 - add CMIN
C   Jeff - Apr 98 - fix conversion/coupling for aerosol number species
C   Jeff - Apr 01 - dyn alloc - Use PINTERB for boundary data - assume the met
C                   data could come from a larger file, but not the conc cata
C   23 Jun 03 J.Young: for layer dependent advection tstep
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   21 Jun 10 J.Young: convert for Namelist redesign
C   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
C                      removed deprecated TRIMLEN
C   08 Jun 12 J.Young: remove full character blank padding for GNU Fortran (GCC) 4.1.2
C   14 Sep 15 J.Young: change Fortran string comparison to logical from advbc_map
C-----------------------------------------------------------------------

      USE GRID_CONF            ! horizontal & vertical domain specifications
      USE CGRID_SPCS           ! CGRID mechanism species
      USE UTILIO_DEFN
      USE RXNS_DATA, ONLY : MECHNAME
      USE AERO_DATA, only : CHECK_AERO_ICBC, N_MODE
	  USE PACKETS
#ifdef parallel
      USE SE_MODULES            ! stenex (using SE_COMM_MODULE, SE_UTIL_MODULE)
      USE SWAP_SANDIA
#else
      USE NOOP_MODULES          ! stenex (using NOOP_COMM_MODULE, NOOP_UTIL_MODULE)
#endif
	  
      IMPLICIT NONE
 
C Includes:
      INCLUDE SUBST_CONST      ! constants
      INCLUDE SUBST_FILES_ID   ! file name parameters
	  INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction
!	  include SUBST_AE_ADV

C Arguments:
      INTEGER, INTENT( IN ) :: JDATE       ! current model date, coded YYYYDDD
      INTEGER, INTENT( IN ) :: JTIME       ! current model time, coded HHMMSS
      INTEGER, INTENT( IN ) :: TSTEP       ! timestep
!      INTEGER, INTENT( IN ) :: LVL         ! layer
!     REAL           BCON( NBNDY,NLAYS,* )  ! boundary concentrations
      REAL, ALLOCATABLE, SAVE :: BCON( :,:,: ) ! boundary concentrations
!      LOGICAL, INTENT( INOUT ) :: LSTAT       ! Flag for identifying Aerosol BC flaws
!	  REAL  :: BCON(:,:) ! boundary conditions	  

C Parameters:
      REAL, PARAMETER :: KGPMG = 1.0E-09  ! Kg / micro-gram
      REAL, PARAMETER :: CMIN = 1.0E-30

      CHARACTER( 16 ) :: CONCMIN

C File variables:
!     REAL        JACOBM_BND( NBNDY,NLAYS )        ! boundary Jacobian
!     REAL        RHOJ_BND  ( NBNDY,NLAYS )        ! mid-layer boundary RhoJ
      REAL, ALLOCATABLE, SAVE :: JACOBM_BND( :,: ) ! boundary Jacobian
      REAL, ALLOCATABLE, SAVE :: RHOJ_BND  ( :,: ) ! mid-layer boundary RhoJ

      REAL, ALLOCATABLE, SAVE :: JACOBM_BUF( :,: ) ! grid Jacobian
      REAL, ALLOCATABLE, SAVE :: RHOJ_BUF  ( :,: ) ! mid-layer grid RhoJ

      REAL, ALLOCATABLE, SAVE :: BBUF      ( :,:,: ) ! bcon file buffer
      REAL, ALLOCATABLE, SAVE :: BBUF_TMP  ( :,: ) ! bcon file buffer

!      REAL :: BCON(1:NBNDY, 1:NLAYS, 1:NSPCSD) ! boundary conditions
!      REAL :: BBUF(1:NBNDY, 1:NLAYS, 1:N_SPC_ADV) ! bcon file buffer	  

C External Functions: None

C Local variables:
      CHARACTER( 16 ) :: PNAME = 'RDBCON'
      CHARACTER( 16 ) :: VNAME

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 16 ) :: BLNK = ' '
      CHARACTER( 16 ) :: BCNAME_TMP = ' '
	  CHARACTER( 16 ) :: PNA
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: BCNAME( : )  ! BC name for adv species
      REAL,            ALLOCATABLE, SAVE :: BCFAC( : )   ! Scale factor for BCs
      LOGICAL,         ALLOCATABLE, SAVE :: ONFILE( : )  ! BC species on the file

      CHARACTER( 96 )  :: XMSG = ' '
      CHARACTER( 199 ) :: XMSG2 = ' '
      INTEGER          :: STAT

      INTEGER, SAVE :: N_SPC_ADV
      INTEGER       :: LMODE    !Identifies the problematic mode from
                                !the BC Check routine
      REAL          :: AER_PAR( 2, N_MODE,5 )  !Modal parameter after the BC 
                                               !check (N, dg, sg)
                                               !      (N, M2, M3) - Before
                                               !      (N, M2, M3) - After

      INTEGER   MDATE             ! mid-advection date
      INTEGER   MTIME             ! mid-advection time
      INTEGER   STEP              ! advection time step in seconds

      INTEGER   BND, VAR, SPC, A2C     ! loop counters
      INTEGER   COL, ROW, IMODE   ! loop counters
      INTEGER   STRT, FINI, IAER, J
      INTEGER   ALLOCSTAT

      INTEGER   COUNT             ! Counter for constructing boundary arrays

      INTEGER   GXOFF, GYOFF      ! global origin offset from file
      LOGICAL, SAVE :: WINDOW = .FALSE. ! posit same file and global
                                        ! processing domain
										
      INTEGER :: C, R, L        ! column, row and layer indices
      INTEGER :: K              ! boundary cell index			
	  
!      INTEGER :: ADV_TO_CGRID_MAP(N_SPC_ADV) ! global ADV to CGRID map
	  INTEGER, ALLOCATABLE, SAVE :: ADV_MAP( : )  ! global adv map to CGRID	  

      INTEGER, SAVE :: GC_ADV_FINI 
      INTEGER, SAVE :: AE_ADV_FINI 
      INTEGER, SAVE :: NR_ADV_FINI 
      INTEGER, SAVE :: TR_ADV_FINI 

      INTEGER, SAVE :: GC_ADV_STRT 
      INTEGER, SAVE :: AE_ADV_STRT 
      INTEGER, SAVE :: NR_ADV_STRT 
      INTEGER, SAVE :: TR_ADV_STRT

 	  INTEGER, SAVE :: ASPC                     ! RHOJ index in CGRID

C for INTERPX
      INTEGER, SAVE :: STRTCOL,   ENDCOL,   STRTROW,   ENDROW
      INTEGER       :: STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC
	  

      INTERFACE
         SUBROUTINE ADVBC_MAP( CONCMIN, BCNAME, BCFAC, ONFILE )
            CHARACTER( 16 ), INTENT( IN )  :: CONCMIN
            CHARACTER( 16 ), INTENT( OUT ) :: BCNAME( : )
            REAL,            INTENT( OUT ) :: BCFAC ( : )
            LOGICAL,         INTENT( OUT ) :: ONFILE( : )
         END SUBROUTINE ADVBC_MAP
      END INTERFACE

#ifdef verbose_rdbcon
      integer, save :: logdev
      real sumbc
      integer s
#endif
C-----------------------------------------------------------------------
 
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
#ifdef verbose_rdbcon
         logdev = init3()
#endif
         WRITE( CONCMIN,'(1PE9.2)' ) CMIN

         N_SPC_ADV = N_GC_TRNS + N_AE_TRNS + N_NR_TRNS + N_TR_ADV + 1
                                                  ! add 1 for advecting RHOJ

 !       write ( logdev,* ) '=0= window,n_spc_adv: ', window, n_spc_adv

         ALLOCATE ( BCNAME( N_SPC_ADV ),
     &              BCFAC ( N_SPC_ADV ),
     &              ONFILE( N_SPC_ADV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating BCNAME, BCFAC, or ONFILE'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C create advected species map to bc`s
 
         CALL ADVBC_MAP ( CONCMIN, BCNAME, BCFAC, ONFILE)

		GC_ADV_FINI = N_GC_TRNS
		AE_ADV_FINI = GC_ADV_FINI + N_AE_TRNS
		NR_ADV_FINI = AE_ADV_FINI + N_NR_TRNS
		TR_ADV_FINI = NR_ADV_FINI + N_TR_ADV

		GC_ADV_STRT = 1
		AE_ADV_STRT = GC_ADV_FINI + 1
		NR_ADV_STRT = AE_ADV_FINI + 1
		TR_ADV_STRT = NR_ADV_FINI + 1
	
	
	     ALLOCATE ( ADV_MAP( N_SPC_ADV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating ADV_MAP'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         ALLOCATE ( BBUF( NBNDY,NLAYS,N_SPC_ADV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating BBUF'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         ALLOCATE ( BCON( NBNDY,NLAYS,NSPCSD ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating BBUF'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
		 
C Create global map from advected species to CGRID species

         SPC = 0
         DO VAR = 1, N_GC_TRNS
            SPC = SPC + 1
!			ADV_MAP(SPC) = 1
            ADV_MAP( SPC ) = GC_STRT - 1 + GC_TRNS_MAP( VAR )
!           write( logdev,* ) spc, adv_map( spc ), gc_spc( gc_adv_map( var ) )
         END DO
         DO VAR = 1, N_AE_TRNS
            SPC = SPC + 1
            ADV_MAP( SPC ) = AE_STRT - 1 + AE_TRNS_MAP( VAR )
!           write( logdev,* ) spc, adv_map( spc ), ae_spc( ae_adv_map( var ) )
         END DO
         DO VAR = 1, N_NR_TRNS
            SPC = SPC + 1
            ADV_MAP( SPC ) = NR_STRT - 1 + NR_TRNS_MAP( VAR )
!           write( logdev,* ) spc, adv_map( spc ), nr_spc( nr_adv_map( var ) )
         END DO
         DO VAR = 1, N_TR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = TR_STRT - 1 + TR_ADV_MAP( VAR )
         END DO
 
         ADV_MAP( N_SPC_ADV ) = ASPC	 
		 
      END IF                    ! if FIRSTIME

      MDATE  = JDATE
      MTIME  = JTIME
      STEP   = TIME2SEC( TSTEP )
      CALL NEXTIME( MDATE, MTIME, SEC2TIME( STEP / 2 ) )

      BCON = CMIN
	  BBUF = CMIN
	  
C Read & interpolate boundary concentrations

      DO SPC = GC_ADV_STRT, GC_ADV_FINI
         A2C = ADV_MAP(SPC)
         IF ( BCNAME(SPC) .NE. BLNK ) THEN
            IF ( .NOT. INTERPB( BNDY_GASC_1, BCNAME(SPC) , PNAME,
     &           MDATE, MTIME, NBNDY*NLAYS, BBUF(1, 1, SPC) ) ) THEN
               XMSG = 'Could not read ' // TRIM(BCNAME(SPC)) //
     &              ' from ' // BNDY_GASC_1
               CALL M3EXIT(PNAME, MDATE, MTIME, XMSG, XSTAT1)
            ELSE                ! found bc's (PPM) on file; convert
               DO L = 1, NLAYS
                  DO K = 1, NBNDY
                     BCON(K, L, A2C) = BCFAC(SPC) * BBUF(K, L, SPC)
                  END DO
               END DO
            END IF
         END IF
      END DO

      DO SPC = AE_ADV_STRT, AE_ADV_FINI
         A2C = ADV_MAP(SPC)
         IF ( BCNAME(SPC) .NE. BLNK ) THEN
            IF ( .NOT. INTERPB ( BNDY_AERO_1, BCNAME(SPC), PNAME,
     &           MDATE, MTIME, NBNDY*NLAYS, BBUF(1, 1, SPC) ) ) THEN
               XMSG = 'Could not read ' // TRIM(BCNAME(SPC)) //
     &              ' from ' // BNDY_AERO_1
               CALL M3EXIT(PNAME, MDATE, MTIME, XMSG, XSTAT1)
            ELSE                ! found bc's (microgram/m**3, m**2/m**3, or number/m**3)
                                ! on file; convert
               IF ( AE_TRNS(SPC-AE_ADV_STRT+1)(1:3) .EQ. 'NUM' ) THEN
                  DO L = 1, NLAYS
                     DO K = 1, NBNDY
                        BCON(K, L, A2C) = BCFAC(SPC) * BBUF(K, L, SPC)
                     END DO
                  END DO
               ELSE IF ( AE_TRNS(SPC-AE_ADV_STRT+1)(1:3) .EQ. 'SRF' ) THEN
                  DO L = 1, NLAYS
                     DO K = 1, NBNDY
                        BCON(K, L, A2C) = BCFAC(SPC) * BBUF(K, L, SPC)
                     END DO
                  END DO
               ELSE
                  DO L = 1, NLAYS
                     DO K = 1, NBNDY
                        BCON(K, L, A2C) = BCFAC(SPC) * BBUF(K, L, SPC) * KGPMG
                     END DO
                  END DO
               END IF
            END IF
         END IF
      END DO

      DO SPC = NR_ADV_STRT, NR_ADV_FINI
         A2C = ADV_MAP(SPC)
         IF ( BCNAME(SPC) .NE. BLNK ) THEN
            IF ( .NOT. INTERPB ( BNDY_NONR_1, BCNAME(SPC), PNAME,
     &           MDATE, MTIME, NBNDY*NLAYS, BBUF(1, 1, SPC) ) ) THEN
               XMSG = 'Could not read ' // TRIM(BCNAME(SPC)) //
     &              ' from ' // BNDY_NONR_1
               CALL M3EXIT(PNAME, MDATE, MTIME, XMSG, XSTAT1)
            ELSE                ! found bc's (PPM) on file; convert
               DO L = 1, NLAYS
                  DO K = 1, NBNDY
                     BCON(K, L, A2C) = BCFAC(SPC) * BBUF(K, L, SPC)
                  END DO
               END DO
            END IF
         END IF
      END DO

      DO SPC = TR_ADV_STRT, TR_ADV_FINI
         A2C = ADV_MAP(SPC)
         IF ( BCNAME(SPC) .NE. BLNK ) THEN
            IF ( .NOT. INTERPB ( BNDY_TRAC_1, BCNAME(SPC), PNAME,
     &           MDATE, MTIME, NBNDY*NLAYS, BBUF(1, 1, SPC) ) ) THEN
               XMSG = 'Could not read ' // TRIM(BCNAME(SPC)) //
     &              ' from ' // BNDY_TRAC_1
               CALL M3EXIT(PNAME, MDATE, MTIME, XMSG, XSTAT1)
            ELSE                ! found bc's (PPM) on file; convert
               DO L = 1, NLAYS
                  DO K = 1, NBNDY
                     BCON(K, L, A2C) = BCFAC(SPC) * BBUF(K, L, SPC)
                  END DO
               END DO
            END IF
         END IF
      END DO

C-----------------------------------------------------------------------
C need changes AP
C modified	  
	  
C-----------------------------------------------------------------------AP	  
C Update Mix for packets in the boundary cells; when a boundary cell is empty,
C create new packet(s) in the cell
      DO L = 1, NLAYS 
         K = 0

C Southern boundary
         R = 0
         DO C = 1, MY_NCOLS + 1
            K = K + 1
            CALL UPDATE_BOUNDARY_CELL(JDATE, JTIME, C, R, L, K, BCON)
         END DO

C Eastern boundary
         C = MY_NCOLS + 1
         DO R = 1, MY_NROWS + 1
            K = K + 1
            CALL UPDATE_BOUNDARY_CELL(JDATE, JTIME, C, R, L, K, BCON)
         END DO

C Northern boundary
         R = MY_NROWS + 1
         DO C = 0, MY_NCOLS
            K = K + 1
            CALL UPDATE_BOUNDARY_CELL(JDATE, JTIME, C, R, L, K, BCON)
         END DO

C Western boundary
         C = 0
         DO R = 0, MY_NROWS
            K = K + 1
            CALL UPDATE_BOUNDARY_CELL(JDATE, JTIME, C, R, L, K, BCON)
         END DO
       END DO	  
	  
	  
      RETURN
      END
