!--------------------------------------------------------------------!
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
C $Header: /project/yoj/arc/CCTM/src/vdiff/acm2/vdiffacm2.F,v 1.13 2012/01/19 14:37:47 yoj Exp $

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     SUBROUTINE VDIFFACMX ( DTSEC, SEDDY, DDEP, ICMP, CNGRD, DDEPJ, DDEPJ_FST )
      SUBROUTINE VDIFFACMX ( JDATE, JTIME, TSTEP, DTSEC, SEDDY, DDEP, ICMP, DDEPJ, DDEPJ_FST, CNGRD)

C-----------------------------------------------------------------------
C Asymmetric Convective Model v2 (ACM2/ACM1) -- Pleim(2006/2014)
C Function:
C   calculates vertical diffusion

C Subroutines and Functions Called:
C   INIT3, SEC2TIME, TIME2SEC, WRITE3, NEXTIME,
C   M3EXIT, EDDYX, TRI, MATRIX, PA_UPDATE_EMIS, PA_UPDATE_DDEP

C Revision History:
C   Analogous to VDIFFACM2
C 11 Apr 13 J.Young: fix double adjustment of conc for DDBF in heterogeneous HONO
C           if-then-else clauses; eliminate some white space
C 13 May 13 J.Young: access met data from VDIFF_MET module
C                    change CRANKP to THBAR, CRANKQ to THETA
C 25 May 13 J.Young: re-do the acm/eddy algorithm for computational efficiency
C 30 Apr 14 J.Young: switch THBAR and THETA
C  2 May 14 J.Pleim, J.Young: replace the banded tridiagonal matrix solver for the
C           convective PBL, with the ACM1 matrix solver followed by the tridiagonal
C           matrix solver
C   30 May 14 J.Young: split vdiff calculation out of vdiff proc.
C   07 Nov 14 J.Bash: Updated for the ASX_DATA_MOD shared data module. 
C-----------------------------------------------------------------------

      USE CGRID_SPCS          ! CGRID mechanism species
      USE GRID_CONF
      USE EMIS_DEFN          
      USE DEPV_DEFN
      USE ASX_DATA_MOD
      USE VDIFF_MAP
      USE UTILIO_DEFN
      
      USE PACKETS
      
      USE BIDI_MOD
      USE LSM_MOD, ONLY: N_LUFRAC
      USE VDIFF_DIAG, NLPCR => NLPCR_MEAN

      IMPLICIT NONE

      INCLUDE SUBST_FILES_ID  ! file name parameters

      CHARACTER( 120 ) :: XMSG = ' '

C Arguments:
      REAL, INTENT( IN )    :: DTSEC                ! model time step in seconds
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS	  
C      REAL, POINTER :: CGRID( :,:,:,: ) ! cell (average) mixing ratios
	  INTEGER    TSTEP(3) ! time step vector (HHMMSS)
      
C--- SEDDY is strictly an input, but it gets modified here
      REAL, INTENT( INOUT ) :: SEDDY    ( :,:,: )   ! flipped EDDYV
      REAL, INTENT( INOUT ) :: DDEP     ( :,:,: )   ! ddep accumulator
      REAL, INTENT( INOUT ) :: ICMP     ( :,:,: )   ! component flux accumlator 
      REAL, INTENT( INOUT ), OPTIONAL :: DDEPJ    ( :,:,:,: ) ! ddep for mosaic
      REAL, INTENT( INOUT ), OPTIONAL :: DDEPJ_FST( :,:,:,: ) ! ddep for stomtal/cuticular pathway
      REAL, INTENT( INOUT ) :: CNGRD    ( :,:,:,: ) ! cgrid replacement

C Parameters:

C explicit, THETA = 0, implicit, THETA = 1     ! Crank-Nicholson: THETA = 0.5
      REAL, PARAMETER :: THETA = 0.5,
     &                   THBAR = 1.0 - THETA

      REAL, PARAMETER :: EPS = 1.0E-06

C External Functions: None

C Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'VDIFFACMX'

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      REAL, ALLOCATABLE, SAVE :: DD_FAC     ( : )   ! combined subexpression
      REAL, ALLOCATABLE, SAVE :: DD_FACJ    ( :,: ) ! combined subexpression for mosaic
      REAL, ALLOCATABLE, SAVE :: DD_FACJ_FST( :,: ) ! combined subexpression for mosaic
      REAL, ALLOCATABLE, SAVE :: DDBF       ( : )   ! secondary DDEP
      REAL, ALLOCATABLE, SAVE :: DDBFJ      ( :,: ) ! secondary DDEP for mosaic
      REAL, ALLOCATABLE, SAVE :: DDBFJ_FST  ( :,: ) ! secondary DDEP for mosaic
      REAl, ALLOCATABLE, SAVE :: CMPF       ( : )   ! intermediate CMP
      REAL, ALLOCATABLE, SAVE :: CONC       ( :,: ) ! secondary CGRID expression
      REAL, ALLOCATABLE, SAVE :: EMIS       ( :,: ) ! emissions subexpression
      REAL        DTDENS1                       ! DT * layer 1 air density
	  
!	  REAL            XI ( N_SPC_DIFF,NLAYS )

C For TG diffusion
      REAL, POINTER :: PMIX(:)  ! packet mixing ratios (alias for ITR%P%Mix)
      TYPE(PacketListMember), POINTER :: ITR
!      REAL           CONC( N_SPC_DIFF,NLAYS ) ! secondary CGRID expression
      REAL            LI ( NLAYS ) ! see under Function, above
      REAL            DI ( NLAYS ) ! see under Function, above
      REAL            UI ( NLAYS ) ! see under Function, above
      REAL            BI ( N_SPC_DIFF,NLAYS ) ! see under Function, above
      REAL            XI ( N_SPC_DIFF,NLAYS ) ! see under Function, above
      REAL            PK ( NLAYS ) ! UI subfactor
      REAL            QK ( NLAYS ) ! LI subfactor
      INTEGER      NSTEPS( NCOLS,NROWS ) ! diffusion time steps
!      REAL         DT    ( NCOLS,NROWS ) ! eddy diff. delta T
      REAL         DELT         ! DT
	  REAL		AVG_MAX	  

      INTEGER      N   ! loop induction variables	  

C For TG drydeposition
      REAL :: ZFF(1:NCOLS, 1:NROWS) ! layer 1 surface ht above ground (m)
	  REAL :: DENS(1:NCOLS, 1:NROWS) ! air density in layer 1 (kg/m**3)
	  REAL :: DT                ! sync time step (sec)
      REAL :: ZFDENS            ! ZF * DENS
      REAL :: Q                 ! mixing ratio at beginning of time step
      REAL :: QBAR              ! mixing ratio at end of time step
      INTEGER :: D2V            ! global DEPV index from global DDEP index	  
	  INTEGER, SAVE :: WSTEP = 0 ! local write counter
	  REAL :: DDEP_CELL(1:N_SPC_DEPV)
	  
      REAL :: QC_MS(1:N_SPC_DIFF, 1:NLAYS) ! predicted mid-step average mixing ratio in cell
      REAL :: QP_MS             ! predicted mid-step packet mixing ratio
	  REAL          RHOJ  ( 0:NCOLS+1,0:NROWS+1,NLAYS ) ! density X Jacobian

C ACM Local Variables
      REAL     :: EDDY  ( NLAYS )               ! local converted eddyv
      REAL        MEDDY                         ! ACM2 intermediate var
      REAL        MBAR                          ! ACM2 mixing rate (S-1)
      REAL     :: MBARKS( NLAYS )               ! by layer
      REAL     :: MDWN  ( NLAYS )               ! ACM down mix rate
      REAL     :: MFAC  ( NLAYS )               ! intermediate loop factor
      REAL     :: AA    ( NLAYS )               ! matrix column one
      REAL     :: BB1   ( NLAYS )               ! diagonal for MATRIX1
      REAL     :: BB2   ( NLAYS )               ! diagonal for TRI
      REAL     :: CC    ( NLAYS )               ! subdiagonal
      REAL     :: EE1   ( NLAYS )               ! superdiagonal for MATRIX1
      REAL     :: EE2   ( NLAYS )               ! superdiagonal for TRI
      REAL, ALLOCATABLE, SAVE :: DD ( :,: )     ! R.H.S
      REAL, ALLOCATABLE, SAVE :: UU ( :,: )     ! returned solution
      REAL        DFACP, DFACQ
      REAL     :: DFSP( NLAYS ), DFSQ( NLAYS )  ! intermediate loop factors
      REAL        DELC, DELP, RP, RQ
      REAL     :: LFAC1( NLAYS )                ! intermediate factor for CONVT
      REAL     :: LFAC2( NLAYS )                ! intermediate factor for CONVT
      REAL     :: LFAC3( NLAYS )                ! intermediate factor for eddy
      REAL     :: LFAC4( NLAYS )                ! intermediate factor for eddy
      REAL, ALLOCATABLE, SAVE :: DEPVCR     ( : )   ! dep vel in one cell
      REAL, ALLOCATABLE, SAVE :: DEPVJCR    ( :,: ) ! dep vel in one cell for each landuse
                                                    ! category
      REAL, ALLOCATABLE, SAVE :: DEPVJCR_FST( :,: ) ! dep vel (stomatal uptake only) in
                                                    ! one cell for each landuse category
      REAL, ALLOCATABLE, SAVE :: EFAC1 ( : )
      REAL, ALLOCATABLE, SAVE :: EFAC2 ( : )
      REAL, ALLOCATABLE, SAVE :: POL   ( : )    ! prodn/lossrate = PLDV/DEPV
      REAL        PLDV_HONO                     ! PLDV for HONO
      REAL        DEPV_NO2                      ! dep vel of NO2
      REAL        DEPV_HNO3                     ! dep vel of HNO3
      REAL        FNL                           ! ACM2 Variable
      INTEGER     NLP, NL, LCBL
      INTEGER, SAVE :: NO2_HIT, HONO_HIT, HNO3_HIT, NO2_MAP, HNO3_MAP
      INTEGER, SAVE :: NH3_HIT
      REAL        DTLIM, DTS, DTACM, RZ

      INTEGER, SAVE :: LOGDEV
      INTEGER     ASTAT
      INTEGER     C, R, L, S, V, I, J, D2C, E2C         ! loop induction variables
      INTEGER     MDATE, MTIME, MSTEP                  ! internal simulation date&time
      INTEGER       FDATE                       ! interpolation date
      INTEGER       FTIME                       ! interpolation time	  

      LOGICAL, SAVE :: XMOSAIC = .FALSE., XFST = .FALSE.

      INTERFACE
         SUBROUTINE RHO_J ( JDATE, JTIME, TSTEP, RHOJ )
            INTEGER, INTENT( IN )  :: JDATE, JTIME, TSTEP( 3 )
            REAL,    INTENT( OUT ) :: RHOJ( :,:,: )
         END SUBROUTINE RHO_J	  
         SUBROUTINE MATRIX1 ( KL, A, B, E, D, X )
            INTEGER,        INTENT( IN )  :: KL
            REAL,           INTENT( IN )  :: A( : ), B( : ), E( : )
            REAL,           INTENT( IN )  :: D( :,: )
            REAL,           INTENT( OUT ) :: X( :,: )
         END SUBROUTINE MATRIX1
         SUBROUTINE TRI ( L, D, U, B, X )
            REAL,           INTENT( IN )  :: L( : ), D( : ), U( : )
            REAL,           INTENT( IN )  :: B( :,: )
            REAL,           INTENT( OUT ) :: X( :,: )
         END SUBROUTINE TRI
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

         MDATE = 0; MTIME = 0

C set auxiliary depv arrays

 
		 
         ALLOCATE ( DD_FAC( N_SPC_DEPV  ),
     &              DDBF  ( N_SPC_DEPV ),
     &              DEPVCR( N_SPC_DEPV ),
     &              EFAC1 ( N_SPC_DEPV ),
     &              EFAC2 ( N_SPC_DEPV ),
     &              POL   ( N_SPC_DEPV ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DD_FAC, DDBF, DEPVCR, EFAC1, EFAC2, or POL'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         IF ( PRESENT ( DDEPJ ) ) XMOSAIC = .TRUE.
         IF ( PRESENT ( DDEPJ_FST ) ) XFST = .TRUE.

         IF ( XMOSAIC ) THEN
            ALLOCATE ( DD_FACJ( N_LUFRAC,N_SPC_DEPV ),
     &                 DDBFJ  ( N_LUFRAC,N_SPC_DEPV ),
     &                 DEPVJCR( N_LUFRAC,N_SPC_DEPV ) )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating DD_FACJ, DDBFJ or DEPVJCR'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
            IF ( XFST ) THEN
               ALLOCATE ( DD_FACJ_FST( N_LUFRAC,N_SPC_DEPV ),
     &                    DDBFJ_FST  ( N_LUFRAC,N_SPC_DEPV ),
     &                    DEPVJCR_FST( N_LUFRAC,N_SPC_DEPV ) )
               IF ( ASTAT .NE. 0 ) THEN
                  XMSG = 'Failure allocating DD_FACJ_FST, DDBFJ_FST or DEPVJCR_FST'
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
            END IF   ! if Fst
         END IF   ! if Mosaic

         ALLOCATE ( CMPF( LCMP ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CMPF'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE ( CONC( N_SPC_DIFF,NLAYS ),
     &              EMIS( N_SPC_DIFF,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CONC or EMIS'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF
         CONC = 0.0; EMIS = 0.0   ! array assignment

         ALLOCATE ( DD( N_SPC_DIFF,NLAYS ),
     &              UU( N_SPC_DIFF,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DD or UU'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF
         DD = 0.0; UU = 0.0   ! array assignment

         HONO_HIT = 0; HNO3_HIT = 0; NO2_HIT  = 0; NH3_HIT = 0
                       HNO3_MAP = 0; NO2_MAP  = 0
         DO V = 1, N_SPC_DEPV
            IF ( DV2DF_SPC( V ) .EQ. 'NO2' ) THEN
               NO2_HIT = V
               NO2_MAP = DV2DF( V )
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'HONO' ) THEN
               HONO_HIT = V
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'HNO3' ) THEN
               HNO3_HIT = V
               HNO3_MAP = DV2DF( V )
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'NH3' ) THEN
               NH3_HIT = V
            END IF
         END DO

 
      END IF   !  if Firstime

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::AP
      
!      CALL UPDATE_MIX_ARRAY(AVG_MIX_1, JDATE, JTIME)
      
      MDATE = JDATE
      MTIME = JTIME
      FDATE = JDATE
      FTIME = JTIME	  
!	  TSTEP (3) = DTSEC
	  CALL RHO_J ( FDATE, FTIME, TSTEP, RHOJ )
!      MSTEP = TIME2SEC( TSTEP( 2 ) )
!      DTSEC = FLOAT( MSTEP )
!      CALL NEXTIME ( MDATE, MTIME, SEC2TIME( MSTEP / 2 ) )

!      CALL CONV_CGRID(AVG_MIX, MDATE, MTIME, CNGRD)

C ------------------------------------------- Row, Col LOOPS -----------
		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'BEFORE VDIFF EMISSION', AVG_MAX
      DO R = 1, NROWS
         DO C = 1, NCOLS

C for ACM time step
         DTLIM = DTSEC

C dt = .75 dzf*dzh / Kz
         DO L = 1, NLAYS - 1
            DTLIM = MIN( DTLIM, 0.75 / ( SEDDY( L,C,R ) * Grid_Data%RDX3F( L ) ) )
         END DO
         MBARKS = 0.0   ! array assignment
         MDWN = 0.0     ! array assignment

C conjoin ACM & EDDY ---------------------------------------------------

         MBAR = 0.0
         FNL = 0.0

         IF ( Met_Data%CONVCT( C,R ) ) THEN   ! Do ACM for this column
            LCBL = Met_Data%LPBL( C,R )
            MEDDY = SEDDY( 1,C,R ) / ( Met_Data%XPBL( C,R ) - X3FACE_GD( 1 ) )
            FNL = 1.0 / ( 1.0 + ( ( KARMAN / ( -Met_Data%HOL( C,R ) ) ) ** 0.3333 )
     &                / ( 0.72 * KARMAN ) )
            MBAR = MEDDY * FNL
            IF ( MBAR .LT. EPS ) THEN
               WRITE( LOGDEV,* ) ' SEDDY, MBAR, FNL, HOL = ',
     &                             SEDDY( 1,C,R ), MBAR, FNL, Met_Data%HOL( C,R )
               XMSG = '*** ACM fails ***'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
            END IF

            IF ( ( FNL .LE. 0.0 ) .OR.            ! never gonna happen for CONVCT
     &           ( LCBL .GE. NLAYS-1 ) .OR.       ! .GT. never gonna happen
     &           ( Met_Data%HOL( C,R ) .GT. -0.00001 ) )   ! never gonna happen
     &         WRITE( LOGDEV,1015 ) LCBL, MBAR, FNL, SEDDY( 1,C,R ), Met_Data%HOL( C,R )
1015           FORMAT( ' LCBL, MBAR, FNL, SEDDY1, HOL:', I3, 1X, 4(1PE13.5) )

            DO L = 1, LCBL - 1
               SEDDY( L,C,R ) = ( 1.0 - FNL ) * SEDDY( L,C,R  )
               MBARKS( L ) = MBAR
               MDWN( L ) = MBAR * ( Met_Data%XPBL( C,R ) - X3FACE_GD( L-1 ) ) * Grid_Data%RDX3F( L )
            END DO
            MBARKS( LCBL ) = MBAR * ( Met_Data%XPBL( C,R ) - X3FACE_GD( LCBL-1 ) ) * Grid_Data%RDX3F( LCBL )
            MDWN( LCBL ) = MBARKS( LCBL )

C Modify Timestep for ACM
            RZ     = ( X3FACE_GD( LCBL ) - X3FACE_GD( 1 ) ) * Grid_Data%RDX3F( 1 )
            DTACM  = 1.0 / ( MBAR * RZ )
            DTLIM  = MIN( 0.75 * DTACM, DTLIM )
         ELSE
            LCBL = 1
         END IF

C-----------------------------------------------------------------------

         NLP = INT( DTSEC / DTLIM + 0.99 )
         IF ( VDIFFDIAG ) NLPCR( C,R ) = REAL( NLP )
         DTS = DTSEC / REAL( NLP )
         DTDENS1 = DTS * Met_Data%DENS1( C,R )

         DFACP = THETA * DTS
         DFACQ = THBAR * DTS

!		PRINT *, 'DFACP   DFACQ', DFACP,DFACQ

#ifdef Verbose_Vdiff
        IF ( R .EQ. MY_NROWS / 2 .AND. C .EQ. MY_NCOLS / 2 )
     &      WRITE( LOGDEV,1021 ) Met_Data%CONVCT( C,R ), DTS, EDDYV( C,R,1 ), MBAR, FNL
1021  FORMAT( ' CONVCT, DTS, EDDYV, MBAR, FNL: ', L3, 1X, 4(1PE13.5) )
#endif

         EMIS = 0.0      ! array assignment
         DO L = 1, EMLAYS
            DO V = 1, N_SPC_DIFF
               EMIS( V,L ) = DTS * VDEMIS( DF2EM( V ),L,C,R )
            END DO
         END DO
			 
C-----------------------------------------------------------------------
            PK(1) = DFACQ * SEDDY(1, C, R) * Grid_Data%RDX3F(1)
            UI(1) = - DFACP * THETA * PK(1)
            LI(1) = 0.0
            DI(1) = 1.0 - UI(1)
            DO L = 2, NLAYS
!               PK(L) = SEDDY(L, C, R) * Grid_Data%RDX3F(L)
!               QK(L) = SEDDY(L-1, C, R) * Grid_Data%RDX3F(L)
!               UI(L) = -THETA * PK(L)
!               LI(L) = -THETA * QK(L)
!6_28 arman
               PK(L) = DFACQ * SEDDY(L, C, R) * Grid_Data%RDX3F(L)
               QK(L) = DFACQ * SEDDY(L-1, C, R) * Grid_Data%RDX3F(L)
               UI(L) = - DFACP * THETA * PK(L)
               LI(L) = - DFACP * THETA * QK(L)			   
               DI(L) = 1.0 - UI(L) - LI(L) ! UI(NLAYS) = 0
            END DO

            CALL TRIDIAG_FACTOR(NLAYS, LI, DI, UI)

            DO L = 1, NLAYS
               DO S = 1, N_SPC_DIFF
                  CONC(S, L) = CNGRD(DIFF_MAP(S), L, C, R)
               END DO
            END DO	

            DO N = 1, NLP

C Calculate the solution of the implicit Eulerian equations for vertical
C diffusion
               L = 1
               BI(:, L) = CONC(:, L) + THBAR *
     &              PK(L) * (CONC(:, L+1) - CONC(:, L))	 

               DO L = 2, NLAYS - 1
                  BI(:, L) = CONC(:, L) + THBAR *
     &                 ( PK(L) * (CONC(:, L+1) - CONC(:, L))
     &                 - QK(L) * (CONC(:, L) - CONC(:, L-1)) )
               END DO

               L = NLAYS
               BI(:, L) = CONC(:, L) - THBAR *
     &              QK(L) * (CONC(:, L) - CONC(:, L-1))

               DO S = 1, N_SPC_DIFF
                  CALL TRIDIAG_SOLVE(NLAYS, LI, DI, UI, BI(S, :), XI(S, :))
               END DO

C Predictor-corrector solution for packet compositions at next time level

C Predict mid-step average mixing ratio in cells
               QC_MS = 0.5 * CONC + 0.5 * XI

C Predict mid-step packet mixing ratio and use it to update the packet mixing
C ratio at next time level
               L = 1
               IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                  ITR => PacketsInCell(C, R, L)%Head
                  DO
                     IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                     PMIX => ITR%P%Mix
                     DO S = 1, N_SPC_DIFF
                        D2C = DIFF_MAP(S)
                        QP_MS  = PMIX(D2C) + 0.5 *
     &                        (CONC(S, L+1) - PMIX(D2C))	 
                        PMIX(D2C) = PMIX(D2C) +  (QC_MS(S, L+1) - QP_MS)
     &                       + EMIS( S,L)						

                     END DO
                     ITR => ITR%Next
                  END DO
               END IF

               DO L = 2, NLAYS - 1
                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                     ITR => PacketsInCell(C, R, L)%Head
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                        PMIX => ITR%P%Mix
                        DO S = 1, N_SPC_DIFF
                           D2C = DIFF_MAP(S)
                           QP_MS  = PMIX(D2C) + 0.5 *
     &                          (  (CONC(S, L+1) - PMIX(D2C))
     &                          -  (PMIX(D2C) - CONC(S, L-1)) )	 
                           PMIX(D2C) = PMIX(D2C) +
     &                          (  (QC_MS(S, L+1) - QP_MS)
     &                          -  (QP_MS - QC_MS(S, L-1)) )
     &                       + EMIS( S,L)	 
                        END DO
                        ITR => ITR%Next
                     END DO
                  END IF
               END DO

               L = NLAYS
               IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                  ITR => PacketsInCell(C, R, L)%Head
                  DO
                     IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                     PMIX => ITR%P%Mix
                     DO S = 1, N_SPC_DIFF
                        D2C = DIFF_MAP(S)
                        QP_MS  = PMIX(D2C) - 0.5 *
     &                       (PMIX(D2C) - CONC(S, L-1))
                        PMIX(D2C) = PMIX(D2C) - (QP_MS - QC_MS(S, L-1))
                     END DO
                     ITR => ITR%Next
                  END DO
               END IF

C Update the average mixing ratio in cells with the solution of the implicit
C equations to prepare for the next time step
               CONC = XI

            END DO              ! end time steps loop

         END DO                 ! end loop on col C
      END DO                    ! end loop on row R	
C-----------------------------------------------------------------------


!C Add emissions to packet compositions
!      DO R = 1, NROWS
!         DO C = 1, NCOLS
!            DO L = 1, EMLAYS
!               IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
!                  ITR => PacketsInCell(C, R, L)%Head
!                  DO
!                     IF ( .NOT. ASSOCIATED(ITR) ) EXIT
!                     DO V = 1, N_SPC_DIFF
!                        E2C = DF2EM( V )
!                        ITR%P%Mix(E2C) =
!     &                  ITR%P%Mix(E2C) + EMIS( V,L)
!                     END DO
!                     ITR => ITR%Next
!                  END DO
!               END IF
!            END DO
!         END DO
!      END DO
!      CALL UPDATE_MIX_ARRAY(AVG_MIX_1, JDATE, JTIME)
!		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
!		WRITE(LOGDEV, '(5X, A, E12.5)') 'AFTER VDIFF EMISSION', AVG_MAX	
!      CALL UPDATE_MIX_ARRAY(AVG_MIX_1, JDATE, JTIME)
!	PRINT * ,'CO CNGRDDDDDDDDDDDDDDDDDD', CNGRD(22,1,50,50)

!C Initialize accumulation
!         WSTEP = 0
!
!C Reset DDEP array for new accumulation
!      IF ( WSTEP .EQ. 0 ) DDEP = 0.0	  
!	  
!      DO R = 1, NROWS
!         DO C = 1, NCOLS
!		 ZFF(C,R) = Met_Data%ZF( 1,R,C )
!		 DENS = Met_Data%DENS1( C,R )
!            IF ( PacketsInCell(C, R, 1)%Count > 0 ) THEN
!               ZFDENS = ZFF(C, R) * DENS(C, R)
!               DDEP_CELL = 0.0
!               ITR => PacketsInCell(C, R, 1)%Head
!               DO
!                  IF ( .NOT. ASSOCIATED(ITR) ) EXIT
!                  DO S = 1, N_SPC_DEPV
!                     D2C = DV2DF(S)
!
!C Calculate exact mixing ratio at the end of the step assuming DEPV/ZF doesn't
!C change
!                     Q = ITR%P%Mix(D2C)
!                     QBAR = EXP(-DTS * DEPV(S, C, R) / ZFF(C, R)) * Q
!                     ITR%P%Mix(D2C) = QBAR
!
!C Convert deposition for this time step to output units (kg/hectare) and add
!C to accumulation for cell
!                     DDEP_CELL(S) = DDEP_CELL(S) +
!     &                    DD_CONV(S) * ZFDENS * (Q - QBAR)
!                  END DO
!                  ITR => ITR%Next
!               END DO
!               DDEP(C, R, :) = DDEP(C, R, :) +
!     &              DDEP_CELL / REAL(PacketsInCell(C, R, 1)%Count) ! average for cell
!            END IF
!         END DO
!      END DO



      RETURN
      END
