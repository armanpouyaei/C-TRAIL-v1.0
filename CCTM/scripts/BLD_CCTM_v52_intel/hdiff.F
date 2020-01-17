
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
C $Header: /project/yoj/arc/CCTM/src/hdiff/multiscale/hdiff.F,v 1.7 2011/10/21 16:11:22 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HDIFF ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Horizontal diffusion with constant eddy diffusivity - gen. coord.
C   The process time step is set equal to TSTEP(2). Boundary concentrations
C   are set using a Dirichlet (no flux) condition
      
C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.  All boundary input variables have the
C   same boundary perimeter structure with a thickness of 1
C   CGRID in ppm units or micro-g/m**3, #/m**3 for aerosols
      
C Subroutines and functions called:
C   INIT3, TIME2SEC, SEC2TIME, CGRID_MAP, NEXTIME, RHO_J, 
C   HCDIFF3D
 
C Revision history:
C   Jeff - 5 Nov 97, 1 Jan 98
C   DWB  - 1 Feb 98, use simple B/C (no conc gradient at domain boundary)

C   David Wong Sep. 1998
C     -- parallelized the code
C     -- removed the intermediate constant CRHOJ_Q and placed the answer of
C        the calculation directly into CGRID. Removed the next immediate
C        loop completely.

C   David Wong 1/19/99
C      -- add a loop_index call
C      -- change loop index ending point to avoid accessing invalid region.
C         (reason to do this is to prevent using boundary data from PINTERP,
C          which sets pseudo-boundary data to 0)
 
C   Daewon Byun 10/10/2000
C      -- generalized 3d horizontal diffusivity module
C      -- accomdates 3d hdiff values

C    15 Dec 00 J.Young: PE_COMM3 -> Dave Wong's f90 stenex COMM 
C     6 Aug 01 J.Young: Use HGRD_DEFN
C    25 Mar 04 G.Hammond: RK11/RK22 ghost cell updates moved outside main loop;
C                         use explicit boundary arrays for CGRID ghost cells;
C                         use SNL's "swap3d".
C    31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module
C    17 Dec 09 J.Young: fix CONC initializing error for sub-cycling timesteps (DO 344)
C                       reported by Talat Odman and Yongtao Hu of GA tech.
C    21 Jun 10 J.Young: convert for Namelist redesign
C    16 Feb 11 S. Roselle: replaced I/O-API include files w/UTILIO_DEFN
C    11 May 11 D.Wong: incorporated twoway model implementation
C-----------------------------------------------------------------------
      
      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID mechanism species
      USE UTILIO_DEFN
	  USE PACKETS
#ifdef parallel
      USE SE_MODULES            ! stenex (using SE_COMM_MODULE, SE_UTIL_MODULE)
      USE SWAP_SANDIA
#else
      USE NOOP_MODULES          ! stenex (using NOOP_COMM_MODULE, NOOP_UTIL_MODULE)
#endif

      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Arguments:
      
      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER, INTENT( IN ) :: JDATE      ! current model date, coded YYYYDDD
      INTEGER, INTENT( IN ) :: JTIME      ! current model time, coded HHMMSS
      INTEGER, INTENT( IN ) :: TSTEP( 3 ) ! time step vector (HHMMSS)
                                          ! TSTEP(1) = local output step
                                          ! TSTEP(2) = sciproc sync. step (chem)
                                          ! TSTEP(3) = twoway model time step w.r.t. wrf time
                                          !            step and wrf/cmaq call frequency

C External Functions: None
           
C Parameters:

C Advected species dimension

      INTEGER, SAVE :: N_SPC_DIFF
 
C File Variables:
 
      REAL			 :: CONC(0:NCOLS+1,0:NROWS+1)   ! conc working array
      REAL          RHOJ  ( 0:NCOLS+1,0:NROWS+1,NLAYS ) ! density X Jacobian
      CHARACTER( 8 ), SAVE :: COMMSTR                 ! for both CONC and RHOJ
      REAL          RK11  (   NCOLS+1,NROWS+1,NLAYS ) ! initially used as RHOJ
                              ! at x1 cell face, then finally as 11 eddy diff. factor
      REAL          RK22  (   NCOLS+1,NROWS+1,NLAYS ) ! initially used as RHOJ
                              ! at x2 cell face, then finally as 22 eddy diff. factor
      REAL          K11BAR3D ( NCOLS+1,NROWS+1,NLAYS ) ! ave. Cx11 eddy diff
      REAL          K22BAR3D ( NCOLS+1,NROWS+1,NLAYS ) ! ave. Cx22 eddy diff
      REAL          DT                          ! diffusion time step
      REAL          CRHOJ_Q                     ! intermediate, coupled conc.
	  

C For TG diffusion
      REAL, ALLOCATABLE, SAVE :: CMIX(:, :, :, :) ! working array
      REAL, POINTER :: PMIX(:)  ! packet mixing ratios (alias for ITR%P%Mix)
      TYPE(PacketListMember), POINTER :: ITR

      REAL :: DTDS2             ! DT/DS**2, DS slightly less than half a cell width
      REAL :: KBAR              ! average eddy diff. factor (not coupled)
      REAL :: SGD_FAC           ! sub-grid diffusion factor
	  INTEGER :: ALLOCSTAT
	  REAL :: AVG_MAX


C Local Variables:

      CHARACTER( 16 ) :: PNAME = 'HDIFF'
      
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      REAL          DX1                         ! dx1 (meters)
      REAL          DX2                         ! dx2 (meters)
      REAL, SAVE :: RDX1S                       ! reciprocal dx1*dx1
      REAL, SAVE :: RDX2S                       ! reciprocal dx2*dx2
      
      REAL          DTDX1S                      ! dt/dx1**2
      REAL          DTDX2S                      ! dt/dx2**2
      REAL          DTSEC                       ! model time step in seconds
      INTEGER       NSTEPS                      ! diffusion time steps
      INTEGER       STEP                        ! FIX dt
      INTEGER       FDATE                       ! interpolation date
      INTEGER       FTIME                       ! interpolation time

#ifdef parallel
      INTEGER      M
      REAL,    ALLOCATABLE, SAVE :: HALO_SOUTH( :,:,: )
      REAL,    ALLOCATABLE, SAVE :: HALO_NORTH( :,:,: )
      REAL,    ALLOCATABLE, SAVE :: HALO_WEST ( :,:,: )
      REAL,    ALLOCATABLE, SAVE :: HALO_EAST ( :,:,: )
      REAL,    ALLOCATABLE, SAVE :: BUF_NS    ( : )
      REAL,    ALLOCATABLE, SAVE :: BUF_EW    ( : )
#endif

      INTEGER, ALLOCATABLE, SAVE :: DIFF_MAP( : )   ! global diff map to CGRID

      INTEGER, SAVE :: LOGDEV

      INTEGER      C, R, L, S, V, N            ! loop counters
      INTEGER      D2C, IOS

      CHARACTER( 96 ) :: XMSG = ' '
     
      INTEGER MY_TEMP
      INTEGER, SAVE :: STARTROW, ENDROW
      INTEGER, SAVE :: STARTCOL, ENDCOL

      INTERFACE
         SUBROUTINE RHO_J ( JDATE, JTIME, TSTEP, RHOJ )
            INTEGER, INTENT( IN )  :: JDATE, JTIME, TSTEP( 3 )
            REAL,    INTENT( OUT ) :: RHOJ( :,:,: )
         END SUBROUTINE RHO_J
         SUBROUTINE HCDIFF3D ( JDATE, JTIME, K11BAR, K22BAR, DT )
            INTEGER, INTENT( IN )  :: JDATE, JTIME
            REAL,    INTENT( OUT ) :: K11BAR( :,:,: ), K22BAR( :,:,: )
            REAL,    INTENT( OUT ) :: DT
         END SUBROUTINE HCDIFF3D
      END INTERFACE
 
C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

C Get dx1 from COORD include file

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD
     &          * COS( PI180*( YORIG_GD + YCELL_GD*FLOAT( GL_NROWS/2 ))) ! in m.
            DX2 = DG2M * YCELL_GD   ! in m.
         ELSE
            DX1 = XCELL_GD          ! in m.
            DX2 = YCELL_GD          ! in m.
         END IF

         RDX1S = 1.0 / ( DX1 * DX1 )
         RDX2S = 1.0 / ( DX2 * DX2 )

         N_SPC_DIFF = N_GC_TRNS + N_AE_TRNS + N_NR_TRNS + N_TR_DIFF

#ifdef parallel
         ALLOCATE ( HALO_SOUTH( NCOLS,NLAYS,N_SPC_DIFF ),
     &              HALO_NORTH( NCOLS,NLAYS,N_SPC_DIFF ),
     &              HALO_WEST ( NROWS,NLAYS,N_SPC_DIFF ),
     &              HALO_EAST ( NROWS,NLAYS,N_SPC_DIFF ),
     &              BUF_NS    ( NCOLS*NLAYS*N_SPC_DIFF ),
     &              BUF_EW    ( NROWS*NLAYS*N_SPC_DIFF ), STAT = IOS )
         IF ( IOS .NE. 0 ) THEN
            XMSG = 'Failure allocating HALO_SOUTH, HALO_NORTH, HALO_WEST, HALO_EAST, BUF_NS, or BUF_EW'
            CALL M3EXIT ( PNAME, FDATE, FTIME, XMSG, XSTAT1 )
         END IF
#endif

         ALLOCATE ( DIFF_MAP( N_SPC_DIFF ), STAT = IOS )
         IF ( IOS .NE. 0 ) THEN
            XMSG = 'Failure allocating DIFF_MAP'
            CALL M3EXIT ( PNAME, FDATE, FTIME, XMSG, XSTAT1 )
         END IF

C Create global map to CGRID
 
         S = 0
         DO V = 1, N_GC_TRNS
            S = S + 1
            DIFF_MAP( S ) = GC_STRT - 1 + GC_TRNS_MAP( V )
         END DO
         DO V = 1, N_AE_TRNS
            S = S + 1
            DIFF_MAP( S ) = AE_STRT - 1 + AE_TRNS_MAP( V )
         END DO
         DO V = 1, N_NR_TRNS
            S = S + 1
            DIFF_MAP( S ) = NR_STRT - 1 + NR_TRNS_MAP( V )
         END DO
         DO V = 1, N_TR_DIFF
            S = S + 1
            DIFF_MAP( S ) = TR_STRT - 1 + TR_DIFF_MAP( V )
         END DO
 
C Get file start and end indices for subdomain
 
         CALL SUBST_LOOP_INDEX ( 'C', 1, NCOLS, 1, MY_TEMP, STARTCOL, ENDCOL )
         CALL SUBST_LOOP_INDEX ( 'R', 1, NROWS, 1, MY_TEMP, STARTROW, ENDROW )

         WRITE( COMMSTR,'(4I2)' )  1, 0, 2, 0
C Allocate memory for working arrays
         ALLOCATE(
!     &        CONC(0:NCOLS+1, 0:NROWS+1),
     &        CMIX(1:N_SPC_DIFF, 0:NCOLS+1, 0:NROWS+1, 1:NLAYS),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT /= 0 ) THEN
            XMSG = 'Failure allocating CONC, CMIX'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

C For no flux boundary conditions - diffusivity across boundary will be zero
         CONC(0, :) = 0.0       ! West boundary
         CONC(NCOLS+1, :) = 0.0 ! East boundary
         CONC(:, 0) = 0.0       ! South boundary
         CONC(:, NROWS+1) = 0.0 ! North boundary

         CMIX(:, 0, :, :) = 0.0 ! West boundary
         CMIX(:, NCOLS+1, :, :) = 0.0 ! East boundary
         CMIX(:, :, 0, :) = 0.0 ! South boundary
         CMIX(:, :, NROWS+1, :) = 0.0 ! North boundary
      END IF                    ! if firstime

		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'HORDIFF START CODE MAX', AVG_MAX
                                     
      DTSEC = FLOAT( TIME2SEC( TSTEP( 2 ) ) )
      FDATE = JDATE
      FTIME = JTIME
 
C Get the computational grid ( rho X Jacobian ) for this step

      CALL RHO_J ( FDATE, FTIME, TSTEP, RHOJ )

      CALL SUBST_COMM ( RHOJ, DSPL_N0_E0_S1_W1, DRCN_S_W, COMMSTR )


C initialize RK11, RK22 with face values for RHOJ (assumes dx1 = dx2)

      RK11 = 0.0   ! array assignment
      RK22 = 0.0   ! array assignment

C Do the gridded computation for horizontal diffusion

C Get the contravariant eddy diffusivities

      CALL HCDIFF3D ( FDATE, FTIME, K11BAR3D, K22BAR3D, DT )

C get number of steps based on eddy time 
 
      NSTEPS = INT ( DTSEC / DT ) + 1
      DT = DTSEC / FLOAT( NSTEPS )
 
      WRITE( LOGDEV,1005 ) DT, NSTEPS

      DTDX1S = DT * RDX1S
      DTDX2S = DT * RDX2S

C Couple eddy diffusivity with RHOJ
      DO L = 1, NLAYS
         DO R = 1, NROWS + 1
            DO C = 1, NCOLS + 1
               RK11(C, R, L) = 0.5 * (RHOJ(C, R, L) + RHOJ(C-1, R, L))
     &              * K11BAR3D(C, R, L)
               RK22(C, R, L) = 0.5 * (RHOJ(C, R, L) + RHOJ(C, R-1, L))
     &              * K22BAR3D(C, R, L)
            END DO
         END DO
      END DO

C No flux boundary conditions (need to do it this way for diffusion between a
C packet and the surrounding cells)
      RK11(1, :, :) = 0.0       ! West boundary
      RK11(NCOLS+1, :, :) = 0.0 ! East boundary
      RK22(:, 1, :) = 0.0       ! South boundary
      RK22(:, NROWS+1, :) = 0.0 ! North boundary
!dbg
!dbg      CALL UPDATE_MIX_ARRAY(MAX_MIX_1, JDATE, JTIME)
!dbg      CALL UPDATE_MIX_ARRAY(MIN_MIX_1, JDATE, JTIME)
!dbg      DO S = 1, N_SPC_DIFF
!dbg         D2C = DIFF_MAP(S)
!dbg         DO L = 1, NLAYS
!dbg            MAX_BD(L, S) = MAXVAL(MAX_MIX(:, :, L, D2C))
!dbg            MIN_BD(L, S) = MINVAL(MIN_MIX(:, :, L, D2C))
!dbg         END DO
!dbg      END DO

C Update AVG_MIX
      CALL UPDATE_MIX_ARRAY(AVG_MIX_1, JDATE, JTIME)

C Reordering greatly improves performance of packet diffusion; boundary is OK
C because all the values are zero
      DO S = 1, N_SPC_DIFF
         CMIX(S, 1:NCOLS, 1:NROWS, :) = AVG_MIX(1:NCOLS, 1:NROWS, :, DIFF_MAP(S))
      END DO

		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'HORDIFF 360 MAX', AVG_MAX

      DO N = 1, NSTEPS

C Update packet mixing ratios with explicit step of Eulerian diffusion
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
                     ITR => PacketsInCell(C, R, L)%Head
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                        PMIX => ITR%P%Mix
                        DO S = 1, N_SPC_DIFF
                           D2C = DIFF_MAP(S)
                           PMIX(D2C) = PMIX(D2C) +
     &                          ( DTDX1S *
     &                          ( RK11(C+1, R, L) * (CMIX(S, C+1, R, L) - PMIX(D2C))
     &                          - RK11(C, R, L) * (PMIX(D2C) - CMIX(S, C-1, R, L)) )
     &                          + DTDX2S *
     &                          ( RK22(C, R+1, L) * (CMIX(S, C, R+1, L) - PMIX(D2C))
     &                          - RK22(C, R, L) * (PMIX(D2C) - CMIX(S, C, R-1, L)) ) )
     &                          / RHOJ(C, R, L)
                        END DO
                        ITR => ITR%Next
                     END DO
                  END IF
               END DO
            END DO
         END DO



C Update cell mixing ratios with explicit step of Eulerian diffusion
         DO S = 1, N_SPC_DIFF
            D2C = DIFF_MAP(S)
            DO L = 1, NLAYS
               CONC(1:NCOLS, 1:NROWS) = AVG_MIX(1:NCOLS, 1:NROWS, L, D2C)
               DO R = 1, NROWS
                  DO C = 1, NCOLS
                     AVG_MIX(C, R, L, D2C) = CONC(C, R) +
     &                    ( DTDX1S *
     &                    ( RK11(C+1, R, L) * (CONC(C+1, R) - CONC(C, R))
     &                    - RK11(C, R, L) * (CONC(C, R) - CONC(C-1, R)) )
     &                    + DTDX2S *
     &                    ( RK22(C, R+1, L) * (CONC(C, R+1) - CONC(C, R))
     &                    - RK22(C, R, L) * (CONC(C, R) - CONC(C, R-1)) ) )
     &                    / RHOJ(C, R, L)
                  END DO
               END DO
            END DO
         END DO

         DO S = 1, N_SPC_DIFF
            CMIX(S, 1:NCOLS, 1:NROWS, :) = AVG_MIX(1:NCOLS, 1:NROWS, :, DIFF_MAP(S))
         END DO

!		 AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
!		 WRITE(LOGDEV, '(5X, A, E12.5)') 'HDIFF MAX', AVG_MAX
      END DO                    ! STEPS

C Sub-grid diffusion:
C 1) pairwise diffusion between packets in a cell is approximated by diffusion
C    between each packet and an imaginary packet with the cell's average mixing
C    ratios; the distance between each packet and the imaginary average packet
C    is taken to be a little less than half a cell width, i.e.,
C            PMIX = PMIX + DTDS2 * KBAR * (CMIX - PMIX)
C 2) to avoid the need for a time step limitation, x=DTDS2*KBAR is approximated
C    by 1-exp(-x), which increases from 0 to 1 as x increases from 0, so each
C    PMIX moves closer to CMIX, without crossing it
C Note that sub-grid diffusion leaves the average over the packets unchanged

		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'HORDIFF 434 MAX', AVG_MAX


      IF ( MAX_SGD_FAC > 0.0 ) THEN
         DTDS2 = 5.0 * DTSEC / MAX(DX1, DX2)**2
         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  IF ( PacketsInCell(C, R, L)%Count > 1 ) THEN
                     KBAR = 0.25 *
     &                    ( K11BAR3D(C, R, L) + K11BAR3D(C+1, R, L)
     &                    + K22BAR3D(C, R, L) + K22BAR3D(C, R+1, L) )
                     SGD_FAC = MIN(1.0 - EXP(-DTDS2 * KBAR), MAX_SGD_FAC)

                     ITR => PacketsInCell(C, R, L)%Head
                     DO
                        IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                        PMIX => ITR%P%Mix
                        DO S = 1, N_SPC_DIFF
                           D2C = DIFF_MAP(S)
                           PMIX(D2C) = PMIX(D2C) +
     &                          SGD_FAC * (CMIX(S, C, R, L) - PMIX(D2C))
                        END DO
                        ITR => ITR%Next
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END IF


!dbg
!dbg      AVG_MIX_ED = AVG_MIX
!dbg      CALL UPDATE_MIX_ARRAY(AVG_MIX_1, JDATE, JTIME)
!dbg
!dbg      DO S = 1, N_SPC_DIFF
!dbg         D2C = DIFF_MAP(S)
!dbg         MAX_AVG_ED(S) = MAX(
!dbg     &        MAXVAL(AVG_MIX(:, :, :, D2C)),
!dbg     &        MAXVAL(AVG_MIX_ED(:, :, :, D2C)),
!dbg     &        CMIN)
!dbg      END DO
!dbg
!dbg      MAX_DIFF_ED = -1.0
!dbg      DO S = 1, N_SPC_DIFF
!dbg         D2C = DIFF_MAP(S)
!dbg         DO L = 1, NLAYS
!dbg            DO R = 1, NROWS
!dbg               DO C = 1, NCOLS
!dbg                  IF ( PacketsInCell(C, R, L)%Count > 0 ) THEN
!dbg                     MAX_DIFF_ED = MAX(
!dbg     &                    ABS(AVG_MIX(C, R, L, D2C) - AVG_MIX_ED(C, R, L, D2C))
!dbg     &                    / MAX_AVG_ED(S),
!dbg     &                    MAX_DIFF_ED)
!dbg                  END IF
!dbg               END DO
!dbg            END DO
!dbg         END DO
!dbg      END DO
!dbg      WRITE(LOGDEV, '(5X, A, E12.5)') 'MAX_HDIFF_ED: ', MAX_DIFF_ED
!dbg
!dbg      CALL UPDATE_MIX_ARRAY(MAX_MIX_1, JDATE, JTIME)
!dbg      CALL UPDATE_MIX_ARRAY(MIN_MIX_1, JDATE, JTIME)
!dbg      DO S = 1, N_SPC_DIFF
!dbg         D2C = DIFF_MAP(S)
!dbg         DO L = 1, NLAYS
!dbg            MAX_AD(L, S) = MAXVAL(MAX_MIX(:, :, L, D2C))
!dbg            MIN_AD(L, S) = MINVAL(MIN_MIX(:, :, L, D2C))
!dbg         END DO
!dbg      END DO
!dbg
!dbg      MAX_VIO = MAXVAL( MAX(
!dbg     &     MAX_AD - ONE_PLUS * MAX_BD,
!dbg     &     MIN_BD - ONE_PLUS * MAX(MIN_AD, CMIN) ) )
!dbg      IF ( MAX_VIO > 0.0 )
!dbg     &     WRITE(LOGDEV, '(5X, A, E12.5)') 'HDIFF VIOLATION', MAX_VIO

		AVG_MAX = MAXVAL(AVG_MIX(:, :, :, 22))
		WRITE(LOGDEV, '(5X, A, E12.5)') 'HORDIFF END CODE MAX', AVG_MAX

      RETURN

 1005 FORMAT(/ 5X, 'H-eddy DT & integration steps: ', 1PE15.7, I8)

      END
