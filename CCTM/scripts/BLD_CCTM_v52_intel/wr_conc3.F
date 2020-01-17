      SUBROUTINE WR_CONC3 ( JDATE, JTIME, TSTEP )
!-------------------------------------------------------------------------------
! Name:     Writes concentration of a specific species of all packets
! Purpose:  Finds the lat/lon/alt of all TG packets and writes concentration of 
!           specific species into files for each timestep 
! Notes:    
! Revised:  25 Jan 2018 Arman Pouyaei  
!-------------------------------------------------------------------------------
C---------------------------------------------------------------------
      USE PACKETS
	  USE GRID_CONF             ! horizontal domain specifications
      USE UTILIO_DEFN
	  USE CGRID_SPCS
	
      IMPLICIT NONE

      INCLUDE SUBST_FILES_ID

	  
	  INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP        ! output timestep (HHMMSS)
	  INTEGER		p_num
	  TYPE(PacketListMember), POINTER :: ITR, PLM
	  TYPE(Packet), POINTER :: P
	  REAL(8) :: X, Y, Z
	  INTEGER :: ID_packet        ! packet ID
	  REAL :: CO        ! CO conc
	  INTEGER :: C, R, L, M        ! column, row and layer indices
	  REAL, ALLOCATABLE, SAVE :: CONC( : )
	  CHARACTER( 16 ) :: PNAME = 'wr_conc3'
      CHARACTER( 96 ) :: XMSG = ' '
	  character(len=70) :: filename
	  CHARACTER( 16 ) :: VNAME
	  INTEGER          ALLOCSTAT
	  REAL		:: xx, yy
	  REAL		:: lat,lon, alt
	  REAL		:: phi, lambda
	  REAL		:: phi1, phi2, phi0, lambda0
	  
	  
	  INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2
	  INTEGER :: LOGDEV
	  LOGICAL, SAVE :: FIRSTIME = .TRUE.
      REAL, ALLOCATABLE, SAVE :: ZHEI( :,:,: ) ! full-layer Jacobian 
	  REAL, ALLOCATABLE, SAVE :: ZHEI2( :,:,: ) ! full-layer Jacobian
	  INTEGER :: CU, RU, LU     ! updated column, row and layer indices
C-----------------------------------------------------------------------

         ALLOCATE ( CONC( NSPCSD ),
     &     STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'CONC memory allocation failed'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG)
         END IF
            ALLOCATE ( ZHEI( NCOLS,NROWS,NLAYS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating ZH'
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG)
            END IF
	CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                     STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3)
			
C--------------------------------------------------------------------
         VNAME = 'ZH'
		 IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                       STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                       1,NLAYS, JDATE, JTIME, ZHEI ) ) THEN
            XMSG = 'Could not interpolate ZH from MET_CRO_3D - '
     &           // 'Using INTERPX!'
            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
			END IF

	WRITE(filename,fmt='(a,i7,i6,a)') 'TG',JDATE,JTIME,'.dat'
	OPEN(100,FILE=filename)
	DO L = 1, NLAYS
         DO R = 0, NROWS + 1
            DO C = 0, NCOLS + 1
               ITR => PacketsInCell(C, R, L)%Head
               DO
                  IF ( .NOT. ASSOCIATED(ITR) ) EXIT
                  P => ITR%P
                  PLM => ITR
                  ITR => ITR%Next
		  X = P%X
		  Y = P%Y
		  Z = P%Z
		  ID_packet = P%ID_GD
		  CONC(:) = P%Mix
		  CO = CONC(22) * 1000
		  xx = X - (-33750)
		  yy = Y - (-33750)
		  phi1 = 30
		  phi2 = 60
		  lambda0 = 126
!		  phi0 = (phi1 + phi2) * 0.5
		  phi0 = 38
		  CALL xy2ll_lam (xx, yy, phi1, phi2, lambda0, phi0, phi, lambda)
		  lat = phi
		  lon = lambda
		  LU=0.0
		  RU=0.0
		  CU=0.0
		CALL GET_CELL(X, Y, Z, CU, RU, LU)
		IF ( 0 <= CU .AND. CU < NCOLS ) THEN
			IF ( 0 <= RU .AND. RU < NROWS ) THEN
				IF (LU .NE. NLAYS) THEN
				alt = 0.5 * (ZHEI(CU,RU,LU) + ZHEI(CU,RU,LU+1))
				ELSE
				alt = ZHEI(CU,RU,LU)
				END IF
			ELSE
				IF (LU .NE. NLAYS) THEN
				alt = 0.5 * (ZHEI(CU,RU-1,LU) + ZHEI(CU,RU-1,LU+1))
				ELSE
				alt = ZHEI(CU,RU-1,LU)
				END IF
			END IF
		ELSE
			IF ( 0 <= RU .AND. RU < NROWS ) THEN
				IF (LU .NE. NLAYS) THEN
				alt = 0.5 * (ZHEI(CU-1,RU,LU) + ZHEI(CU-1,RU,LU+1))
				ELSE
				alt = ZHEI(CU-1,RU,LU)
				END IF
			ELSE
				IF (LU .NE. NLAYS) THEN
				alt = 0.5 * (ZHEI(CU-1,RU-1,LU) + ZHEI(CU-1,RU-1,LU+1))
				ELSE
				alt = ZHEI(CU-1,RU-1,LU)
				END IF
			END IF	
		END IF

		  WRITE(100,'(I10,1x,I10,1x,I10,1x,F20.3,1x,
     &		F20.3,1x,F20.3,1x,F20.3)')
     &		JDATE,JTIME,ID_packet,lat,lon,alt,CO
               END DO
            END DO
         END DO
      END DO

	CLOSE(100)
	RETURN

      END SUBROUTINE WR_CONC3

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	SUBROUTINE xy2ll_lam (xx, yy, phi1, phi2, lambda0, phi0, phi, lambda)

!-------------------------------------------------------------------------------
! Name:     (X,Y) to Latitude-Longitude for Lambert Conformal Projection
! Purpose:  Calcluates latitude-longitude for a given (X,Y) pair from origin
!           and Lambert conformal projection information.
! Notes:    Equations adapted from http://mathworld.wolfram.com.
! Revised:  12 Dec 2007  Original version.  (T. Otte)
!           18 Sep 2009  Added reference latitude (PHI0) as an argument.
!                        Converted to double-precision.  Corrected comments
!                        associated with RHO0.  Corrected calculation of PSI
!                        (with no impact on results).  (T. Otte)
!           02 Oct 2009  Changed algorithm to follow Wolfram to eliminate the
!                        divide-by-zero condition for computing latitude along
!                        the standard longitude.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!-------------------------------------------------------------------------------

	IMPLICIT NONE

	REAL(8)                      :: deg2rad    ! convert degrees to radians
	REAL(8)                      :: drearth    ! earth radius [m]
	REAL,          INTENT(OUT)   :: lambda     ! longitude [deg]
	REAL(8)                      :: lambdarad  ! longitude [rad]
	REAL,          INTENT(IN)    :: lambda0    ! standard longitude [deg]
	REAL(8)                      :: lambda0rad ! standard longitude [rad]
	REAL,          INTENT(OUT)   :: phi        ! latitude [deg]
	REAL(8)                      :: phirad     ! latitude [rad]
	REAL,          INTENT(IN)    :: phi0       ! reference latitude [deg]
	REAL(8)                      :: phi0rad    ! reference latitude [rad]
	REAL,          INTENT(IN)    :: phi1       ! true latitude 1 [deg]
	REAL(8)                      :: phi1rad    ! true latitude 1 [rad]
	REAL,          INTENT(IN)    :: phi2       ! true latitude 2 [deg]
	REAL(8)                      :: phi2rad    ! true latitude 2 [rad]
	REAL(8)                      :: pi
	REAL(8)                      :: piover2    ! pi/2
	REAL(8)                      :: piover4    ! pi/4
	REAL(8)                      :: psi        ! auxiliary function
	REAL(8)                      :: rad2deg
	REAL(8)                      :: rho
	REAL(8)                      :: rho0       ! polar radius to origin
	REAL(8)                      :: term0
	REAL(8)                      :: term1
	REAL(8)                      :: term2
	REAL(8)                      :: theta      ! polar angle
	REAL(8)                      :: sinphi0    ! cone constant
	REAL(8)                      :: sinphi0inv ! 1/sinphi0
	REAL(8)					   :: rearth
	REAL,          INTENT(IN)    :: xx         ! X-coordinate from origin
	REAL,          INTENT(IN)    :: yy         ! Y-coordinate from origin

!-------------------------------------------------------------------------------
! Compute constants.
!-------------------------------------------------------------------------------

	piover4 = DATAN(1.0d0)
	pi      = 4.0d0 * piover4
	piover2 = 2.0d0 * piover4
	deg2rad = pi / 1.8d2
	rad2deg = 1.8d2 / pi
	rearth  = 6371000
	drearth = DBLE(rearth)

!-------------------------------------------------------------------------------
! Compute cone constant, SINPHI0.
!-------------------------------------------------------------------------------

	phi0rad = DBLE(phi0) * deg2rad  ! convert PHI0 from degrees to radians
	phi1rad = DBLE(phi1) * deg2rad  ! convert PHI1 from degrees to radians
	phi2rad = DBLE(phi2) * deg2rad  ! convert PHI2 from degrees to radians
	term0 = DTAN ( piover4 + phi0rad/2.0d0 )
	term1 = DTAN ( piover4 + phi1rad/2.0d0 )
	term2 = DTAN ( piover4 + phi2rad/2.0d0 )

	sinphi0 = DLOG ( DCOS(phi1rad) / DCOS(phi2rad) )
	sinphi0 = sinphi0 / DLOG (term2 / term1)
	sinphi0inv = 1.0d0 / sinphi0

!-------------------------------------------------------------------------------
! Compute polar radius to origin, RHO0, where origin is at PHI0.
!-------------------------------------------------------------------------------

	psi  = drearth * DCOS(phi1rad) * sinphi0inv * (term1**sinphi0)
	rho0 = psi / (term0**sinphi0)

!-------------------------------------------------------------------------------
! Compute longitude, LAMBDA.
!-------------------------------------------------------------------------------

	lambda0rad = lambda0 * deg2rad
	theta     = DATAN( DBLE(xx) / (rho0 - DBLE(yy)) )
	lambdarad = lambda0rad + theta * sinphi0inv
	lambda    = REAL(lambdarad * rad2deg)

!-------------------------------------------------------------------------------
! Compute latitude, PHI.
!-------------------------------------------------------------------------------

	rho = DSQRT( DBLE(xx)*DBLE(xx) + (rho0-DBLE(yy))*(rho0-DBLE(yy)) )
	rho = DSIGN(1.0d0, sinphi0) * rho

	phirad = (psi / rho)**sinphi0inv
	phirad = 2.0d0 * DATAN(phirad) - piover2
	phi    = REAL(phirad * rad2deg)

	END SUBROUTINE xy2ll_lam
