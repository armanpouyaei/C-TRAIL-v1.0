
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
C $Header: /project/yoj/arc/CCTM/src/init/yamo/opconc.F,v 1.4 2011/10/21 16:11:24 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPCONC (MIX_FILE_1, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Create the IO/API netCDF header and open the output CONC file

C Revision history:
C   Jeff - Dec 00 - split out from initscen.F
C                 - move CGRID_MAP into f90 module
C   Jeff - Feb 01 - assumed shape arrays
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN
C
C    3 Sep 01 David Wong
C     -- let PE 0 open CTM_CONC_1 as new and later on let the rest open
C        it for read and write
C     -- put an explicit barrier before opening a new netCDF file to avoid
C        NCOPEN error
C    7 May 03 J.Young: open and close conc file in processor 0; spin-wait to
C                      deal with nfs network latency for mpich cluster
C   28 Aug 03 J.Young: following Zion Wang at CERT, remove the spin-wait and
C                      simplify opening and closing CONC file (see initscen)
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   30 May 05 J.Young: optional save derived vert. vel. to conc file
C   21 Jan 06 J.Young: add subset non-reactives capability with NR_CONC.EXT
C                      and subset layers saved to conc file; no longer save
C                      RHOJ to conc file
C   22 Aug 06 J.Young: enable capturing IC's as 0 timestep
C   11 Jun 08 J.Young: clean up
C   30 Mar 10 David Wong: removed unnecessary barrier
C   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN;
C                      removed deprecated TRIMLEN
C   Aug 15 D.Wong:     Replaced MYPE with IO_PE_INCLUSIVE for parallel I/O
C                        implementation
C                      Added a block of code to let non I/O processor to open
C                        CTM_CONC_1 file
C   Dec 15 D.Wong:     Removed all the data writing parts and moved it into a
C                      new file, wr_init.F
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE WVEL_DEFN             ! derived vertical velocity component
      USE STD_CONC              ! standard CONC
      USE UTILIO_DEFN
#ifdef parallel
      USE SE_MODULES            ! stenex (using SE_UTIL_MODULE)
#else
      USE NOOP_MODULES          ! stenex (using NOOP_UTIL_MODULE)
#endif
      USE CGRID_SPCS
      USE PACKETS

      IMPLICIT NONE

      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      INTEGER      JDATE        ! starting date (YYYYDDD)
      INTEGER      JTIME        ! starting time (HHMMSS)
      INTEGER      TSTEP        ! output timestep (HHMMSS)
      CHARACTER(16), INTENT(IN) :: MIX_FILE_1

C Local Variables:

      REAL, POINTER :: MGRID( :,:,:,: ) ! alias for a gridded mixing ratio array
      REAL, ALLOCATABLE :: DBUFF( :,:,: )
      INTEGER      ALLOCSTAT
      CHARACTER( 16 ) :: PNAME = 'OPCONC'
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 28 ) :: SSTR = ' species saved to CONC file:'

C environment variable for no. of layers from bottom to save on CONC file
      CHARACTER( 16 ) :: NLAYS_CONC = 'NLAYS_CONC'

C  environment variable description
      CHARACTER( 80 ) :: VARDESC

      INTEGER      LOGDEV

      INTEGER      K, KD, L, SPC, V   ! loop counters
      INTEGER      STRT, FINI         ! loop counters
      INTEGER      INDX

      LOGICAL, EXTERNAL :: FLUSH3

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

C Make MGRID an alias for the gridded mixing ratio array corresponding to
C the file
      SELECT CASE ( MIX_FILE_1 )

         CASE ( CLS_MIX_1 )
            MGRID => CLS_MIX

         CASE ( AVG_MIX_1 )
            MGRID => AVG_MIX

         CASE ( MIN_MIX_1 )
            MGRID => MIN_MIX

         CASE ( MAX_MIX_1 )
            MGRID => MAX_MIX

         CASE ( OLD_MIX_1 )
            MGRID => OLD_MIX

         CASE DEFAULT
            XMSG =  '"' // MIX_FILE_1 // '"' //
     &           ' is not a valid file selection'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)

      END SELECT
      
C Set output file characteristics based on COORD.EXT and open it

      FTYPE3D = GRDDED3
      SDATE3D = JDATE
      STIME3D = JTIME
      TSTEP3D = TSTEP
      NVARS3D = N_CSPCS
      IF ( W_VEL ) NVARS3D = NVARS3D + 1   ! for W_YAMO
      NCOLS3D = GL_NCOLS
      NROWS3D = GL_NROWS
      NLAYS3D = C_NLAYS
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
!     VGTPUN3D = VGTPUN_GD ! currently, not defined
      DO L = 1, NLAYS3D + 1
         VGLVS3D( L ) = VGLVS_GD( L )
      END DO
      GDNAM3D = GRID_NAME  ! from HGRD_DEFN

      FDESC3D( 1 ) = 'Concentration file output'
      FDESC3D( 2 ) = 'From CMAQ model dyn alloc version CTM'
      FDESC3D( 3 ) = 'Set of variables (possibly) reduced from CGRID'
      FDESC3D( 4 ) = 'For next scenario continuation runs,'
      FDESC3D( 5 ) = 'use the "one-step" CGRID file'
      KD = 5
      V = CONC_BLEV
      L = 0
      DO K = KD + 1, MIN ( C_NLAYS + KD, MXDESC3 )
         L = L + 1
         WRITE( FDESC3D( K ),'( "Layer", I3, " to", I3, " " )' )
     &   V + L - 1, L
      END DO
      IF ( ( KD + 1 + L ) .LT. MXDESC3 ) THEN
         DO K = KD + 1 + L, MXDESC3
            FDESC3D( K ) = ' '
         END DO
      END IF

      WRITE( LOGDEV,* ) ' '
      WRITE( LOGDEV,* ) '      Conc File Header Description:'
      DO K = 1, KD + L
         WRITE( LOGDEV,* ) '    => ', TRIM( FDESC3D( K ) )
      END DO

      V = 0
      STRT = 1
      FINI = N_C_GC_SPC
      DO SPC = STRT, FINI
         V = V + 1
         VTYPE3D( SPC ) = M3REAL
         VNAME3D( SPC ) = C_GC_SPC( V )
         UNITS3D( SPC ) = 'ppmV'
         VDESC3D( SPC ) = 'Variable ' // VNAME3D( SPC )
      END DO

      V = 0
      STRT = FINI + 1
      FINI = FINI + N_C_AE_SPC
      DO SPC = STRT, FINI
         V = V + 1
         VTYPE3D( SPC ) = M3REAL
         VNAME3D( SPC ) = C_AE_SPC( V )   ! from STD_CONC module
         IF ( VNAME3D( SPC )(1:3) .EQ. 'NUM' ) THEN
            UNITS3D( SPC ) = 'number/m**3'
         ELSE IF ( VNAME3D( SPC )(1:3) .EQ. 'SRF' ) THEN
            UNITS3D( SPC ) = 'm**2/m**3'
         ELSE
            UNITS3D( SPC ) = 'micrograms/m**3'
         END IF
         VDESC3D( SPC ) = 'Variable ' // VNAME3D( SPC )
      END DO

      V = 0
      STRT = FINI + 1
      FINI = FINI + N_C_NR_SPC
      DO SPC = STRT, FINI
         V = V + 1
         VTYPE3D( SPC ) = M3REAL
         VNAME3D( SPC ) = C_NR_SPC( V )   ! from STD_CONC module
         UNITS3D( SPC ) = 'ppmV'
         VDESC3D( SPC ) = 'Variable ' // VNAME3D( SPC )
      END DO

      V = 0
      STRT = FINI + 1
      FINI = FINI + N_C_TR_SPC ! write all TR species
      DO SPC = STRT, FINI
         V = V + 1
         VTYPE3D( SPC ) = M3REAL
         VNAME3D( SPC ) = C_TR_SPC( V )   ! from STD_CONC module
         UNITS3D( SPC ) = 'ppmV'
         VDESC3D( SPC ) = 'Variable ' // VNAME3D( SPC )
      END DO

      IF ( W_VEL ) THEN   ! for W_YAMO
         SPC = FINI + 1 
         VTYPE3D( SPC ) = M3REAL
         VNAME3D( SPC ) = 'W_VEL'
         UNITS3D( SPC ) = 'm/s'
         VDESC3D( SPC ) = 'Derived vertical velocity component'
      END IF

C create header

      IF ( IO_PE_INCLUSIVE ) THEN   ! open new

!         IF ( .NOT. OPEN3( CTM_CONC_1, FSNEW3, PNAME ) ) THEN
!            XMSG = 'Could not open ' // CTM_CONC_1
!            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
         
         IF ( .NOT. OPEN3( MIX_FILE_1, FSNEW3, PNAME ) ) THEN
            XMSG = 'Could not open ' // MIX_FILE_1 // ' file'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

!         IF ( .NOT. FLUSH3 ( CTM_CONC_1 ) ) THEN
!            XMSG = 'Could not sync to disk ' // CTM_CONC_1
!            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

      END IF

C write the initial concentrations as step 0 on the conc file
C (inital data assumed to be in correct output units)

      ALLOCATE ( DBUFF( MY_NCOLS,MY_NROWS,NLAYS ), STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'Failure allocating DBUFF'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

      DO SPC = 1, N_GC_CONC
         V = GC_STRT - 1 + GC_CONC_MAP( SPC )

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,V )

         IF ( .NOT. WRITE3( MIX_FILE_1, GC_CONC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' //
     &           TRIM( GC_CONC( SPC ) ) //
     &           ' to ' // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      DO SPC = 1, N_AE_CONC
         V = AE_STRT - 1 + AE_CONC_MAP( SPC )

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,V )

         IF ( .NOT. WRITE3( MIX_FILE_1, AE_CONC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' //
     &           TRIM( AE_CONC( SPC ) ) //
     &           ' to ' // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      V = NR_STRT - 1
      DO SPC = 1, N_NR_SPC
         V = V + 1

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,V )

         IF ( .NOT. WRITE3( MIX_FILE_1, NR_SPC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' //
     &           TRIM( NR_SPC( SPC ) ) //
     &           ' to ' // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      V = TR_STRT - 1
      DO SPC = 1, N_TR_SPC
         V = V + 1

         DBUFF = MGRID( 1:MY_NCOLS,1:MY_NROWS,1:NLAYS,V )

         IF ( .NOT. WRITE3( MIX_FILE_1, TR_SPC( SPC ),
     &        JDATE, JTIME, DBUFF ) ) THEN
            XMSG = 'Could not write ' //
     &           TRIM( TR_SPC( SPC ) ) //
     &           ' to ' // MIX_FILE_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END DO

      DEALLOCATE ( DBUFF )

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &     'Timestep written to', TRIM( MIX_FILE_1 ),
     &     'for date and time', JDATE, JTIME
      WRITE( LOGDEV, '(  5X,  A, 1X, I8, ":", I6.6 )' )
     &     'from timestep on initial data files for date and time',
     &     JDATE, JTIME


      END SUBROUTINE OPCONC
