
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
C $Header: /project/yoj/arc/CCTM/src/plrise/smoke/write3_distr.f,v 1.4 2011/10/21 16:11:31 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      FUNCTION WRITE3_DISTR ( FNAME, VNAME, JDATE, JTIME, DIM1, DIM2, DATA )
     &                        RESULT ( SUCCESS )

C-----------------------------------------------------------------------
C Function: Use I/O-API's WRITE3 and avoid using PWRITE3 to write layer
C           fractions to processor unique fake-gridded files

C Note: This function must live in a file that has a ".f" extension!

C Revision History:
C     5 Dec 2007 J.Young: initial implementation
C    16 Feb 2011 S.Roselle: replaced I/O API include files with M3UTILIO
C-----------------------------------------------------------------------

      USE M3UTILIO

      IMPLICIT NONE

      CHARACTER( * ), INTENT( IN ) :: FNAME
      CHARACTER( * ), INTENT( IN ) :: VNAME
      INTEGER,        INTENT( IN ) :: JDATE, JTIME
      INTEGER,        INTENT( IN ) :: DIM1, DIM2
      REAL,           INTENT( IN ) :: DATA( DIM1,DIM2 )
      LOGICAL SUCCESS

      INTEGER, EXTERNAL :: SETUP_LOGDEV

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER LOGDEV
      CHARACTER( 96 ) :: XMSG = ' '

C-----------------------------------------------------------------------

      SUCCESS = .TRUE.

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = SETUP_LOGDEV()
      END If

      IF ( .NOT. WRITE3( FNAME, VNAME, JDATE, JTIME, DATA ) ) THEN
         XMSG = 'Error writing ' // 'to file "' // TRIM( FNAME ) // '."'
         WRITE( LOGDEV,* ) XMSG
         SUCCESS = .FALSE.
      END IF

      RETURN
      END FUNCTION WRITE3_DISTR
