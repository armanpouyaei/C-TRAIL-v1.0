!------------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in           !
!  continuous development by various groups and is based on information        !
!  from these groups: Federal Government employees, contractors working        !
!  within a United States Government contract, and non-Federal sources         !
!  including research institutions.  These groups give the Government          !
!  permission to use, prepare derivative works of, and distribute copies       !
!  of their work in the CMAQ system to the public and to permit others         !
!  to do so.  The United States Environmental Protection Agency                !
!  therefore grants similar permission to use the CMAQ system software,        !
!  but users are requested to provide copies of derivative works or            !
!  products designed to operate in the CMAQ system to the United States        !
!  Government without restrictions as to use by others.  Software              !
!  that is used with the CMAQ system but distributed under the GNU             !
!  General Public License or the GNU Lesser General Public License is          !
!  subject to their copyright restrictions.                                    !
!------------------------------------------------------------------------------!

SUBROUTINE close_files

!-------------------------------------------------------------------------------
! Name:     Close Files
! Purpose:  Close I/O API files.
! Revised:  10 Sep 2001  Original version.  (T. Otte)
!           09 Jan 2002  Changed "stop" statements to calls to "m3exit" for
!                        graceful shut-down of I/O API files.  (T. Otte)
!           29 Aug 2011  Replaced module IODECL3 with I/O API module M3UTILIO.
!                        Replaced F77 character declarations with F90 standard.
!                        Improved error handling.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!-------------------------------------------------------------------------------

  USE m3utilio

  IMPLICIT NONE

  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'CLOSE_FILES'

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   COULD NOT CLOSE I/O API OUTPUT FILES', &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Use I/O API routine to close I/O API files.
!-------------------------------------------------------------------------------

  IF ( .NOT. shut3() ) THEN
    WRITE (*,f9000) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

END SUBROUTINE close_files
