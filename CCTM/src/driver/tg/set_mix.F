C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Set the composition of a packet created in an interior cell by a fill method
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SET_MIX(JDATE, JTIME, P, C, R, L, MIX_METHOD)

      USE PACKETS               ! data and functions for TG packets
      USE CGRID_SPCS            ! CGRID species number and offsets
	  USE GRID_CONF
	  USE UTILIO_DEFN

      IMPLICIT NONE

C Includes:
!      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)
	  INCLUDE SUBST_CONST
	  INCLUDE SUBST_FILES_ID

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      TYPE(Packet) :: P
      INTEGER, INTENT(IN) :: C, R, L ! column, row and layer indices
      CHARACTER(16), INTENT(IN) :: MIX_METHOD ! see PACKETS.F

C Local Variables:
      CHARACTER(16) :: PNAME = 'SET_MIX'
      CHARACTER(96) :: XMSG = ' '

      TYPE(PacketList) :: NEARBY_PACKETS
      TYPE(PacketListMember), POINTER :: ITR
      TYPE(Packet), POINTER :: CLOSEST

      REAL :: MIX(1:NSPCSD)     ! mixing ratios for packet

      INTEGER :: NCX            ! max number of cells from given cell to consider

C-----------------------------------

      SELECT CASE ( MIX_METHOD )

      CASE ( CLOSEST_MIX )

         CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L)
         P%Mix = CLOSEST%Mix

      CASE ( PU_SUN_MIX )

         NCX = 1
         CLOSEST => CLOSEST_PACKET(JDATE, JTIME, C, R, L, NCX, NEARBY_PACKETS)

         CALL DELETE_PACKET_LIST(JDATE, JTIME, NEARBY_PACKETS)

      CASE DEFAULT
         XMSG =  '"' // TRIM(MIX_METHOD) // '"' //
     &        ' is not a valid selection'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)

      END SELECT

      END SUBROUTINE SET_MIX

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Given a target packet and a list of other packet in the same layer, set the
C target packet's composition equal to the list's distance weighted average
C composition
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DIST_WGT_AVG_MIX(P_TARG, PACKET_LIST)

      USE PACKETS               ! data and functions for TG packets

      IMPLICIT NONE

C Arguments:
      TYPE(Packet) :: P_TARG
      TYPE(PacketList) :: PACKET_LIST

C Local Variables:
      TYPE(PacketListMember), POINTER :: ITR
      TYPE(Packet), POINTER :: P

      REAL(8) :: WGT            ! weight (inverse distance squared)
      REAL(8) :: SUM_WGTS       ! sum of weights

C-----------------------------------

      P_TARG%Mix = 0.0
      SUM_WGTS = 0.0

      ITR => PACKET_LIST%Head
      DO
         IF ( .NOT. ASSOCIATED(ITR) ) EXIT
         P => ITR%P

         WGT = 1.0 / ((P%X - P_TARG%X)**2 + (P%Y - P_TARG%Y)**2)

         P_TARG%Mix = P_TARG%Mix + (WGT * P%Mix)
         SUM_WGTS = SUM_WGTS + WGT

         ITR => ITR%Next
      END DO

      P_TARG%Mix = P_TARG%Mix / SUM_WGTS

      END SUBROUTINE DIST_WGT_AVG_MIX
