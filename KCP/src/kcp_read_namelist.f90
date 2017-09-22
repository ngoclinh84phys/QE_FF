!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE kcp_read_namelists_module
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE input_parameters
  USE read_namelists_module 
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: kcp_read_namelists
  !
  ! ... modules needed by read_xml.f90
  !
  PUBLIC :: kcp_defaults, kcp_bcast, kcp_defaults 
  !
  !  ... end of module-scope declarations
  !
  !  ----------------------------------------------
  !
  CONTAINS
     !
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist KCP_VARS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE kcp_defaults( prog )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       which_orbdep = 'none'
       !
       nkscalfact   = 0.0
       !
       nkscalfact_odd = .FALSE.
       !
       do_innerloop = .FALSE.
       !
       innerloop_step = 1
       !
       RETURN
       ! 
     ENDSUBROUTINE

     !=----------------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist KCP
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE kcp_bcast()
       !-----------------------------------------------------------------------
       !
       USE io_global, ONLY: ionode_id
       USE mp,        ONLY: mp_bcast
       USE mp_images, ONLY : intra_image_comm
       !
       IMPLICIT NONE
       !
       CALL mp_bcast( which_orbdep,    ionode_id, intra_image_comm )
       CALL mp_bcast( nkscalfact,      ionode_id, intra_image_comm )
       CALL mp_bcast( nkscalfact_odd,  ionode_id, intra_image_comm )
       CALL mp_bcast( do_innerloop,    ionode_id, intra_image_comm )
       CALL mp_bcast( innerloop_step,  ionode_id, intra_image_comm )
       !
       RETURN
       ! 
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Check input values for Namelist KCP
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE kcp_checkin( prog )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)  :: prog   ! ... specify the calling program
       CHARACTER(LEN=20) :: sub_name = 'kcp_checkin'
       !
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       !
       IF ( LEN_TRIM( which_orbdep ) > 0 ) THEN
          !
          DO i = 1, SIZE( which_orbdep_allowed )
             !
             IF( TRIM(which_orbdep) == which_orbdep_allowed(i) ) allowed = .TRUE.
             !
          END DO
          ! 
          IF ( .NOT. allowed ) &
             ! 
             CALL errore( sub_name, ' which_orbdep '''// &
                          & TRIM(which_orbdep)//''' not allowed ',1)
          !
       ENDIF
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     SUBROUTINE kcp_read_namelists( unit )
       !-----------------------------------------------------------------------
       !
       !  this routine reads data from standard input and puts them into
       !  module-scope variables (accessible from other routines by including
       !  this module, or the one that contains them)
       !  ----------------------------------------------
       !
       ! ... declare modules
       !
       USE io_global, ONLY : ionode, ionode_id
       USE mp,        ONLY : mp_bcast
       USE mp_images, ONLY : intra_image_comm
       !
       IMPLICIT NONE
       !
       ! ... declare variables
       !
       INTEGER, INTENT(IN), optional :: unit
       !
       ! ... declare other variables
       !
       INTEGER :: ios
       !
       INTEGER :: unit_loc=5
       !
       IF(PRESENT(unit)) unit_loc = unit
       !
       ! ... default settings for all namelists
       !
       CALL control_defaults( prog )
       CALL system_defaults( prog )
       CALL electrons_defaults( prog )
       CALL ions_defaults( prog )
       CALL cell_defaults( prog )
       !
       ! ... Here start reading standard input file
       !
       ! ... CONTROL namelist
       !
       ios = 0
       IF( ionode ) THEN
          READ( unit_loc, control, iostat = ios )
       END IF
       CALL check_namelist_read(ios, unit_loc, "control")
       !
       CALL control_bcast( )
       CALL control_checkin( prog )
       !
       ! ... fixval changes some default values according to the value
       ! ... of "calculation" read in CONTROL namelist
       !
       CALL fixval( prog )
       !
       ! ... SYSTEM namelist
       !
       ios = 0
       IF( ionode ) THEN
          READ( unit_loc, system, iostat = ios )
       END IF
       CALL check_namelist_read(ios, unit_loc, "system")
       !
       CALL system_bcast( )
       !
       CALL system_checkin( prog )
       !
       ! ... ELECTRONS namelist
       !
       ios = 0
       IF( ionode ) THEN
          READ( unit_loc, electrons, iostat = ios )
       END IF
       CALL check_namelist_read(ios, unit_loc, "electrons")
       !
       CALL electrons_bcast( )
       CALL electrons_checkin( prog )
       !
       ! ... IONS namelist
       !
       ios = 0
       IF ( ionode ) THEN
          !
          IF ( TRIM( calculation ) == 'relax'    .OR. &
               TRIM( calculation ) == 'md'       .OR. &
               TRIM( calculation ) == 'vc-relax' .OR. &
               TRIM( calculation ) == 'vc-md'    .OR. &
               TRIM( calculation ) == 'cp'       .OR. &
               TRIM( calculation ) == 'vc-cp'    .OR. &
               TRIM( calculation ) == 'smd'      .OR. &
               TRIM( calculation ) == 'cp-wf-nscf' .OR. &
               TRIM( calculation ) == 'vc-cp-wf'   .OR. &
               TRIM( calculation ) == 'cp-wf' ) READ( unit_loc, ions, iostat = ios )
  
       END IF
       CALL check_namelist_read(ios, unit_loc, "ions")
       !
       CALL ions_bcast( )
       CALL ions_checkin( prog )
       !
       ! ... CELL namelist
       !
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'vc-relax' .OR. &
              TRIM( calculation ) == 'vc-cp'    .OR. &
              TRIM( calculation ) == 'vc-md'    .OR. &
              TRIM( calculation ) == 'vc-md'    .OR. & 
              TRIM( calculation ) == 'vc-cp-wf') THEN
             READ( unit_loc, cell, iostat = ios )
          END IF
       END IF
       CALL check_namelist_read(ios, unit_loc, "cell")
       !
       CALL cell_bcast()
       CALL cell_checkin( prog )
       !
       ios = 0
       IF( ionode ) THEN
          if (tabps) then
             READ( unit_loc, press_ai, iostat = ios )
          end if
       END IF
       CALL check_namelist_read(ios, unit_loc, "press_ai")
       !
       CALL press_ai_bcast()
       !
       ! ... KOOPMANS CP NAMELIST
       !
       CALL kcp_defaults( prog )
       !
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'cp' ) THEN
             READ( unit_loc, kcp_vars, iostat = ios )
          END IF
       END IF
       CALL check_namelist_read(ios, unit_loc, "kcp_vars")
       !
       CALL kcp_bcast()
       CALL kcp_checkin( prog )
       ! 
       ! ... WANNIER NAMELIST
       !
       CALL wannier_defaults( prog )
       ios = 0
       IF( ionode ) THEN
          IF( TRIM( calculation ) == 'cp-wf'       .OR. &
              TRIM( calculation ) == 'vc-cp-wf'    .OR. &
              TRIM( calculation ) == 'cp-wf-nscf') THEN
             READ( unit_loc, wannier, iostat = ios )
          END IF
       END IF
       CALL check_namelist_read(ios, unit_loc, "wannier")
       !
       CALL wannier_bcast()
       CALL wannier_checkin( prog )
       !
       ! ... WANNIER_NEW NAMELIST
       !
       CALL wannier_ac_defaults( prog )
       ios = 0
       IF( ionode ) THEN
          IF( use_wannier ) THEN
             READ( unit_loc, wannier_ac, iostat = ios )
          END IF
       END IF
       CALL check_namelist_read(ios, unit_loc, "wannier_ac")
       !
       CALL wannier_ac_bcast()
       CALL wannier_ac_checkin( prog )
       !
       RETURN
       !
     END SUBROUTINE kcp_read_namelists
     !
END MODULE kcp_read_namelists_module
