!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! written by Carlo Cavazzoni

MODULE kcp_interfaces
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: kcp_move_electrons
  PUBLIC :: kcp_printout_new
  PUBLIC :: kcp_runcp_uspp
  PUBLIC :: kcp_eigs
  !
  INTERFACE kcp_move_electrons
    SUBROUTINE kcp_move_electrons_x( &
         nfi, tfirst, tlast, b1, b2, b3, fion, c0_bgrp, cm_bgrp, phi_bgrp, enthal, enb, &
            &  enbi, fccc, ccc, dt2bye, stress,l_cprestart )
         USE kinds,         ONLY: DP
         IMPLICIT NONE
         INTEGER,  INTENT(IN)    :: nfi
         LOGICAL,  INTENT(IN)    :: tfirst, tlast
         REAL(DP), INTENT(IN)    :: b1(3), b2(3), b3(3)
         REAL(DP)                :: fion(:,:)
         COMPLEX(DP)             :: c0_bgrp(:,:), cm_bgrp(:,:), phi_bgrp(:,:)
         REAL(DP), INTENT(IN)    :: dt2bye
         REAL(DP)                :: fccc, ccc
         REAL(DP)                :: enb, enbi
         REAL(DP)                :: enthal
         REAL(DP)                :: stress(3,3)
         LOGICAL, INTENT(in)     :: l_cprestart
     END SUBROUTINE
  END INTERFACE
  !
  INTERFACE kcp_printout_new
     SUBROUTINE kcp_printout_new_x &
         ( nfi, tfirst, tfilei, tprint, tps, h, stress, tau0, vels, &
           fion, ekinc, temphc, tempp, temps, etot, enthal, econs, econt, &
           vnhh, xnhh0, vnhp, xnhp0, vnhe, xnhe0, atot, ekin, epot, print_forces, print_stress,tstdout )
         USE kinds,          ONLY: DP
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: nfi
         LOGICAL, INTENT(IN)  :: tfirst, tfilei, tprint
         REAL(DP), INTENT(IN) :: tps
         REAL(DP), INTENT(IN) :: h( 3, 3 )
         REAL(DP), INTENT(IN) :: stress( 3, 3 )
         REAL(DP), INTENT(IN) :: tau0( :, : )  ! real positions
         REAL(DP), INTENT(IN) :: vels( :, : )  ! scaled velocities
         REAL(DP), INTENT(IN) :: fion( :, : )  ! real forces
         REAL(DP), INTENT(IN) :: ekinc, temphc, tempp, etot, enthal, econs, econt
         REAL(DP), INTENT(IN) :: temps( : ) ! partial temperature for different ionic species
         REAL(DP), INTENT(IN) :: vnhh( 3, 3 ), xnhh0( 3, 3 ), vnhp( 1 ), xnhp0( 1 ), vnhe, xnhe0
         REAL(DP), INTENT(IN) :: atot! enthalpy of system for c.g. case
         REAL(DP), INTENT(IN) :: ekin
         REAL(DP), INTENT(IN) :: epot ! ( epseu + eht + exc )
         LOGICAL, INTENT(IN) :: print_forces, print_stress, tstdout
     END SUBROUTINE
  END INTERFACE
  !
  INTERFACE kcp_runcp_uspp
     SUBROUTINE kcp_runcp_uspp_x &
         ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec_bgrp, c0_bgrp, cm_bgrp, fromscra, restart )
         USE kinds,             ONLY: DP
         IMPLICIT NONE
         integer, intent(in) :: nfi
         real(DP) :: fccc, ccc
         real(DP) :: ema0bg(:), dt2bye
         real(DP) :: rhos(:,:)
         real(DP) :: bec_bgrp(:,:)
         complex(DP) :: c0_bgrp(:,:), cm_bgrp(:,:)
         logical, optional, intent(in) :: fromscra
         logical, optional, intent(in) :: restart
     END SUBROUTINE
  END INTERFACE
  !
  INTERFACE kcp_eigs
     SUBROUTINE kcp_eigs_x( nfi, lambdap, lambda, desc )
         USE kinds,            ONLY: DP
         USE descriptors,      ONLY: la_descriptor
         IMPLICIT NONE
         INTEGER :: nfi
         REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
         TYPE(la_descriptor), INTENT(IN) :: desc( : )
     END SUBROUTINE
  END INTERFACE

END MODULE
