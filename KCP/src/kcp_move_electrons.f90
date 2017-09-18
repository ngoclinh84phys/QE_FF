!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE kcp_move_electrons_x ( nfi, tfirst, tlast, b1, b2, b3, fion, c0_bgrp, &
            cm_bgrp, phi_bgrp, enthal, enb, enbi, fccc, ccc, dt2bye, stress, l_cprestart )
  !----------------------------------------------------------------------------
  !
  ! ... this routine updates the electronic degrees of freedom
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE control_flags,        ONLY : lwf, tfor, tprnfor, thdyn
  USE cg_module,            ONLY : tcg
  USE cp_main_variables,    ONLY : eigr, irb, eigrb, rhog, rhos, rhor, drhor, &
                                   drhog, sfac, ema0bg, bec_bgrp, becdr_bgrp, &
                                   taub, lambda, lambdam, lambdap, vpot, dbec, descla
  USE cell_base,            ONLY : omega, ibrav, h, press
  USE uspp,                 ONLY : becsum, vkb, nkb, nlcc_any
  USE energies,             ONLY : ekin, enl, entropy, etot
  USE electrons_base,       ONLY : nbsp, nspin, f, nudx, nupdwn, nbspx_bgrp, nbsp_bgrp
  USE electrons_base,       ONLY : ispin_bgrp, f_bgrp, nspin, nupdwn_bgrp, iupdwn_bgrp
  USE core,                 ONLY : rhoc
  USE ions_positions,       ONLY : tau0
  USE ions_base,            ONLY : nat
  USE dener,                ONLY : detot, denl, dekin6
  USE efield_module,        ONLY : tefield, ipolp, qmat, gqq, evalue, &
                                   tefield2, ipolp2, qmat2, gqq2, evalue2
  !
  USE wannier_subroutines,  ONLY : get_wannier_center, wf_options, &
                                   write_charge_and_exit, ef_tune
  USE ensemble_dft,         ONLY : compute_entropy2
  USE efield_module,        ONLY : berry_energy, berry_energy2
  USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                   interpolate_lambda
  USE gvecw,                ONLY : ngw
  USE orthogonalize_base,   ONLY : calphi_bgrp
  USE control_flags,        ONLY : force_pairing
  USE cp_interfaces,        ONLY : rhoofr, compute_stress, vofrho, nlfl_bgrp, prefor, nlfq_bgrp
  USE electrons_module,     ONLY : distribute_c, collect_c, distribute_b
  USE gvect,                ONLY : eigts1, eigts2, eigts3 
  USE control_flags,        ONLY : lwfpbe0nscf  ! exx_wf related
  USE wavefunctions_module, ONLY : cv0 ! Lingzhu Kong
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
  !
      use nksic,                    only : do_orbdep, do_innerloop, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, innerloop_init_n, innerloop_cg_ratio, &
                                           vsicpsi, vsic, wtot, fsic, fion_sic, deeq_sic, f_cutoff, &
                                           pink, do_wxd, sizwtot, do_bare_eigs, innerloop_until, &
                                           valpsi, odd_alpha, eodd
  use kcp_electrons_module,     only : wfc_spreads, wfc_centers, icompute_spread
  USE kcp_interfaces,       ONLY : kcp_runcp_uspp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: nfi
  LOGICAL,  INTENT(IN)    :: tfirst, tlast
  REAL(DP), INTENT(IN)    :: b1(3), b2(3), b3(3)
  REAL(DP)                :: fion(:, :)
  COMPLEX(DP)             :: c0_bgrp(:,:), cm_bgrp(:,:), phi_bgrp(:,:)
  REAL(DP), INTENT(IN)    :: dt2bye
  REAL(DP)                :: fccc, ccc
  REAL(DP)                :: enb, enbi
  REAL(DP)                :: enthal
  REAL(DP)                :: stress(3,3)
  LOGICAL, INTENT(in)     :: l_cprestart
  !
  INTEGER      :: i, j, is, n2
  INTEGER      :: ninner
  INTEGER,SAVE :: nouter = 0
  REAL(DP)     :: ei_unp
  !
  CALL start_clock('kcp_move_electrons')
  !
  electron_dynamic: IF ( tcg ) THEN
     !
     CALL kcp_runcg_uspp( nfi, tfirst, tlast, eigr, bec_bgrp, irb, eigrb, &
                      rhor, rhog, rhos, rhoc, eigts1, eigts2, eigts3, sfac, &
                      fion, ema0bg, becdr_bgrp, lambdap, lambda, SIZE(lambda,1), vpot, c0_bgrp, &
                      cm_bgrp, phi_bgrp, dbec, l_cprestart  )
     !
     CALL compute_stress( stress, detot, h, omega )
     !
  ELSE
     !
     IF ( lwf ) &
          CALL get_wannier_center( tfirst, cm_bgrp, bec_bgrp, eigr, &
                                   eigrb, taub, irb, ibrav, b1, b2, b3 )
     !
     CALL rhoofr( nfi, c0_bgrp, irb, eigrb, bec_bgrp, dbec, becsum, rhor, &
                  drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6 )
     !
!=================================================================
!exx_wf related
     IF ( dft_is_hybrid().AND.exx_is_active() ) THEN
        !
        IF ( lwfpbe0nscf ) THEN
           !
           CALL start_clock('exact_exchange')
           CALL exx_es(nfi, c0_bgrp, cv0)
           CALL stop_clock('exact_exchange')
           !
        ELSE
           !
           CALL start_clock('exact_exchange')
           CALL exx_gs(nfi, c0_bgrp)
           CALL stop_clock('exact_exchange')
           !
        ENDIF
        !
     ENDIF
!=================================================================
     ! ... put core charge (if present) in rhoc(r)
     !
     IF ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
     !
     IF ( lwf ) THEN
        !
        CALL write_charge_and_exit( rhog )
        CALL ef_tune( rhog, tau0 )
        !
     ENDIF
     !
     vpot = rhor
     !
     CALL vofrho( nfi, vpot, drhor, rhog, drhog, rhos, rhoc, tfirst, tlast,&
                     eigts1, eigts2, eigts3, irb, eigrb, sfac, &
                     tau0, fion )
     !
     ! compute auxiliary potentials
     !
     IF ( do_orbdep ) THEN
        !
     !   if (odd_nkscalfact) then
           !
     !      odd_alpha(:) = 0.0_DP
           !
     !      call odd_alpha_routine(c0, nbsp, nbspx, lgam, .false.)
           ! 
     !   endif
        !
        !IF ( tens .or. tsmear) THEN
           !
        !   fsic = fmat0_diag
           !
        !ELSE
           !
           fsic = f_bgrp
           !
       ! ENDIF
        !
        IF ( tlast ) THEN
           !
           icompute_spread=.true.
           !
        ENDIF
        !
        CALL nksic_potential( nbsp_bgrp, nbspx_bgrp, c0_bgrp, fsic, bec_bgrp, becsum, deeq_sic, &
                    ispin_bgrp, iupdwn_bgrp, nupdwn_bgrp, rhor, rhoc, wtot, sizwtot, vsic, do_wxd, pink, nudx, &
                    wfc_centers, wfc_spreads, icompute_spread, .false.)
        !
        eodd = SUM (pink(1:nbsp))
        !
        nouter = nouter + 1
        !
        ninner = 0
        !
        IF ( do_innerloop .and. ( nouter == 1 ) ) THEN
           !
           CALL nksic_rot_emin_cg_general(nouter, innerloop_init_n, ninner, etot, innerloop_cg_ratio, &
                                      nbsp_bgrp, nbspx_bgrp, nudx, iupdwn_bgrp, nupdwn_bgrp, ispin_bgrp, c0_bgrp, becsum, bec_bgrp, rhor, rhoc, &
                         vsic, pink, deeq_sic, wtot, fsic, sizwtot, do_wxd, wfc_centers, wfc_spreads, .false.)

           !CALL nksic_rot_emin_cg(nouter, innerloop_init_n, ninner, etot, Omattot, &
           !                       esic_conv_thr, lgam)
           !  
           eodd = SUM (pink(:))
           !
        ENDIF
        !
        WRITE(stdout,'(2I10,2F24.13)') ninner, nouter,etot,sum(pink(:))
        !
        etot = etot + eodd
        !
     ENDIF !if( do_orbdep )
     !
     IF ( lwf ) CALL wf_options( tfirst, nfi, cm_bgrp, becsum, bec_bgrp, dbec, &
                                 eigr, eigrb, taub, irb, ibrav, b1,   &
                                 b2, b3, vpot, drhor, rhog, drhog, rhos, enl, ekin  )
     !
     CALL compute_stress( stress, detot, h, omega )
     !
     enthal = etot + press * omega
     !
     IF ( tefield )  THEN
        !
        CALL berry_energy( enb, enbi, bec_bgrp, c0_bgrp, fion )
        !
        etot = etot + enb + enbi
        !
     ENDIF
     !
     IF ( tefield2 )  THEN
        !
        CALL berry_energy2( enb, enbi, bec_bgrp, c0_bgrp, fion )
        !
        etot = etot + enb + enbi
        !
     ENDIF

     !
     !=======================================================================
     !
     !              verlet algorithm
     !
     !     loop which updates electronic degrees of freedom
     !     cm=c(t+dt) is obtained from cm=c(t-dt) and c0=c(t)
     !     the electron mass rises with g**2
     !
     !=======================================================================
     !
     CALL newd( vpot, irb, eigrb, becsum, fion )
     !
     CALL prefor( eigr, vkb )
     !
     IF( force_pairing ) THEN
        !
        CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, &
                                       rhos, bec_bgrp, c0_bgrp, cm_bgrp, ei_unp )
        !
     ELSE
        !
        CALL kcp_runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec_bgrp, c0_bgrp, cm_bgrp )
        !
     ENDIF
     !
     !----------------------------------------------------------------------
     !                 contribution to fion due to lambda
     !----------------------------------------------------------------------
     !
     ! ... nlfq needs deeq bec
     !
     IF ( tfor .OR. tprnfor ) THEN
        !
        CALL nlfq_bgrp( c0_bgrp, eigr, bec_bgrp, becdr_bgrp, fion )
        !
     ENDIF
     !
     IF ( (tfor.or.tprnfor) .AND. tefield ) &
        CALL bforceion( fion, .TRUE. , ipolp, qmat, bec_bgrp, becdr_bgrp, gqq, evalue )
     !
     IF ( (tfor.or.tprnfor) .AND. tefield2 ) &
        CALL bforceion( fion, .TRUE. , ipolp2, qmat2, bec_bgrp, becdr_bgrp, gqq2, evalue2 )
     !
     IF ( force_pairing ) THEN
        !
        lambda( :, :, 2 ) =  lambda(:, :, 1 )
        lambdam( :, :, 2 ) = lambdam(:, :, 1 )
        !
     ENDIF
     ! 
     IF ( tfor .OR. thdyn ) then
        !
        CALL interpolate_lambda( lambdap, lambda, lambdam )
        !
     ELSE
        ! 
        ! take care of the otherwise uninitialized lambdam
        !
        lambdam = lambda
        !
     ENDIF
     !
     ! ... calphi calculates phi
     ! ... the electron mass rises with g**2
     !
     CALL calphi_bgrp( c0_bgrp, ngw, bec_bgrp, nkb, vkb, phi_bgrp, nbspx_bgrp, ema0bg )
     !
     ! ... begin try and error loop (only one step!)
     !
     ! ... nlfl and nlfh need: lambda (guessed) becdr
     !
     IF ( tfor .OR. tprnfor ) THEN
        !
        CALL nlfl_bgrp( bec_bgrp, becdr_bgrp, lambda, descla, fion )
        !
     ENDIF
     !
  ENDIF electron_dynamic
  !
  CALL stop_clock('kcp_move_electrons')
  !
  RETURN
  !
END SUBROUTINE kcp_move_electrons_x
