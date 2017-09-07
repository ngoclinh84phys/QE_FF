!
!
!
module kc_mod
      !
      use kinds,                      only : dp
      !
      SAVE
      !
      logical :: ln_kc = .FALSE.
      logical :: ln_ki = .FALSE.
      logical :: ln_kipz = .TRUE.
      logical :: ln_kc_proj = .FALSE.
      logical :: ln_kc_diag = .FALSE.
      logical :: use_wannier_orbs = .FALSE.
      character(len=256) :: localized_orb_file 
      !
      real(dp),    allocatable :: occ_orbs(:) 
      complex(dp), allocatable :: vsicpsi0(:,:) 
      complex(dp), allocatable :: psi0(:,:) 
      !
endmodule kc_mod
!
subroutine allocate_kc_correction ()
      !
      use wvfct,                only : nbnd, npwx
      use kc_mod,               only : vsicpsi0
      use kc_mod,               only : psi0
      use kc_mod,               only : occ_orbs
      !
      allocate (psi0(npwx, nbnd), vsicpsi0(npwx, nbnd))
      allocate (occ_orbs(nbnd))
      !   
endsubroutine allocate_kc_correction
!
subroutine deallocate_kc_correction ()
      !
      use kc_mod,               only : vsicpsi0
      use kc_mod,               only : psi0
      use kc_mod,               only : occ_orbs
      !
      deallocate (psi0, vsicpsi0, occ_orbs)
      !   
endsubroutine deallocate_kc_correction
!
subroutine build_nksic_projection ( nbnd, npw, current_spin, evc, first_call ) 
      !
      use kinds,                      only : dp
      use io_global,                  only : stdout
      use constants,                  only : rytoev
      use fft_base,                   only : dfftp, dffts
      use fft_interfaces,             only : fwfft, invfft
      use gvect,                      only : ngm, nl, nlm
      use scf,                        only : rho, rho_core, rhog_core
      use scf,                        only : scf_type
      use gvect,                      only : gstart
      use lsda_mod,                   only : nspin
      use control_flags,              only : gamma_only
      use wavefunctions_module,       only : psic
      use wvfct,                      only : wg, et, current_k, npwx
      use realus,                     only : fwfft_orbital_gamma, fwfft_orbital_k
      use mp_bands,                   only : intra_bgrp_comm
      use mp,                         only : mp_sum
      use kc_mod,                     only : ln_ki, ln_kipz
      use kc_mod,                     only : vsicpsi0 
      use kc_mod,                     only : psi0 
      use kc_mod,                     only : use_wannier_orbs 
      use kc_mod,                     only : occ_orbs
      !
      implicit none
      !
      integer, intent(in)     :: nbnd, npw
      integer, intent(in)     :: current_spin
      logical, intent(in)     :: first_call
      complex(dp), intent(in) :: evc(npwx, nbnd)
      ! 
      ! in/out vars
      !
      ! local variables
      !
      integer  :: ibnd, jbnd, band_index 
      real(dp) :: focc, focc_i, focc_j, shart, etxc_tot, vtxc_tot
      real(dp) :: delta_et, proj_1, proj_2
      real(dp),    allocatable :: vsic(:,:), pink(:), list_alpha(:)
      real(dp),    allocatable :: vxc_tot(:,:)
      complex(dp), allocatable :: evc_r1(:), evc_r2(:)
      type (scf_type)          :: rhoref
      type (scf_type)          :: rhobar
      type (scf_type)          :: orb_rho
      logical                  :: is_empty
      logical, parameter       :: odd_nkscalfact = .True. 
      real(dp), external       :: get_clock
      !
      ! main body
      !
      WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
      !
      CALL start_clock( 'nksic_drv' )
      !
      ! Fist main part
      !
      IF (first_call) THEN
         !
         allocate(vsic(dfftp%nnr,2), pink(nbnd), list_alpha(nbnd))
         allocate(orb_rho%of_r(dfftp%nnr,2))
         allocate(orb_rho%of_g(ngm,2))
         allocate(rhoref%of_r(dfftp%nnr,nspin))
         allocate(rhoref%of_g(ngm,nspin))
         allocate(rhobar%of_r(dfftp%nnr,nspin))
         allocate(rhobar%of_g(ngm,nspin))
         allocate(evc_r1(dfftp%nnr), evc_r2(dfftp%nnr))
         allocate(vxc_tot(dfftp%nnr,nspin))
         !
         ! reset kc vars in the first call
         !
         occ_orbs(:)   = 0.0_DP
         psi0(:,:)     = (0.0_DP, 0.0_DP)
         vsicpsi0(:,:) = (0.0_DP, 0.0_DP)
         !
         IF (use_wannier_orbs) THEN
            !
            CALL read_orbs_spin ( nbnd, npwx, psi0(:,:), current_spin, '.occ.emp' )
            !
         ELSE
            !
            CALL read_orbs_spin ( nbnd, npwx, psi0(:,:), current_spin, '.ks.occ.emp' )
            !
         ENDIF
         !
         ! read odd alpha from file wrt spin
         !
         CALL read_odd_alpha(nbnd, current_spin, list_alpha)           
         !
         CALL v_xc( rho, rho_core, rhog_core, etxc_tot, vtxc_tot, vxc_tot )
         !
         ! loop over bands (2 ffts at the same time)
         !
         vsicpsi0(:,:) = (0.0_DP, 0.0_DP) 
         !
         DO ibnd = 1, nbnd, 2
            !
            ! compute orbital densities
            !
            CALL get_orbital_density ( npwx, nbnd, ibnd, psi0(:,ibnd), psi0(:,ibnd+1), evc_r1, evc_r2, orb_rho)
            !
            ! compute orbital potentials
            !
            inner_loop: DO jbnd = 1, 2
               !
               band_index = ibnd + jbnd - 1
               !
               WRITE( stdout, 9002 ) band_index
               !
               ! this condition is important when n is odd
               !
               if ( band_index > nbnd ) exit inner_loop
               !
               focc = wg(band_index, current_k)*DBLE(nspin)/2.0d0
               occ_orbs(band_index) = focc
               !
               if (focc < 0.0001) then
                  !
                  is_empty = .true.
                  focc = 0.0_DP
                  !
               else
                  !
                  is_empty = .false.
                  focc = 1.0_DP
                  !
               endif
               !
               ! define rhoref and rhobar
               !
               rhoref%of_r = 0.0
               rhoref%of_g = (0.0,0.0)
               rhobar%of_r = 0.0
               rhobar%of_g = (0.0,0.0)
               !
               CALL nksic_get_rhoref ( current_spin, focc, &
                                       orb_rho%of_r(:,jbnd), orb_rho%of_g(:,jbnd), rhoref, rhobar )
               !
               IF (ln_ki .OR. ln_kipz) THEN
                  !
                  CALL nksic_correction_nki ( focc, current_spin, orb_rho%of_r(:,jbnd), orb_rho%of_g(:,jbnd), &
                                              rhoref, rhobar, etxc_tot, vxc_tot, &
                                              vsic(:,jbnd), pink(band_index), is_empty )
                  !
               ENDIF
               !
               IF (ln_kipz) THEN
                  ! 
                  CALL nksic_correction_nkipz( focc, current_spin, orb_rho%of_r(:,jbnd), orb_rho%of_g(:,jbnd), &
                                               vsic(:,jbnd), pink(band_index), is_empty )
                  !
               ENDIF
               !
               pink(band_index) = 2.d0*pink(band_index)/nspin
               !
               WRITE(stdout, 9005) band_index, 2.d0*pink(band_index)/nspin*rytoev
               !
            ENDDO inner_loop
            !
            ! vsic|psi0> in r
            !
            IF (gamma_only) THEN
               !
               psic(:) = (0.0_dp, 0.0_dp) 
               psic(:) = CMPLX(vsic(:,1) * DBLE(evc_r1(:)), vsic(:,2) * DBLE(evc_r2(:)))
               !
               ! ... transform psic back in reciprocal space
               !
               CALL fwfft_orbital_gamma(vsicpsi0, ibnd, nbnd)
               !
            ELSE
               !
               psic(:) = (0.0_dp, 0.0_dp) 
               psic(:) = vsic(:,1) * evc_r1(:)
               !
               CALL fwfft_orbital_k (vsicpsi0, ibnd, nbnd)
               ! 
               IF (ibnd < nbnd) THEN
                  !
                  psic(:) = (0.0_dp, 0.0_dp) 
                  psic(:) = vsic(:,2) * evc_r2(:)
                  !
                  CALL fwfft_orbital_k(vsicpsi0, ibnd+1, nbnd)
                  !
               ENDIF
               !
            ENDIF 
            !
         ENDDO
         !
         ! computing orbital dependent alpha
         !
         DO ibnd = 1, nbnd
            !
            vsicpsi0(:,ibnd) = vsicpsi0(:,ibnd)*list_alpha(ibnd)
            !
         ENDDO
         !
         ! save vsicpsi to files for restarting ?
         ! 
         DEALLOCATE(vsic, pink, list_alpha)
         DEALLOCATE(rhoref%of_r, rhoref%of_g)
         DEALLOCATE(rhobar%of_r, rhobar%of_g)
         DEALLOCATE(orb_rho%of_r, orb_rho%of_g)
         DEALLOCATE(vxc_tot)
         DEALLOCATE(evc_r1,evc_r2)
         !
      ENDIF
      !
      ! Second main part
      ! 
      WRITE( stdout, 9003 ) 
      !
      DO ibnd = 1, nbnd ! for evc
         !
         focc_i = wg(ibnd, current_k)*DBLE(nspin)/2.0d0
         ! 
         delta_et = 0.0_DP
         !
         kc_loop:DO jbnd = 1, nbnd  ! for evc0
            !
            focc_j = occ_orbs(jbnd)
            !
            IF (focc_j == focc_i) THEN
               !
               ! <vsicpsi0|evc> 
               !
               IF (gamma_only) THEN
                  ! 
                  proj_1 = 2.0 * DBLE ( DOT_PRODUCT( evc(1:npwx,ibnd), vsicpsi0(1:npwx,jbnd)))
                  !
                  IF ( gstart == 2 ) THEN
                     !
                     proj_1 = proj_1  -  DBLE ( CONJG( evc(1,ibnd) ) * vsicpsi0(1,jbnd) )
                     !
                  ENDIF
                  !  
               ELSE
                  !
                  proj_1 = DBLE ( DOT_PRODUCT( evc(1:npwx,ibnd),  vsicpsi0(1:npwx,jbnd)) )
                  !
               ENDIF
               !
               ! <psi0|evc>
               !
               IF (gamma_only) THEN
                  ! 
                  proj_2 = 2.0 * DBLE ( DOT_PRODUCT( evc(1:npwx,ibnd), psi0(1:npwx,jbnd)))
                  !
                  IF ( gstart == 2 ) THEN
                     !
                     proj_2 = proj_2  -  DBLE ( CONJG( evc(1,ibnd) ) * psi0(1,jbnd) )
                     !
                  ENDIF
                  !  
               ELSE
                  !
                  proj_2 = DBLE ( DOT_PRODUCT( evc(1:npwx,ibnd), psi0(1:npwx,jbnd)))
                  !
               ENDIF 
               !
               CALL mp_sum( proj_1, intra_bgrp_comm )
               CALL mp_sum( proj_2, intra_bgrp_comm )
               !
               delta_et = delta_et + proj_1*proj_2
               !
            ENDIF
            !
         ENDDO kc_loop
         !    
         WRITE(stdout, 9004) ibnd, et(ibnd, current_k)*rytoev, (et(ibnd, current_k) + delta_et)*rytoev
         !
         et(ibnd, current_k) = et(ibnd, current_k) + delta_et
         !
      ENDDO
      ! 
      RETURN
      !
      CALL stop_clock( 'nksic_drv' )
      ! 
9000 FORMAT(/'     Applying Koopmans correction,  total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Correcting for band # ', I5)
9003 FORMAT(/'     Eigenvalue    KS   and  KC  ')
9004 FORMAT(/'     ', I5 , 5x,  F12.8, 5x,  F12.8)
9005 FORMAT(/'    Self-Hatree', 5x, I5 , 5x,  F12.8, ' eV')
      !
      RETURN
      !
end subroutine build_nksic_projection
!
!
!
subroutine nksic_potential( nbnd, npw, current_spin, psi, hpsi) 
      !
      use kinds,                      only : dp
      use io_global,                  only : stdout
      use constants,                  only : rytoev
      use fft_base,                   only : dfftp, dffts
      use fft_interfaces,             only : fwfft, invfft
      use gvect,                      only : ngm, nl, nlm
      use scf,                        only : rho, rho_core, rhog_core
      use scf,                        only : scf_type
      use lsda_mod,                   only : nspin
      use control_flags,              only : gamma_only
      use wavefunctions_module,       only : psic
      use wvfct,                      only : wg, current_k
      use realus,                     only : fwfft_orbital_gamma, fwfft_orbital_k
      use kc_mod,                     only : ln_ki, ln_kipz
      !
      implicit none
      !
      integer, intent(in)     :: nbnd, npw
      integer, intent(in)     :: current_spin
      complex(dp), intent(in) :: psi(npw,nbnd)
      complex(dp), intent(inout) :: hpsi(npw,nbnd)
      ! 
      ! in/out vars
      !
      ! local variables
      !
      integer  :: ibnd, jbnd, band_index 
      real(dp) :: focc, shart, etxc_tot, vtxc_tot
      real(dp),  allocatable   :: vsic(:,:), pink(:), vxc_tot(:,:)
      real(dp),  allocatable   :: list_alpha(:)
      complex(dp), allocatable :: evc_r1(:), evc_r2(:)
      complex(dp), allocatable :: vsicpsi(:,:)
      type (scf_type)          :: rhoref
      type (scf_type)          :: rhobar
      type (scf_type)          :: orb_rho 
      logical                  :: is_empty
      logical, parameter       :: odd_nkscalfact = .True. 
      real(dp), external       :: get_clock
      !
      allocate(vsic(dfftp%nnr,2), pink(nbnd))
      allocate(orb_rho%of_r(dfftp%nnr,2))
      allocate(orb_rho%of_g(ngm,2))
      allocate(rhoref%of_r(dfftp%nnr,2))
      allocate(rhoref%of_g(ngm,2))
      allocate(rhobar%of_r(dfftp%nnr,2))
      allocate(rhobar%of_g(ngm,2))
      allocate(evc_r1(dfftp%nnr), evc_r2(dfftp%nnr))
      allocate(vsicpsi(npw,nbnd))
      allocate(vxc_tot(dfftp%nnr,nspin))
      allocate(list_alpha(nbnd))
      !
      ! main body
      !
      WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
      !
      CALL start_clock( 'nksic_drv' )
      !
      CALL v_xc( rho, rho_core, rhog_core, etxc_tot, vtxc_tot, vxc_tot )
      !
      ! loop over bands (2 ffts at the same time)
      !
      vsicpsi(:,:) = (0.0_DP, 0.0_DP) 
      !
      DO ibnd = 1, nbnd, 2
         !
         ! compute orbital densities
         !
         CALL get_orbital_density ( npw, nbnd, ibnd, psi(:,ibnd), psi(:,ibnd+1), evc_r1, evc_r2, orb_rho)
         !
         ! compute orbital potentials
         !
         inner_loop: DO jbnd = 1, 2
           !
           band_index = ibnd + jbnd - 1
           !
           !WRITE( stdout, 9002 ) band_index
           !
           ! this condition is important when n is odd
           !
           if ( band_index > nbnd ) exit inner_loop
           !
           ! note: iupdwn(2) is set to zero if nspin = 1
           !
           focc = wg(band_index, current_k)*DBLE(nspin)/2.0d0
           !
           if (focc < 0.0001) then
              !
              is_empty = .true.
              focc = 0.0_DP
              !
           else
              !
              is_empty = .false.
              focc = 1.0_DP
              !
           endif
           !
           ! define rhoref and rhobar
           !
           rhoref%of_r = 0.0
           rhoref%of_g = (0.0,0.0)
           rhobar%of_r = 0.0
           rhobar%of_g = (0.0,0.0)
           !
           CALL nksic_get_rhoref ( current_spin, focc, &
                                   orb_rho%of_r(:,jbnd), orb_rho%of_g(:,jbnd), rhoref, rhobar )
           !
           IF (ln_ki .OR. ln_kipz) THEN
              !
              CALL nksic_correction_nki ( focc, current_spin, orb_rho%of_r(:,jbnd), orb_rho%of_g(:,jbnd), &
                                          rhoref, rhobar, etxc_tot, vxc_tot, &
                                          vsic(:,jbnd), pink(band_index), is_empty )
              !
           ENDIF
           !
           IF (ln_kipz) THEN
              ! 
              CALL nksic_correction_nkipz( focc, current_spin, orb_rho%of_r(:,jbnd), orb_rho%of_g(:,jbnd), &
                                           vsic(:,jbnd), pink(band_index), is_empty )
              !
           ENDIF
           !
           pink(band_index) = 2.d0*pink(band_index)/nspin
           !
           WRITE( stdout, 9002 ) band_index, 2.d0*pink(band_index)/nspin*rytoev
           !
         ENDDO inner_loop
         !
         ! vsic|psi> in r
         !
         IF (gamma_only) THEN
            !
            psic(:) = (0.0_dp, 0.0_dp) 
            psic(:) = CMPLX(vsic(:,1) * DBLE(evc_r1(:)), vsic(:,2) * DBLE(evc_r2(:)))
            !
            ! ... transform psic back in reciprocal space and assign it to hpsi
            !
            CALL fwfft_orbital_gamma(vsicpsi, ibnd, nbnd)
            !
         ELSE
            !
            psic(:) = (0.0_dp, 0.0_dp) 
            psic(:) = vsic(:,1) * evc_r1(:)
            !
            CALL fwfft_orbital_k (vsicpsi, ibnd, nbnd)
            ! 
            IF (ibnd < nbnd) THEN
               !
               psic(:) = (0.0_dp, 0.0_dp) 
               psic(:) = vsic(:,2) * evc_r2(:)
               !
               CALL fwfft_orbital_k(vsicpsi, ibnd+1, nbnd)
               !
            ENDIF
            !
         ENDIF 
         !
      ENDDO
      !
      ! read odd alpha from file wrt spin
      !
      CALL read_odd_alpha(nbnd, current_spin, list_alpha)
      ! 
      ! computing orbital dependent alpha
      !
      DO ibnd = 1, nbnd
         !
         vsicpsi(:,ibnd) = vsicpsi(:,ibnd)*list_alpha(ibnd)
         !
      ENDDO
      !
      hpsi(:,:) = hpsi(:,:) + vsicpsi(:,:)
      !
      DEALLOCATE(vsic, pink, list_alpha)
      DEALLOCATE(rhoref%of_r, rhoref%of_g)
      DEALLOCATE(rhobar%of_r, rhobar%of_g)
      DEALLOCATE(orb_rho%of_r, orb_rho%of_g)
      DEALLOCATE(vsicpsi, vxc_tot)
      DEALLOCATE(evc_r1, evc_r2)
      !
      CALL stop_clock( 'nksic_drv' )
      ! 
9000 FORMAT(/'     Applying Koopmans correction,  total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Correcting for band # ', I5, 5x, 'SH', 3x,  F12.8, ' eV')
      !
      RETURN
      !
end subroutine nksic_potential
!
!
!
subroutine get_orbital_density ( npw, nbnd, i1, c1, c2, rc1, rc2, orb_rho)
      !
      use kinds,                      only : dp
      use cell_base,                  only : omega
      use lsda_mod,                   only : nspin
      use fft_base,                   only : dfftp, dffts
      use fft_interfaces,             only : fwfft, invfft
      use gvect,                      only : ngm, nl, nlm
      use gvecs,                      only : nls, nlsm, doublegrid
      use wavefunctions_module,       only : psic
      use control_flags,              only : gamma_only
      use scf,                        only : scf_type
      !
      implicit none
      !
      ! input/output vars
      !
      integer,     intent(in)  :: npw,nbnd,i1
      complex(dp), intent(in)  :: c1(npw),c2(npw)
      complex(dp), intent(out) :: rc1(dfftp%nnr),rc2(dfftp%nnr)
      type (scf_type), intent(inout) :: orb_rho
      !
      ! local vars
      !
      character(20) :: subname='get_orbital_density'
      real(dp)      :: sa1
      complex(dp), allocatable :: psic_aux(:)
      !
      call start_clock( 'nksic_orbrho' )
      !
      sa1 = 1.0d0 / omega
      !
      ! This case should be the one when using NCPP
      !
      IF ( gamma_only ) THEN
         !
         IF (i1 < nbnd) THEN
            !
            psic(:) = ( 0.D0, 0.D0 )
            !
            psic(nls(1:npw))  = c1(1:npw) + ( 0.D0, 1.D0 ) * c2(1:npw)
            psic(nlsm(1:npw)) = CONJG(c1(1:npw) - ( 0.D0, 1.D0 ) * c2(1:npw))
            ! 
         ELSE
            !
            psic(nls (1:npw)) = c1(1:npw) 
            psic(nlsm(1:npw)) = CONJG( c1(1:npw) )
            !    
         ENDIF
         !
         CALL invfft ('Wave', psic, dffts)
         !
         orb_rho%of_r(:,1) = sa1 *  DBLE(psic(:))*DBLE(psic(:)) 
         orb_rho%of_r(:,2) = sa1 * AIMAG(psic(:))*AIMAG(psic(:))
         !
         rc1(:) = CMPLX(DBLE(psic(:)), 0.0_DP )
         rc2(:) = CMPLX(AIMAG(psic(:)), 0.0_DP )
         !
         psic(:) = (0.D0, 0.D0)
         psic(:) = CMPLX( orb_rho%of_r(:,1), 0.D0) 
         CALL fwfft ('Dense', psic, dfftp)
         orb_rho%of_g(1:ngm, 1) = psic(nl(1:ngm))
         !
         psic(:) = (0.D0, 0.D0)
         psic(:) = CMPLX( orb_rho%of_r(:,2), 0.D0)
         CALL fwfft ('Dense', psic, dfftp)
         orb_rho%of_g(1:ngm, 2) = psic(nl(1:ngm))
         !
      ELSE
         !
         ALLOCATE (psic_aux(dfftp%nnr))
         !
         psic(:)     = ( 0.D0, 0.D0 )
         psic_aux(:) = ( 0.D0, 0.D0 )
         !
         psic(nls(1:npw))     = c1(1:npw)
         psic_aux(nls(1:npw)) = c2(1:npw)
         !
         CALL invfft ('Wave', psic, dffts)
         CALL invfft ('Wave', psic_aux, dffts)
         ! 
         orb_rho%of_r(:,1) = sa1*DBLE(psic(:))*DBLE(psic(:)) + &
                            & AIMAG(psic(:))*AIMAG(psic(:)) 
         orb_rho%of_r(:,2) = sa1*DBLE(psic_aux(:))*DBLE(psic_aux(:)) + &
                            & AIMAG(psic_aux(:))*AIMAG(psic_aux(:)) 
         !
         rc1(:) = psic(:)
         rc2(:) = psic_aux(:)
         !
         psic(:) = (0.D0, 0.D0)
         psic(:) = CMPLX( orb_rho%of_r(:,1), 0.D0)
         CALL fwfft ('Dense', psic, dfftp)
         orb_rho%of_g(1:ngm, 1) = psic(nl(1:ngm))
         !
         psic(:) = (0.D0, 0.D0)
         psic(:) = CMPLX( orb_rho%of_r(:,2), 0.D0)
         CALL fwfft ('Dense', psic, dfftp)
         orb_rho%of_g(1:ngm, 2) = psic(nl(1:ngm))
         !
         DEALLOCATE (psic_aux)
         ! 
      ENDIF
      !
      call stop_clock('nksic_orbrho')
      !
      RETURN
      !
end subroutine get_orbital_density
!
!
!
subroutine nksic_get_rhoref(  ispin, f, &
                              orb_rhor, orb_rhog, rhoref, rhobar)
      !
      ! Computes rhoref and rhobar
      !
      use kinds,                      only : dp
      use fft_base,                   only : dfftp,dffts
      use fft_interfaces,             only : fwfft
      use control_flags,              only : gamma_only
      use gvect,                      only : ngm, nl, nlm
      use lsda_mod,                   only : nspin
      use scf,                        only : rho, scf_type
      use wavefunctions_module,       only : psic
      use gvect,                      only : nl
      use cell_base,                  only : omega
      !
      implicit none
      !
      ! input/output vars
      !
      integer,       intent(in)    :: ispin
      real(dp),      intent(in)    :: f
      real(dp),      intent(in)    :: orb_rhor(dfftp%nnr)
      complex(dp),   intent(in)    :: orb_rhog(ngm)
      type (scf_type), intent(inout) :: rhoref
      type (scf_type), intent(inout) :: rhobar
      !
      ! local vars
      !
      integer :: is
      real(dp), parameter :: fref = 1.0_DP 
      !
      ! main body
      !
      call start_clock( 'nksic_get_rhoref' )
      !
      ! rhobar_i = rho - f_i * rho_i
      !
      IF (nspin == 1 ) THEN
         !
         rhobar%of_r(:,1) = rho%of_r(:,1)*0.5D0
         rhobar%of_r(:,2) = rho%of_r(:,1)*0.5D0
         rhobar%of_g(:,1) = rho%of_g(:,1)*0.5D0
         rhobar%of_g(:,2) = rho%of_g(:,1)*0.5D0
         ! 
      ELSE
         !
         rhobar%of_r(:,:) = rho%of_r(:,:)
         rhobar%of_g(:,:) = rho%of_g(:,:)
         !
      ENDIF
      ! 
      rhobar%of_r(:,ispin) = rhobar%of_r(:,ispin) - f * orb_rhor(:)
      rhobar%of_g(:,ispin) = rhobar%of_g(:,ispin) - f * orb_rhog(:)
      !
      ! rhoref = rho + (f_ref - f_i) rho_i = rhobar_i + f_ref * rho_i
      !
      rhoref%of_r(:,:) = rhobar%of_r(:,:)
      rhoref%of_g(:,:) = rhobar%of_g(:,:)
      ! 
      rhoref%of_r(:,ispin) = rhoref%of_r(:,ispin) + fref * orb_rhor(:)
      rhoref%of_g(:,ispin) = rhoref%of_g(:,ispin) + fref * orb_rhog(:)
      !
      CALL stop_clock( 'nksic_get_rhoref' )
      !
      return
      !
endsubroutine nksic_get_rhoref
!
!
!
subroutine nksic_correction_nki( f, ispin, orb_rhor, orb_rhog, &
                                 rhoref, rhobar, etxc_tot, vxc_tot, &
                                 vsic, pink, is_empty )
      !
      ! ... calculate the non-Koopmans (integrated, NKI)
      !     potential from the orbital density
      !
      !     then  rho_ref = rho - rho_i + n_i
      !           rho_bar = rho - rho_i
      !
      use kinds,                      only : dp
      use io_global,                  only : stdout
      use constants,                  only : e2, fpi, hartree_si, electronvolt_si
      use cell_base,                  only : tpiba2, omega
      use lsda_mod,                   only : nspin      
      use fft_base,                   only : dfftp, dffts
      use fft_interfaces,             only : fwfft, invfft
      use gvect,                      only : ngm, nl, nlm, gg, gstart
      use gvecs,                      only : nls, nlsm, doublegrid
      use wavefunctions_module,       only : psic
      use control_flags,              only : gamma_only      
      use mp_bands,                   only : intra_bgrp_comm
      use mp,                         only : mp_sum
      use martyna_tuckerman,          only : wg_corr_h, do_comp_mt
      use scf,                        only : scf_type 
      !
      implicit none
      !
      integer,     intent(in)     :: ispin
      real(dp),    intent(in)     :: f, orb_rhor(dfftp%nnr)
      complex(dp), intent(in)     :: orb_rhog(ngm)
      type (scf_type), intent(in) :: rhoref
      type (scf_type), intent(in) :: rhobar
      real(dp),    intent(in)     :: etxc_tot
      real(dp),    intent(in)     :: vxc_tot(dfftp%nnr,nspin)
      real(dp),    intent(out)    :: vsic(dfftp%nnr)
      real(dp),    intent(out)    :: pink
      logical,     intent(in)     :: is_empty
      !
      integer       :: ig, nspin_save
      real(dp)      :: fact, ehele, etmp, eh_corr, w2cst
      real(dp)      :: etxc_ref, vtxc_ref, etxc_bar, vtxc_bar
      real(dp), allocatable :: vxc_ref(:,:), vxc_bar(:,:), rhor_core_aux(:)
      complex(dp), allocatable :: rhog_core_aux(:)
      complex(dp), allocatable :: vtmp(:), vcorr(:)
      logical, parameter :: do_comp_mt_kc = .true.
      !
      !==================
      ! main body
      !==================
      !
      CALL start_clock( 'nki_corr' )
      CALL start_clock( 'nki_corr_h' )
      !
      fact=omega/DBLE(dffts%nr1*dffts%nr2*dffts%nr3)
      !
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      !
      pink    = 0.0_dp
      vsic(:) = 0.0_dp
      !
      ! compute self-hartree contributions
      !
      do ig = gstart, ngm
         !
         vtmp(ig) = e2*orb_rhog(ig)*fpi/( tpiba2*gg(ig) )
         !  
      enddo
      !
      if ( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      !
      ! compute periodic corrections
      !
      IF (do_comp_mt) THEN
         !
         CALL wg_corr_h (omega, ngm, orb_rhog, vcorr, eh_corr)
         !
         vtmp(:) = vtmp(:) + vcorr(:)
         !
      endif
      !
      psic(:) = (0.0_dp, 0.0_dp)
      !
      do ig = 1, ngm
         !
         psic(nl(ig)) = vtmp(ig)  
         !
      enddo
      ! 
      if (gamma_only) then
         !
         do ig = 1, ngm
            !
            psic(nlm(ig)) = conjg(psic(nl(ig)))
            !
         enddo
         !
      endif
      !
      call invfft('Dense',psic, dfftp)
      !
      ! this is just the self-hartree potential
      !
      vsic(:) = (1.0_dp - f) * DBLE( psic(:) )
      !
      ! self-hartree contrib to pink
      ! and w2cst for vsic
      !
      if (gamma_only) THEN
         ! 
         ehele = 2.0 * DBLE ( DOT_PRODUCT( vtmp(1:ngm), orb_rhog(1:ngm)))
         !
         if ( gstart == 2 ) then
            !
            ehele = ehele  -  DBLE ( CONJG( vtmp(1) ) * orb_rhog(1))
            !
         endif
         !  
      else
         !
         ehele = DBLE ( DOT_PRODUCT( vtmp(1:ngm), orb_rhog(1:ngm)))
         !
      endif
      !
      ! -self-hartree energy to be added to the vsic potential
      ! the scalar Hatree term of both empty and occupied states is 
      ! in the same form: -E_H[n_i]
      !
      w2cst = -0.5_dp * ehele * omega
      !
      CALL mp_sum( w2cst, intra_bgrp_comm )
      !
write(stdout,*) 'w2cst ha', w2cst
      vsic(:) = vsic(:) + w2cst
      !
      IF (is_empty) THEN
         !
         ehele = 0.5_dp * ehele * omega / fact
         !
         CALL mp_sum( ehele, intra_bgrp_comm )
         !
      ENDIF
      !
      deallocate(vtmp)
      deallocate(vcorr)
      !
      CALL stop_clock( 'nki_corr_h' )
      !
      CALL start_clock( 'nki_corr_vxc' )
      !
      !   add self-xc contributions
      !
      allocate(rhor_core_aux(dffts%nnr))
      allocate(rhog_core_aux(ngm))
      !
      rhor_core_aux(:) = 0.0_DP
      rhog_core_aux(:) = 0.0_DP
      !
      if (.not. is_empty) then
         !
         allocate(vxc_ref(dffts%nnr,nspin))
         !
         etxc_ref     = etxc_tot
         vxc_ref(:,:) = vxc_tot(:,:)
         !
         nspin_save = nspin
         nspin = 2
         !
         allocate(vxc_bar(dffts%nnr,nspin))
         etxc_bar=0.0_dp
         vtxc_bar=0.0_dp
         vxc_bar =0.0_dp
         !
         call v_xc( rhobar, rhor_core_aux, rhog_core_aux, etxc_bar, vtxc_bar, vxc_bar )
         !
         nspin = nspin_save
         !
      else
         !
         allocate(vxc_bar(dffts%nnr,nspin))
         etxc_bar     = etxc_tot
         vxc_bar(:,:) = vxc_tot(:,:)
         !
         nspin_save = nspin
         nspin = 2
         !
         allocate(vxc_ref(dffts%nnr,nspin))
         etxc_ref=0.0_dp
         vtxc_ref=0.0_dp
         vxc_ref =0.0_dp
         !
         call v_xc( rhoref, rhor_core_aux, rhog_core_aux, etxc_ref, vtxc_ref, vxc_ref )
         !
         nspin = nspin_save
         !
      endif 
      !
      ! update potential and pink (including other constant terms)
      !
      IF (.not. is_empty) THEN
         !
         etmp  = sum(vxc_ref(:,ispin) * orb_rhor(:))*fact
         ! 
         CALL mp_sum( etmp , intra_bgrp_comm )
         !
         w2cst = ( etxc_ref - etxc_bar ) - etmp
         !
         vsic(:) = vsic(:) + w2cst
         !
      ELSE
         !
         etmp  = sum(vxc_ref(:, ispin) * orb_rhor(:))*fact
         !
         CALL mp_sum( etmp , intra_bgrp_comm )
         !
         w2cst = ( etxc_ref - etxc_bar ) - etmp
         !
         etmp  = sum(vxc_tot(:,ispin) * orb_rhor(:)) *fact
         !
         CALL mp_sum( etmp , intra_bgrp_comm )
         !
         pink = etxc_ref - etxc_bar - etmp + ehele
         !
         vsic(:) = vsic(:) &
                        + vxc_ref(:,ispin) - vxc_tot(:,ispin) + w2cst
      ENDIF
      !
write(stdout,*) 'w2cst, etxc_ref, etxc_bar, etmp', w2cst, etxc_ref, etxc_bar,  etmp
      call stop_clock( 'nki_corr_vxc' )
      !
      deallocate(vxc_bar)
      deallocate(vxc_ref)
      deallocate(rhor_core_aux)
      deallocate(rhog_core_aux)
      !
      CALL stop_clock( 'nki_corr' )
      !
      return
      !
endsubroutine nksic_correction_nki
!
!
!
subroutine nksic_correction_nkipz( f, ispin, orb_rhor, orb_rhog, vsic, pink, is_empty )
      !
      ! ... calculate the non-Koopmans PZ term in KIPZ 
      !     potential from the orbital density
      !
      use kinds,                      only : dp
      use io_global,                  only : stdout
      use constants,                  only : e2, fpi, hartree_si, electronvolt_si
      use cell_base,                  only : tpiba2, omega
      use lsda_mod,                   only : nspin      
      use fft_base,                   only : dfftp, dffts
      use fft_interfaces,             only : fwfft, invfft
      use gvect,                      only : ngm, nl, nlm, gg, gstart
      use gvecs,                      only : nls, nlsm, doublegrid
      use wavefunctions_module,       only : psic
      use control_flags,              only : gamma_only      
      use mp_bands,                   only : intra_bgrp_comm
      use mp,                         only : mp_sum
      use martyna_tuckerman,          only : wg_corr_h, do_comp_mt
      use scf,                        only : scf_type 
      !
      implicit none
      !
      integer,     intent(in)     :: ispin
      real(dp),    intent(in)     :: f, orb_rhor(dfftp%nnr)
      complex(dp), intent(in)     :: orb_rhog(ngm)
      real(dp),    intent(inout)  :: vsic(dfftp%nnr)
      real(dp),    intent(inout)  :: pink
      logical,     intent(in)     :: is_empty
      !
      integer         :: ig, nspin_save
      real(dp)        :: fact, ehele, etmp, eh_corr, w2cst
      real(dp)        :: etxc_orb, vtxc_orb 
      type (scf_type) :: rho_orb
      real(dp), allocatable    :: vxc_orb(:,:)
      real(dp), allocatable    :: rhor_core_aux(:)
      complex(dp), allocatable :: rhog_core_aux(:)
      complex(dp), allocatable :: vtmp(:), vcorr(:)
      !
      !==================
      ! main body
      !==================
      !
      CALL start_clock( 'nkipz_corr' )
      CALL start_clock( 'nkipz_corr_h' )
      !
      fact=omega/DBLE(dffts%nr1*dffts%nr2*dffts%nr3)
      !
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      !
      ! compute self-hartree contributions
      !
      do ig = gstart, ngm
         !
         vtmp(ig) = e2*orb_rhog(ig)*fpi/( tpiba2*gg(ig) )
         !  
      enddo
      !
      if ( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      !
      ! compute periodic corrections
      !
      IF (do_comp_mt) THEN
         !
         CALL wg_corr_h (omega, ngm, orb_rhog, vcorr, eh_corr)
         !
         vtmp(:) = vtmp(:) + vcorr(:)
         !
      endif
      !
      psic(:) = (0.0_dp, 0.0_dp)
      !
      do ig = 1, ngm
         !
         psic(nl(ig)) = vtmp(ig)  
         !
      enddo
      ! 
      if (gamma_only) then
         !
         do ig = 1, ngm
            !
            psic(nlm(ig)) = conjg(psic(nl(ig)))
            !
         enddo
         !
      endif
      !
      call invfft('Dense',psic, dfftp)
      !
      ! this is just the self-hartree potential
      !
      vsic(:) = vsic(:) - DBLE( psic(:) )
      !
      ! self-hartree contrib to pink
      ! and w2cst for vsic
      !
      if (gamma_only) THEN
         ! 
         ehele = 2.0 * DBLE ( DOT_PRODUCT( vtmp(1:ngm), orb_rhog(1:ngm)))
         !
         if ( gstart == 2 ) then
            !
            ehele = ehele  -  DBLE ( CONJG( vtmp(1) ) * orb_rhog(1))
            !
         endif
         !  
      else
         !
         ehele = DBLE ( DOT_PRODUCT( vtmp(1:ngm), orb_rhog(1:ngm)))
         !
      endif
      !
      ! +self-hartree energy to be added to the vsic potential
      ! the scalar Hatree term of both empty and occupied states is 
      ! in the same form: +E_H[n_i]
      !
      w2cst = 0.5_dp * ehele * omega
      !
      CALL mp_sum( w2cst, intra_bgrp_comm )
      !
      vsic(:) = vsic(:) + w2cst
      !
      pink = w2cst
      !
      deallocate(vtmp)
      deallocate(vcorr)
      !
      CALL stop_clock( 'nkipz_corr_h' )
      !
      CALL start_clock( 'nkipz_corr_vxc' )
      !
      !   add self-xc contributions
      !
      allocate(rho_orb%of_r(dfftp%nnr,1))
      allocate(rho_orb%of_g(ngm,1))
      allocate(vxc_orb(dfftp%nnr,1))
      allocate(rhor_core_aux(dffts%nnr))
      allocate(rhog_core_aux(ngm))
      !
      rhor_core_aux(:) = 0.0_DP
      rhog_core_aux(:) = 0.0_DP
      !
      nspin_save = nspin
      nspin = 1
      !
      rho_orb%of_r(:,1) = orb_rhor(:) 
      rho_orb%of_g(:,1) = orb_rhog(:) 
      !
      etxc_orb=0.0_dp
      vtxc_orb=0.0_dp
      vxc_orb =0.0_dp
      !
      CALL v_xc( rho_orb, rhor_core_aux, rhog_core_aux, etxc_orb, vtxc_orb, vxc_orb )
      !
      nspin = nspin_save
      ! 
      vsic(:) = vsic(:) - vxc_orb (:,1)
      !
      etmp  = sum(vxc_orb(:,1) * orb_rhor(:))*fact
      ! 
      CALL mp_sum( etmp , intra_bgrp_comm )
      !
      w2cst = etxc_orb - etmp
      !
      vsic(:) = vsic(:) - w2cst
      !
      call stop_clock( 'nkipz_corr_vxc' )
      !
      deallocate(rho_orb%of_r,rho_orb%of_g)
      deallocate(vxc_orb)
      deallocate(rhor_core_aux)
      deallocate(rhog_core_aux)
      !
      CALL stop_clock( 'nkipz_corr' )
      !
      return
      !
endsubroutine nksic_correction_nkipz
!
!
SUBROUTINE read_odd_alpha(nbnd, ispin, list_alpha)
  !
  USE kinds,          ONLY : dp
  USE io_global,     ONLY : stdout, ionode, ionode_id 
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE mp_world,      ONLY : world_comm
  !
  IMPLICIT NONE
  !
  integer, intent(in)   :: nbnd, ispin
  real(dp), intent(out) :: list_alpha(nbnd,1)
  !
  integer :: ibnd, num_states
  character(len=4) :: my_spin
  !
  ! read odd alpha from file
  !
  call mp_barrier(world_comm)
  !
  if (ionode) then 
     !
     WRITE(my_spin,'(i1)') ispin
     open(unit = 99,file ='alpha_list_ispin_'//trim(my_spin)//'.dat',form = 'formatted',status = 'old')
     !
     read(99,*), num_states
     ! 
  endif 
  ! 
  call mp_bcast (num_states,ionode_id,world_comm )
  !
  if (num_states .ne. nbnd) then
     ! 
     call errore ('read_odd_alpha', 'number_states in file not equal nbnd in pw.in', num_states) 
     !
  endif
  ! 
  if (ionode) then
     !  
     do ibnd = 1, num_states
        ! 
        read (99, * ) list_alpha(ibnd,1)
        !
     enddo
     !
     close (99)
     !     
  endif
  !
  call mp_bcast(list_alpha,ionode_id,world_comm)
  !
  return
  !
ENDSUBROUTINE read_odd_alpha
!
!
!
SUBROUTINE read_orbs_spin ( ne, npw, c_emp, ispin, extension )
  ! 
  ! ... This subroutine reads wannier orbital from unit emptyunit
  !
  USE kinds,              ONLY: DP
  USE io_global,          ONLY: stdout, ionode, ionode_id
  USE io_files,           ONLY: tmp_dir, prefix
  USE mp,                 ONLY: mp_bcast, mp_sum
  USE lsda_mod,           ONLY: nspin
  USE gvect,              ONLY: ig_l2g
  USE mp_images,          ONLY: intra_image_comm
  USE mp_wave,            ONLY: splitwf
  USE mp,                 ONLY: mp_get, mp_size, mp_rank, mp_sum
  USE klist,              ONLY: ngk, igk_k
  USE wvfct,              ONLY: npwx
  !
  IMPLICIT none
  !
  INTEGER,     INTENT (IN)   :: ne, npw, ispin
  COMPLEX(DP), INTENT(INOUT) :: c_emp(npw,ne)
  CHARACTER(LEN=*),INTENT(IN):: extension 
  !
  LOGICAL :: exst
  INTEGER :: ierr, ig, i, iss, is
  INTEGER :: ngw_rd, ne_rd, ngw_l, ngw_g, emptyunit
  INTEGER :: io_in_parent, nproc_in_parent, me_in_parent
  !
  CHARACTER(LEN=256) :: fileempty, dirname
  COMPLEX(DP), ALLOCATABLE :: ctmp(:)
  INTEGER, ALLOCATABLE :: igk_l2g(:)
  CHARACTER(LEN=4) :: my_spin
  !
  ! ... Subroutine Body
  !
  ALLOCATE ( igk_l2g( npwx ) )
  !
  ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
  !
  igk_l2g = 0
  DO ig = 1, ngk(1)
     igk_l2g(ig) = ig_l2g(igk_k(ig,1))
  ENDDO
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//'/wannier_orbs'//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  !
  WRITE(stdout,*) 'Reading wannier orbs from: ',fileempty
  !
  me_in_parent    = mp_rank( intra_image_comm )
  nproc_in_parent = mp_size( intra_image_comm )
  !    
  io_in_parent = 0
  IF( ionode ) io_in_parent = me_in_parent
  CALL mp_sum( io_in_parent, intra_image_comm )
  !
  emptyunit = 100
  !
  ngw_g    = npw
  ngw_l    = npw
  !
  CALL mp_sum( ngw_g, intra_image_comm )
  !
  ALLOCATE( ctmp(ngw_g) )
  !
  IF ( ionode ) THEN
     !
     INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )
     !
     IF ( EXST ) THEN
        !
        OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
        !
        READ(emptyunit) ngw_rd, ne_rd
        !
        IF ( ne .ne. ne_rd ) THEN
           !   
           EXST = .false.
           ! 
           WRITE( stdout,10)  TRIM(fileempty) 
           WRITE( stdout,20)  ngw_rd, ne_rd
           WRITE( stdout,20)  ngw_g, ne
           !
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
 10  FORMAT('*** WANNIER STATES : wavefunctions dimensions changed  ', A )
 20  FORMAT('*** NGW = ', I8, ' NE = ', I4)
  !
  CALL mp_bcast(exst,   ionode_id, intra_image_comm)
  CALL mp_bcast(ne_rd,  ionode_id, intra_image_comm)
  CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)
  !
  IF (.NOT.exst) THEN 
     CALL errore( 'read_orbs_spin', 'the wannier orbital file is not exist', 1)
  ENDIF
  ! 
  IF ( exst ) THEN
     !     
     DO i = 1, MIN( ne, ne_rd )
        !      
        IF ( ionode ) THEN
           !
           READ(emptyunit) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
           !
        ENDIF
        !
        IF ( i <= ne ) THEN
           !
           CALL splitwf(c_emp(:,i), ctmp, ngw_l, igk_l2g, & 
                       me_in_parent, nproc_in_parent, io_in_parent, intra_image_comm)
           !
        ENDIF
        ! 
     ENDDO
     !
  ENDIF
  !
  IF ( ionode .AND. EXST ) CLOSE(emptyunit)
  !
  DEALLOCATE(ctmp, igk_l2g)
  !
  RETURN
  !
END SUBROUTINE read_orbs_spin
!
!
!
SUBROUTINE read_orbs_from_file ( ne, npw, c_emp )
  ! 
  ! ... This subroutine reads wannier orbital from unit emptyunit
  !
  USE kinds,              ONLY: DP
  USE io_global,          ONLY: stdout, ionode, ionode_id
  USE io_files,           ONLY: tmp_dir, prefix
  USE mp,                 ONLY: mp_bcast, mp_sum
  USE lsda_mod,           ONLY: nspin
  USE gvect,              ONLY: ig_l2g
  USE mp_images,          ONLY: intra_image_comm
  USE mp_wave,            ONLY: splitwf
  USE mp,                 ONLY: mp_get, mp_size, mp_rank, mp_sum
  USE klist,              ONLY: ngk, igk_k
  USE wvfct,              ONLY: npwx
  USE kc_mod,             ONLY: localized_orb_file
  !
  IMPLICIT none
  !
  INTEGER,     INTENT (IN)   :: ne, npw
  COMPLEX(DP), INTENT(INOUT) :: c_emp(npw,ne)
  !
  LOGICAL :: exst
  INTEGER :: ierr, ig, i, iss, is
  INTEGER :: ngw_rd, ne_rd, ngw_l, ngw_g, emptyunit
  INTEGER :: io_in_parent, nproc_in_parent, me_in_parent
  !
  CHARACTER(LEN=256) :: fileempty, dirname
  COMPLEX(DP), ALLOCATABLE :: ctmp(:)
  INTEGER, ALLOCATABLE :: igk_l2g(:)
  CHARACTER(LEN=4) :: my_spin
  !
  ! ... Subroutine Body
  !
  ALLOCATE ( igk_l2g( npwx ) )
  !
  ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
  !
  igk_l2g = 0
  DO ig = 1, ngk(1)
     igk_l2g(ig) = ig_l2g(igk_k(ig,1))
  ENDDO
  !
  fileempty = TRIM(localized_orb_file)
  !
  WRITE(stdout,*) 'Reading localized orbs from: ',fileempty
  !
  me_in_parent    = mp_rank( intra_image_comm )
  nproc_in_parent = mp_size( intra_image_comm )
  !    
  io_in_parent = 0
  IF( ionode ) io_in_parent = me_in_parent
  CALL mp_sum( io_in_parent, intra_image_comm )
  !
  emptyunit = 100
  !
  ngw_g    = npw
  ngw_l    = npw
  !
  CALL mp_sum( ngw_g, intra_image_comm )
  !
  ALLOCATE( ctmp(ngw_g) )
  !
  IF ( ionode ) THEN
     !
     INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )
     !
     IF ( EXST ) THEN
        !
        OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
        !
        READ(emptyunit) ngw_rd, ne_rd
        !
        IF ( ne .ne. ne_rd ) THEN
           !   
           EXST = .false.
           ! 
           WRITE( stdout,10)  TRIM(fileempty) 
           WRITE( stdout,20)  ngw_rd, ne_rd
           WRITE( stdout,20)  ngw_g, ne
           !
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
 10  FORMAT('*** WANNIER STATES : wavefunctions dimensions changed  ', A )
 20  FORMAT('*** NGW = ', I8, ' NE = ', I4)
  !
  CALL mp_bcast(exst,   ionode_id, intra_image_comm)
  CALL mp_bcast(ne_rd,  ionode_id, intra_image_comm)
  CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)
  !
  IF (.NOT.exst) THEN 
     CALL errore( 'read_orbs_spin', 'the wannier orbital file is not exist', 1)
  ENDIF
  ! 
  IF ( exst ) THEN
     !     
     DO i = 1, MIN( ne, ne_rd )
        !      
        IF ( ionode ) THEN
           !
           READ(emptyunit) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
           !
        ENDIF
        !
        IF ( i <= ne ) THEN
           !
           CALL splitwf(c_emp(:,i), ctmp, ngw_l, igk_l2g, & 
                       me_in_parent, nproc_in_parent, io_in_parent, intra_image_comm)
           !
        ENDIF
        ! 
     ENDDO
     !
  ENDIF
  !
  IF ( ionode .AND. EXST ) CLOSE(emptyunit)
  !
  DEALLOCATE(ctmp, igk_l2g)
  !
  RETURN
  !
END SUBROUTINE read_orbs_from_file
