!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE finite_field_mod
  !
  USE kinds,                ONLY : DP
  !
  LOGICAL                  :: l_ff_bse = .FALSE.
  LOGICAL                  :: l_ff_kc  = .FALSE.
  LOGICAL                  :: finite_field_pert = .FALSE.
  INTEGER                  :: orb_index_i
  INTEGER                  :: orb_index_j
  REAL(DP)                 :: scale_factor= 1.0
  REAL(DP)                 :: overlap_thr = 0.0
  REAL(DP)                 :: uPi
  REAL(DP)                 :: rPi
  REAL(DP), ALLOCATABLE    :: rho_init(:,:)
  COMPLEX(DP), ALLOCATABLE :: vcdrho0 (:)
  !
END MODULE finite_field_mod
!
!
!
SUBROUTINE initialize_finite_field_method (status_this_calculation)
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : restart
  USE io_global,            ONLY : stdout 
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : lsda, nspin
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm
  USE finite_field_mod,     ONLY : rho_init
  USE finite_field_mod,     ONLY : orb_index_i
  USE finite_field_mod,     ONLY : orb_index_j
  USE finite_field_mod,     ONLY : vcdrho0
  USE finite_field_mod,     ONLY : l_ff_bse, l_ff_kc
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(OUT) :: status_this_calculation
  !
  ! local variable
  !
  LOGICAL :: do_this_calc
  COMPLEX(DP), ALLOCATABLE :: drhog (:)! (ngm)
  REAL(DP), ALLOCATABLE :: drhor (:)! (dftp%nnr)
  !
  IF (.not. restart) THEN
     !
     CALL errore( 'init_finite_field', 'Finite_field_method need a restart calculation ' ,1 )
     !
  ENDIF
  !
  IF (l_ff_bse) THEN
     !
     ALLOCATE (drhog(ngm))
     !
     CALL orbital_density_bse (orb_index_i, orb_index_j, drhog, do_this_calc)
     !
  ENDIF
  ! 
  IF (l_ff_kc) THEN
     !
     ALLOCATE (drhog(ngm), drhor(dfftp%nnr))
     !
     CALL orbital_density_kc (orb_index_i, orb_index_j, drhog, drhor)
     !
     do_this_calc = .TRUE. 
     !
  ENDIF
  ! 
  IF (.NOT.do_this_calc) THEN
     !
     status_this_calculation = .false.
     WRITE(stdout,*) " "
     WRITE(stdout,*) " ********** IMPORTANT INFORMATION ********* "
     WRITE(stdout,*) " OVERLAP WANNIER ORBITAL : ZERO"
     WRITE(stdout,*) " No finite field pertubation for this orbital couple ! "
     WRITE(stdout,*) " "
     !
     DEALLOCATE (drhog)
     !
     RETURN
     !
  ELSE 
     !
     status_this_calculation = .true.
     WRITE(stdout,*) " "
     WRITE(stdout,*) " ********** IMPORTANT INFORMATION ********* "
     WRITE(stdout,*) " OVERLAP WANNIER ORBITAL : ONE"
     WRITE(stdout,*) " Finite field pertubation is starting ... "
     WRITE(stdout,*) " "
     ! 
  ENDIF   
  !
  ! allocate the global vars
  !
  ALLOCATE (rho_init(dfftp%nnr, nspin))
  ALLOCATE (vcdrho0(dfftp%nnr))
  !
  rho_init(:,:) = 0.0_DP
  rho_init(:,:) = rho%of_r(:,:)
  vcdrho0(:)    = 0.0_DP
  !
  IF (l_ff_bse) THEN
     ! 
     CALL compute_vcdrho_general_bse (drhog, vcdrho0, .true.) 
     !
     DEALLOCATE (drhog)
     !
  ENDIF
  !
  IF (l_ff_kc) THEN
     !
     CALL compute_vcdrho_general_kc (drhog, drhor, vcdrho0, .true.) 
     !
     DEALLOCATE (drhog, drhor)
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE initialize_finite_field_method


SUBROUTINE end_finite_field_method ()
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout 
  USE cell_base,            ONLY : omega, tpiba2
  USE scf,                  ONLY : rho
  USE gvect,                ONLY : ngm, nl
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE finite_field_mod,     ONLY : rho_init
  USE finite_field_mod,     ONLY : vcdrho0
  USE finite_field_mod,     ONLY : orb_index_i
  USE finite_field_mod,     ONLY : orb_index_j
  USE finite_field_mod,     ONLY : scale_factor
  USE finite_field_mod,     ONLY : l_ff_bse, l_ff_kc
  USE finite_field_mod,     ONLY : rPi, uPi
  USE io_bsepot_xml,        ONLY : pdep_merge_and_write_G 
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE lsda_mod,             ONLY : lsda, nspin
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ir
  COMPLEX(DP), ALLOCATABLE :: drhog(:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  COMPLEX(DP), ALLOCATABLE :: vcdrho(:)
  CHARACTER(LEN=6)         :: my_labeli
  CHARACTER(LEN=6)         :: my_labelj
  CHARACTER(LEN=256)       :: fname
  !
  ALLOCATE (aux(dfftp%nnr))
  ! 
  aux(:) = (0.0_DP, 0.0_DP)
  aux(:) = CMPLX((rho%of_r(:,1) - rho_init(:,1)), 0.0_DP) 
  IF (nspin==2) aux(:) = aux(:) + CMPLX ((rho%of_r(:,2) - rho_init(:,2)), 0.0_DP)    
  !
  aux(:) = aux(:)/scale_factor
  !
  IF (l_ff_bse) THEN
     !
     ALLOCATE (drhog(ngm), vcdrho(dfftp%nnr))
     !
     CALL fwfft ('Dense', aux, dfftp)
     !
     drhog(:) = (0.0_DP, 0.0_DP) 
     !
     drhog(:) = aux(nl(:))
     !
     CALL compute_vcdrho_general_bse (drhog, vcdrho, .false.)
     !
     aux(:) = (0.0_DP, 0.0_DP)
     !
     aux(:) = vcdrho(:) + vcdrho0(:)
     !
     WRITE(stdout,*) 'band', orb_index_i, orb_index_j
     WRITE(stdout,*) 'ir, vc*rhoij,     vcXvc*rhoij'
     DO ir=1, 50
        WRITE(stdout,*) ir, dble(vcdrho0(ir)), dble(vcdrho(ir))
     ENDDO
     ! 
     ! bring back to G space
     ! 
     CALL fwfft ('Dense', aux, dfftp)
     !
     ! save to file  
     !
     WRITE(my_labeli,'(i6.6)') orb_index_i
     WRITE(my_labelj,'(i6.6)') orb_index_j
     ! 
     fname = "PWSCF"//TRIM(ADJUSTL(my_labeli))//"_"//TRIM(ADJUSTL(my_labelj))//".dat"
     !
     CALL pdep_merge_and_write_G(fname,aux(nl(:)))
     !
     DEALLOCATE (drhog, vcdrho)
     !
  ENDIF
  !
  IF (l_ff_kc) THEN
     !
     rPi = 0.0_DP
     rPi = SUM(vcdrho0(:)*aux(:))*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     !
     CALL mp_sum( rPi , intra_bgrp_comm ) 
     !
     WRITE(stdout,*) "Relaxed Koopmans rPi=", rPi
     WRITE(stdout,*) ""
     !
     WRITE(stdout,*) " Koopmans screeening coefficience 1+rPi/uPi= ", 1+rPi/uPi
     !
  ENDIF
  !
  CALL deallocate_finite_field_vars ()
  !
  RETURN 
  !
ENDSUBROUTINE end_finite_field_method
!
!
SUBROUTINE deallocate_finite_field_vars ()
  !
  USE finite_field_mod, ONLY : rho_init
  USE finite_field_mod, ONLY : vcdrho0
  !
  ! Save on_off_matrix to file
  !
  IF (ALLOCATED (rho_init)) DEALLOCATE (rho_init) 
  IF (ALLOCATED (vcdrho0 )) DEALLOCATE (vcdrho0) 
  !
ENDSUBROUTINE deallocate_finite_field_vars
!
!
!
SUBROUTINE compute_vcdrho_general_bse (drhog, vcdrho, resert_vloc)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout 
  USE control_flags,        ONLY : gamma_only
  USE constants,            ONLY : fpi, e2
  USE cell_base,            ONLY : omega, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, nl, nlm, gg, gstart
  USE scf,                  ONLY : vltot, v, vrs
  USE lsda_mod,             ONLY : lsda, nspin
  USE martyna_tuckerman,    ONLY : wg_corr_h, do_comp_mt
  USE finite_field_mod,     ONLY : scale_factor
  !
  IMPLICIT NONE
  !
  LOGICAL,     INTENT(IN)    :: resert_vloc 
  COMPLEX(DP), INTENT(IN)    :: drhog(ngm)
  COMPLEX(DP), INTENT(INOUT) :: vcdrho(dfftp%nnr)
  !
  ! local variables
  !
  INTEGER  :: ig 
  REAL(DP) :: eh_corr, fac
  COMPLEX(DP), ALLOCATABLE :: aux_g(:), dvhart(:)
  COMPLEX(DP), ALLOCATABLE :: dvaux_mt(:)
  !
  ALLOCATE (dvhart(dfftp%nnr))
  !
  dvhart(:) = (0.0_DP, 0.0_DP)
  DO ig = gstart, ngm
     !
     fac = e2*fpi/tpiba2
     !
     dvhart(nl(ig)) = fac*drhog(ig)/gg(ig)
     !
  ENDDO
  !
  IF (do_comp_mt) THEN
     !
     ALLOCATE(dvaux_mt(ngm))
     ALLOCATE(aux_g(ngm))
     !
     aux_g(:) = (0.0_DP, 0.0_DP)
     aux_g(:) = drhog(:) 
     ! 
     CALL wg_corr_h (omega, ngm, aux_g, dvaux_mt, eh_corr)
     !
     DO ig = 1, ngm
        !
        dvhart(nl(ig))  = dvhart(nl(ig))  + dvaux_mt(ig)
        !
     ENDDO
     !
     DEALLOCATE(aux_g)
     DEALLOCATE(dvaux_mt)
     !
  ENDIF
  ! 
  IF (gamma_only) THEN
     !
     DO ig = 1, ngm
        !
        dvhart(nlm(ig)) = conjg(dvhart(nl(ig)))
        !
     ENDDO
     !
  ENDIF  
  !
  CALL invfft ('Dense', dvhart, dfftp) 
  !
  vcdrho(:) = DBLE(dvhart(:))
  !
  DEALLOCATE (dvhart)
  !
  IF (resert_vloc) THEN 
     !
     vltot(:) = vltot(:) + DBLE(vcdrho(:))*scale_factor
     !
     CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE compute_vcdrho_general_bse
!
!
!
SUBROUTINE compute_vcdrho_general_kc (drhog, drhor, vcdrho, resert_vloc)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout 
  USE control_flags,        ONLY : gamma_only
  USE constants,            ONLY : fpi, e2
  USE cell_base,            ONLY : omega, tpiba2, alat
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, nl, nlm, gg, g, gstart
  USE scf,                  ONLY : vltot, v, vrs
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : lsda, nspin
  USE wavefunctions_module, ONLY : psic
  USE martyna_tuckerman,    ONLY : wg_corr_h, do_comp_mt
  USE finite_field_mod,     ONLY : scale_factor, uPi, rPi
  USE finite_field_mod,     ONLY : orb_index_i
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE lsda_mod,             ONLY : lsda, nspin
  USE qpoint,               ONLY : xq
  USE noncollin_module,     ONLY : nspin_mag, nspin_gga, nspin_lsda
  USE eqv,                  ONLY : dmuxc
  USE gc_lr,                ONLY : grho, gmag, dvxc_rr, dvxc_sr, &
                                   dvxc_ss, dvxc_s, vsgga, segni
  USE funct,                ONLY : dft_is_gradient
  !
  IMPLICIT NONE
  !
  LOGICAL,     INTENT(IN)    :: resert_vloc 
  COMPLEX(DP), INTENT(IN)    :: drhog(ngm)
  REAL(DP),    INTENT(IN)    :: drhor(dfftp%nnr)
  COMPLEX(DP), INTENT(INOUT) :: vcdrho(dfftp%nnr)
  !
  ! local variables
  !
  INTEGER  :: ir, ig, ispin, is
  REAL(DP) :: eh_corr, fac
  COMPLEX(DP), ALLOCATABLE :: dvaux(:,:),drhoraux(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux_g(:), dvhart(:)
  COMPLEX(DP), ALLOCATABLE :: dvaux_mt(:)
  !
  ! ff method is done at xq(:) = 0
  !
  xq(:) = 0.0_DP
  !
  ! define spin index
  !
  CALL ispin_index(orb_index_i, ispin)
  !
  ALLOCATE (dvhart(dfftp%nnr))
  !
  dvhart(:) = (0.0_DP, 0.0_DP)
  DO ig = gstart, ngm
     !
     fac = e2*fpi/tpiba2
     !
     dvhart(nl(ig)) = fac*drhog(ig)/gg(ig)
     !
  ENDDO
  !
  IF (do_comp_mt) THEN
     !
     ALLOCATE(dvaux_mt(ngm))
     ALLOCATE(aux_g(ngm))
     !
     aux_g(:) = (0.0_DP, 0.0_DP)
     aux_g(:) = drhog(:) 
     ! 
     CALL wg_corr_h (omega, ngm, aux_g, dvaux_mt, eh_corr)
     !
     DO ig = 1, ngm
        !
        dvhart(nl(ig))  = dvhart(nl(ig))  + dvaux_mt(ig)
        !
     ENDDO
     !
     DEALLOCATE(aux_g)
     DEALLOCATE(dvaux_mt)
     !
  ENDIF
  ! 
  IF (gamma_only) THEN
     !
     DO ig = 1, ngm
        !
        dvhart(nlm(ig)) = conjg(dvhart(nl(ig)))
        !
     ENDDO
     !
  ENDIF  
  !
  CALL invfft ('Dense', dvhart, dfftp) 
  !
  vcdrho(:) = DBLE(dvhart(:))
  !
  DEALLOCATE (dvhart)
  !
  ! ... Add xc contribution
  !
  ALLOCATE (dmuxc (dfftp%nnr , nspin_mag , nspin_mag))
  !
  CALL setup_dmuxc()
  !
  CALL setup_dgc()
  !
  fac = 1.d0 / DBLE (nspin_lsda)
  !
  DO is = 1, nspin_mag
     !
    vcdrho(:) = vcdrho(:) + dmuxc(:, ispin, is) * drhor(:)
     !
  ENDDO
  !
  ! Add gradient correction to the response XC potential.
  ! NB: If nlcc=.true. we need to add here its contribution. 
  ! grho contains already the core charge
  !
  IF ( dft_is_gradient() ) THEN 
     ! 
     ALLOCATE (dvaux(dfftp%nnr, nspin_mag))
     ALLOCATE (drhoraux(dfftp%nnr, nspin_mag))
     !
     dvaux(:,:)    = (0.0_DP,0.0_DP)
     drhoraux(:,:) = (0.0_DP,0.0_DP)
     drhoraux(:,ispin) = CMPLX(drhor(:),0.0_DP)
     !
     CALL dgradcorr (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
                     drhoraux, dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, &
                     alat, dvaux)
     !
     vcdrho(:) = vcdrho(:) + DBLE(dvaux(:,ispin))
     !
     DEALLOCATE (dvaux, drhoraux)
     !
  ENDIF
  !
  DEALLOCATE(dmuxc)
  !
  IF (ALLOCATED(segni))   DEALLOCATE (segni)
  IF (ALLOCATED(vsgga))   DEALLOCATE (vsgga)
  IF (ALLOCATED(gmag))    DEALLOCATE (gmag)
  IF (ALLOCATED(dvxc_rr)) DEALLOCATE (dvxc_rr)
  IF (ALLOCATED(dvxc_sr)) DEALLOCATE (dvxc_sr)
  IF (ALLOCATED(dvxc_ss)) DEALLOCATE (dvxc_ss)
  IF (ALLOCATED(dvxc_s))  DEALLOCATE (dvxc_s)
  IF (ALLOCATED(grho))    DEALLOCATE (grho)
  !
  ! ... compute unrelax Pi energy
  !
  uPi=0.D0
  DO ir = 1, dfftp%nnr
     ! 
     uPi = uPi + ( drhor(ir) * DBLE(vcdrho(ir)) )
     !
  ENDDO
  !  
  uPi = uPi *omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  !
  CALL mp_sum(  uPi , intra_bgrp_comm )
  ! 
  WRITE(stdout,*) "unrelaxed Koopmans uPi=", uPi
  !
  IF (resert_vloc) THEN
     !
     vltot(:) = vltot(:) + DBLE(vcdrho(:))*scale_factor
     !
     CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE compute_vcdrho_general_kc
!
! 
!----------------------------------------------------------------------------
SUBROUTINE orbital_density_bse (fixed_band_i, fixed_band_j, drhog, do_this_calc)
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE constants,            ONLY : degspin
  USE fft_base,             ONLY : dfftp, dffts, dtgs
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, nl, nlm
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k, &
                                   nelec, nelup, neldw
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions_module, ONLY : psic
  USE wvfct,                ONLY : npwx, nbnd
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE wavefunctions_module, ONLY : evc
  USE kc_mod,               ONLY : use_wannier_orbs
  USE finite_field_mod,     ONLY : l_ff_bse
  USE mp_images,            ONLY : intra_image_comm
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: fixed_band_i 
  INTEGER, INTENT(IN) :: fixed_band_j 
  COMPLEX(DP), INTENT(OUT) :: drhog (ngm)
  LOGICAL    , INTENT(OUT) :: do_this_calc
  !
  ! ... local variables
  !
  LOGICAL      :: wannierization
  INTEGER      :: ibnd, ik, npw, nbnd_all, max_nocc
  INTEGER      :: ispin, jspin
  INTEGER      :: ir
  INTEGER      :: sorted_fixed_band_i, sorted_fixed_band_j
  INTEGER      :: nbnd_val_i, nbnd_val_j
  REAL(DP)     :: sa1, sb1, overlap
  COMPLEX(DP), ALLOCATABLE :: evc_init(:,:), aux_rho(:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  !
  CALL start_clock( 'orbital_density' )
  !
  sa1 = 1.0_DP/omega
  sb1 = 1.0_DP/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  ik  = 1 ! sofar, single k point only
  npw = ngk(ik)
  !
  max_nocc = MAX(NINT( nelec / degspin ), NINT(nelup), NINT(neldw))
  !
  IF (nbnd .ne. max_nocc) THEN
     CALL errore( 'orbital_density', ' ff_bse need nbnd = nelec/2', nbnd)
  ENDIF
  !
  IF (nspin == 2) THEN
     !
     IF ((fixed_band_i > nelec) .OR. &
         (fixed_band_j > nelec)) &
         CALL errore( 'orbital_density', ' nspin == 2 and fixed_band_i/j > nelec', 1)
     !
     nbnd_all = nbnd*2 
     !
     CALL ispin_index (fixed_band_i, ispin)
     CALL ispin_index (fixed_band_j, jspin)
     !
     ! for orb_i
     !  
     IF (ispin == 1) THEN
        nbnd_val_i = nbnd
        sorted_fixed_band_i =  fixed_band_i
     ENDIF
     !
     IF (ispin == 2) THEN
        nbnd_val_i = nbnd
        sorted_fixed_band_i = fixed_band_i - nbnd
     ENDIF
     !
     ! for orb_j
     !
     IF (jspin == 1) THEN
        nbnd_val_j = nbnd
        sorted_fixed_band_j =  fixed_band_j
     ENDIF
     ! 
     IF (jspin == 2) THEN
         nbnd_val_j = nbnd
         sorted_fixed_band_j = fixed_band_j - nbnd
     ENDIF
     !
  ENDIF 
  !
  IF (nspin == 1) THEN
     !
     ispin = 1
     jspin = 1
     !
     nbnd_all   = nbnd
     !
     nbnd_val_i = nbnd
     nbnd_val_j = nbnd
     !
     IF ((fixed_band_i > nbnd_val_i) .OR. &
         (fixed_band_j > nbnd_val_j)) &
         CALL errore( 'orbital_density', ' nspin == 1 and fixed_band_i/j > nbnd', nbnd)
     !
     sorted_fixed_band_i = fixed_band_i
     sorted_fixed_band_j = fixed_band_j
     !
  ENDIF
  !
write(stdout,*) 'orbital_i', sorted_fixed_band_i, fixed_band_i 
write(stdout,*) 'orbital_j', sorted_fixed_band_j, fixed_band_j
write(stdout,*) 'ij', nbnd_all, nbnd_val_i, nbnd_val_j
  !
  CALL read_overlap_matrix (nbnd_all, sorted_fixed_band_i, sorted_fixed_band_j, &
                            do_this_calc, ispin, jspin, '.wan.occ')
  !
  IF (.NOT. do_this_calc) RETURN
  !
  !
  ALLOCATE (evc_init (npw, nbnd_all)) 
  !
  evc_init(:,:) = (0.0_DP, 0.0_DP)
  !
  IF (use_wannier_orbs) THEN
     !
     IF (nspin == 1) THEN
        !
        CALL read_orbs_spin ( nbnd, npw, evc_init(:,1:nbnd), 1, '.wan.occ' )
        !
     ENDIF
     !
     IF (nspin == 2) THEN
        !
        CALL read_orbs_spin ( nbnd, npw, evc_init(:,1:nbnd), 1, '.wan.occ' )
        CALL read_orbs_spin ( nbnd, npw, evc_init(:,nbnd+1:nbnd_all), 2, '.wan.occ' )
        !
     ENDIF
     !
  ELSE
     !
     IF (nspin == 1) THEN
        !
        evc_init(:,:) = evc(:,:)
        !
     ENDIF
     !
     IF (nspin == 2) THEN
        !
        CALL get_buffer ( evc, nwordwfc, iunwfc, 1 )
        !
        evc_init(:,1:nbnd) = evc(:,1:nbnd)
        ! 
        CALL get_buffer ( evc, nwordwfc, iunwfc, 2 )
        !
        evc_init(:,nbnd+1:nbnd_all) = evc(:,1:nbnd)
        !
     ENDIF
     !
  ENDIF
  !
  ! here is the main part
  !
  ALLOCATE ( aux_rho  (dfftp%nnr))
  ALLOCATE ( aux      (dfftp%nnr))
  !
  IF ( gamma_only ) THEN
     !
     ! orb 1: G -> R
     !
     psic(:) = ( 0.D0, 0.D0 )
     !
     IF ( MOD(sorted_fixed_band_i,2)==0 ) THEN
        !
        psic(nls(1:npw))  = evc_init(1:npw,fixed_band_i-1) + &
                            ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i)
        psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_i-1) - &
                            ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i))
        !
     ELSE
        !
        IF (sorted_fixed_band_i < nbnd_val_i) THEN 
           ! 
           psic(nls(1:npw))  = evc_init(1:npw,fixed_band_i) + &
                               ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i+1)
           psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_i) - &
                               ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i+1))
        ELSE
           !
           psic(nls (1:npw))  = evc_init(1:npw,fixed_band_i)
           psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_i) ) 
           !    
        ENDIF
        !
     ENDIF
     !
     CALL invfft ('Wave', psic, dffts)
     !
     aux(:) = ( 0.D0, 0.D0 )
     !
     aux(:) = psic(:)
     !
     ! orb 2: G -> R
     ! 
     psic(:) = ( 0.D0, 0.D0 )
     !
     IF ( MOD(sorted_fixed_band_j,2)==0 ) THEN
        !
        psic(nls(1:npw))  = evc_init(1:npw,fixed_band_j-1) + &
                            ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_j)
        psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_j-1) - &
                            ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_j))
        !
     ELSE
        !
        IF (sorted_fixed_band_j < nbnd_val_j) THEN
           ! 
           psic(nls(1:npw))  = evc_init(1:npw,fixed_band_j) + &
                               ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_j+1)
           psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_j) - &
                               ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_j+1))
           !
        ELSE
           !
           psic(nls (1:npw))  = evc_init(1:npw,fixed_band_j)
           psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_j) )
           !
        ENDIF
        !
     ENDIF
     !
     CALL invfft ('Wave', psic, dffts)
     !
     ! compute rho(:) = evc_i(:) * evc_j(:)
     !
     aux_rho(:) = (0.0_DP, 0.0_DP)
     !   
     IF ((MOD(sorted_fixed_band_i, 2)==0) .AND. (MOD(sorted_fixed_band_j, 2)==1))  THEN
        ! 
        aux_rho(:)   = sa1 * ( AIMAG(aux(:)) * DBLE(psic(:)) )
        !
     ENDIF 
     !
     IF ((MOD(sorted_fixed_band_i, 2)==1) .AND. (MOD(sorted_fixed_band_j, 2)==1))  THEN 
        !  
        aux_rho(:)   = sa1 * ( DBLE(aux(:))  * DBLE(psic(:)) )
        !
     ENDIF
     !  
     IF ((MOD(sorted_fixed_band_i, 2)==0) .AND. (MOD(sorted_fixed_band_j, 2)==0))  THEN
        ! 
        aux_rho(:)   = sa1 * ( AIMAG(aux(:)) * AIMAG(psic(:)) )
        !
     ENDIF
     !
     IF ((MOD(sorted_fixed_band_i, 2)==1) .AND. (MOD(sorted_fixed_band_j, 2)==0))  THEN
        !  
        aux_rho(:)   = sa1 * ( DBLE(aux(:))  * AIMAG(psic(:)) )
        !
     ENDIF
     !
  ELSE
     !
     ! orb 1: G->R
     ! 
     psic(:) = ( 0.D0, 0.D0 )
     !
     psic(nls(igk_k(1:npw,1))) = evc_init(1:npw,fixed_band_i)
     !
     CALL invfft ('Wave', psic, dffts)
     !
     aux(:) = ( 0.D0, 0.D0 )
     !
     aux(:) = psic(:)
     ! 
     ! orb 2: G->R   
     !
     psic(:) = ( 0.D0, 0.D0 )
     !
     psic(nls(igk_k(1:npw,ik))) = evc_init(1:npw,fixed_band_j)
     !
     CALL invfft ('Wave', psic, dffts)
     !
     ! compute rho(:) = evc_i(:) * evc_j(:)
     ! 
     aux_rho(:) = (0.0_DP, 0.0_DP)
     !
     aux_rho(:) = sa1 * (DBLE(aux(:)) * DBLE(psic(:)) + AIMAG(aux(:)) * AIMAG(psic(:)))
     !
  ENDIF
  !
  ! R->G 
  !
  CALL fwfft ('Dense', aux_rho, dfftp)
  !
  drhog(:) = (0.0_DP, 0.0_DP)
  drhog(:) = aux_rho(nl(:))
  !
  DEALLOCATE ( evc_init ) 
  DEALLOCATE ( aux_rho  )
  DEALLOCATE ( aux )
  !
  CALL stop_clock( 'orbital_density' )
  !
  RETURN
  !
ENDSUBROUTINE
!
!
SUBROUTINE read_io_bsepot_xml (fixed_band_i, fixed_band_j, dvg)
  !
  USE kinds,               ONLY : DP 
  USE io_bsepot_xml,       ONLY : pdep_read_G_and_distribute
  USE gvect,               ONLY : ngm
  !
  INTEGER, INTENT(IN) :: fixed_band_i
  INTEGER, INTENT(IN) :: fixed_band_j
  COMPLEX(DP), INTENT(OUT) :: dvg(ngm)
  !
  CHARACTER(LEN=6)    :: my_labeli
  CHARACTER(LEN=6)    :: my_labelj
  CHARACTER(LEN=256)  :: file_base
  !
  WRITE(my_labeli,'(i6.6)') fixed_band_i
  WRITE(my_labelj,'(i6.6)') fixed_band_j
  !
  file_base = "PWSCF"//TRIM(ADJUSTL(my_labeli))//"_"//TRIM(ADJUSTL(my_labelj))//".dat"
  CALL pdep_read_G_and_distribute(file_base,dvg(:))
  !
  RETURN
  !
END SUBROUTINE
!
!
!
SUBROUTINE read_overlap_matrix ( num_wan, orb_index_i, orb_index_j, do_this_calc, ispin, jspin, extension) 
  !
  !
  ! ...   This subroutine read overlap orb values  
  ! 
  USE kinds,                ONLY: DP
  USE io_global,            ONLY: ionode, ionode_id, stdout
  USE mp_images,            ONLY: intra_image_comm
  USE io_files,             ONLY: prefix, tmp_dir
  USE mp,                   ONLY: mp_bcast
  USE kc_mod,               ONLY: use_wannier_orbs
  USE finite_field_mod,     ONLY: overlap_thr
  !
  IMPLICIT NONE
  ! 
  integer,           intent(in)  :: num_wan, ispin, jspin
  integer,           intent(in)  :: orb_index_i, orb_index_j
  logical,           intent(out) :: do_this_calc
  CHARACTER(LEN=*),  INTENT(IN)  :: extension
  !
  LOGICAL :: exst
  INTEGER :: iw, jw, emptyunit
  REAL(DP) :: overlap_value
  character(len=256) :: fileempty, dirname
  character(len=4)   :: my_spin
  REAL(DP), ALLOCATABLE :: u_matrix(:,:)
  !
  ! ... Subroutine Body
  !
  ! allocate for first time, deallocate in the final
  !
  do_this_calc = .true.
  !
  IF (.NOT.use_wannier_orbs) THEN
     RETURN
  ENDIF
  !
  IF (ispin .ne. jspin) THEN
     CALL errore( 'read_overlap_matrix', 'Sofar only singlet, ispin == jspin, is implemented', 1)
  ENDIF
  !
  ! the main code
  !
  ALLOCATE (u_matrix(num_wan, num_wan))
  !
  u_matrix = 0.0_DP
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//'/overlap_matrix'//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  !
  WRITE(stdout,*) 'Reading wannier orbs from: ',fileempty
  ! 
  emptyunit = 100
  !
  IF ( ionode ) THEN
     !
     INQUIRE( FILE = TRIM(fileempty), EXIST = EXST )
     !
     IF ( EXST ) THEN
        !
        OPEN( UNIT=emptyunit, FILE=TRIM(fileempty), STATUS='OLD', FORM='UNFORMATTED' )
        !
        READ(emptyunit) (( u_matrix(iw,jw), iw=1, num_wan), jw=1, num_wan )
        !
     ENDIF
     !
  ENDIF
  !
  CALL mp_bcast(u_matrix, ionode_id, intra_image_comm)
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  overlap_value = u_matrix(orb_index_i, orb_index_j)   
  !
  IF (overlap_value < overlap_thr) do_this_calc = .false. 
  !
  DEALLOCATE (u_matrix)
  ! 
  RETURN
  !
END SUBROUTINE
!
!
!
SUBROUTINE ispin_index (orb_i, ispin)
  !
  USE lsda_mod,         ONLY : lsda, nspin
  USE wvfct,            ONLY : nbnd
  USE noncollin_module, ONLY : noncolin
  !
  implicit none
  !
  integer,  intent(in)    :: orb_i
  integer,  intent(inout) :: ispin
  !
  IF (nspin == 1) THEN
     !
     ispin = 1
     !
  ENDIF
  !
  IF (nspin == 2) THEN
     !
     ! number states: 2*nbnd 
     ! 1..nbnd for uu 
     ! nbnd+1..2*nbnd for dd 
     !
     IF (orb_i<=nbnd) THEN
        !
        ispin = 1 
        !
     ENDIF
     !
     IF ((nbnd+1 <= orb_i) .AND. (orb_i <= (2*nbnd))) THEN
        !
        ispin = 2 
        !
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
ENDSUBROUTINE ispin_index
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE orbital_density_kc (fixed_band_i, fixed_band_j, drhog, drhor )
  !----------------------------------------------------------------------------
  !
  ! In case of computing screening factor of KCF for each fixed_band_i orbitals,
  ! the fixed_band_j will be used with different meaning, i.e. it defines the 
  ! total number of orbitals which are either wannier orbs or ks orbs, computed
  ! and stored in the save file. 
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE constants,            ONLY : degspin
  USE fft_base,             ONLY : dfftp, dffts, dtgs
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, nl, nlm
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k, &
                                   nelec, nelup, neldw
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions_module, ONLY : psic
  USE wvfct,                ONLY : npwx, nbnd
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE wavefunctions_module, ONLY : evc
  USE kc_mod,               ONLY : use_wannier_orbs
  USE mp_images,            ONLY : intra_image_comm
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: fixed_band_i 
  INTEGER, INTENT(IN) :: fixed_band_j 
  REAL(DP),    INTENT(OUT) :: drhor (dfftp%nnr)
  COMPLEX(DP), INTENT(OUT) :: drhog (ngm)
  !
  ! ... local variables
  !
  LOGICAL      :: wannierization
  INTEGER      :: ibnd, ik, npw, nbnd_all, nbnd_all2, max_nocc
  INTEGER      :: ispin, jspin
  INTEGER      :: ir
  INTEGER      :: sorted_fixed_band_i, sorted_fixed_band_j
  INTEGER      :: nbnd_val_i, nbnd_val_j
  INTEGER      :: num_band_computed
  REAL(DP)     :: sa1, sb1, overlap
  COMPLEX(DP), ALLOCATABLE :: evc_init(:,:), aux_rho(:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  !
  CALL start_clock( 'orbital_density' )
  !
  sa1 = 1.0_DP/omega
  sb1 = 1.0_DP/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  ik  = 1 ! sofar, single k point only
  npw = ngk(ik)
  !
  num_band_computed = fixed_band_j
  !
  max_nocc = MAX(NINT( nelec / degspin ), NINT(nelup), NINT(neldw))
  !
  IF (num_band_computed < max_nocc) THEN
     CALL errore( 'orbital_density_kc', 'check README, fixed_band_j should ge nbnd_occ computed with nscf',1)
  ENDIF  
  ! 
  nbnd_all = num_band_computed*nspin
  !
  ALLOCATE (evc_init(npw, nbnd_all))
  evc_init(:,:) = (0.0_DP, 0.0_DP)
  ! 
  IF (use_wannier_orbs) THEN
     !
     IF (nspin == 1) THEN
        !
        CALL read_orbs_spin ( nbnd_all, npw, evc_init(:,1:nbnd_all), 1, '.wan.occ.emp' )
        !
     ENDIF
     !
     IF (nspin == 2) THEN
        !
        nbnd_all2 =  nbnd_all/2
        CALL read_orbs_spin ( nbnd_all2, npw, evc_init(:,1:nbnd_all2), 1, '.wan.occ.emp' )
        CALL read_orbs_spin ( nbnd_all2, npw, evc_init(:,nbnd_all2+1:nbnd_all), 2, '.wan.occ.emp' )
        !
     ENDIF
     !
  ELSE
     !
     IF (nspin == 1) THEN
        !
        CALL read_orbs_spin ( nbnd_all, npw, evc_init(:,1:nbnd_all), 1, '.kc.occ.emp' )
        !
     ENDIF
     !
     IF (nspin == 2) THEN
        !
        nbnd_all2 =  nbnd_all/2
        CALL read_orbs_spin ( nbnd_all2, npw, evc_init(:,1:nbnd_all2), 1, '.kc.occ.emp' )
        CALL read_orbs_spin ( nbnd_all2, npw, evc_init(:,nbnd_all2+1:nbnd_all), 2, 'kc.occ.emp' )
        !
     ENDIF
     !
  ENDIF 
  !
  IF (nspin == 2) THEN
     !
     CALL ispin_index (fixed_band_i, ispin)
     !
     ! for orb_i
     !  
     IF (ispin == 1) THEN
        nbnd_val_i = nbnd_all/2
        sorted_fixed_band_i =  fixed_band_i
     ENDIF
     !
     IF (ispin == 2) THEN
        nbnd_val_i = nbnd_all/2
        sorted_fixed_band_i = fixed_band_i - nbnd_all/2
     ENDIF
     !
  ENDIF 
  !
  IF (nspin == 1) THEN
     !
     ispin = 1
     !
     nbnd_all   = nbnd_all
     !
     nbnd_val_i = nbnd_all
     !
     sorted_fixed_band_i = fixed_band_i
     !
  ENDIF
  !
  ! here is the main part
  !
  ALLOCATE ( aux_rho  (dfftp%nnr))
  !
  IF ( gamma_only ) THEN
     !
     ! orb 1: G -> R
     !
     psic(:) = ( 0.D0, 0.D0 )
     !
     IF ( MOD(sorted_fixed_band_i,2)==0 ) THEN
        !
        psic(nls(1:npw))  = evc_init(1:npw,fixed_band_i-1) + &
                            ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i)
        psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_i-1) - &
                            ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i))
        !
     ELSE
        !
        IF (sorted_fixed_band_i < nbnd_val_i) THEN 
           ! 
           psic(nls(1:npw))  = evc_init(1:npw,fixed_band_i) + &
                               ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i+1)
           psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_i) - &
                               ( 0.D0, 1.D0 ) * evc_init(1:npw,fixed_band_i+1))
        ELSE
           !
           psic(nls (1:npw))  = evc_init(1:npw,fixed_band_i)
           psic(nlsm(1:npw)) = CONJG( evc_init(1:npw,fixed_band_i) ) 
           !    
        ENDIF
        !
     ENDIF
     !
     CALL invfft ('Wave', psic, dffts)
     !
     ! compute rho(:) = evc_i(:) * evc_i(:)
     !
     aux_rho(:) = (0.0_DP, 0.0_DP)
     !   
     IF ((MOD(sorted_fixed_band_i, 2)==1)) THEN
        !  
        aux_rho(:)   = sa1 * ( DBLE(psic(:))  * DBLE(psic(:)) )
        !
     ENDIF
     !  
     IF ((MOD(sorted_fixed_band_i, 2)==0)) THEN
        ! 
        aux_rho(:)   = sa1 * ( AIMAG(psic(:)) * AIMAG(psic(:)) )
        !
     ENDIF
     !
  ELSE
     !
     ! orb 1: G->R
     ! 
     psic(:) = ( 0.D0, 0.D0 )
     !
     psic(nls(igk_k(1:npw,1))) = evc_init(1:npw,fixed_band_i)
     !
     CALL invfft ('Wave', psic, dffts)
     !
     ! compute rho(:) = evc_i(:) * evc_i(:)
     !
     aux_rho(:) = (0.0_DP, 0.0_DP)
     aux_rho(:) = sa1 * (DBLE(psic(:)) * DBLE(psic(:)) + AIMAG(psic(:)) * AIMAG(psic(:)))
     !
  ENDIF
  !
  drhor(:) = 0.0_DP
  drhor(:) = DBLE(aux_rho(:))
  !
  ! R->G 
  !
  CALL fwfft ('Dense', aux_rho, dfftp)
  !
  drhog(:) = (0.0_DP, 0.0_DP)
  drhog(:) = aux_rho(nl(:))
  !
  DEALLOCATE ( evc_init ) 
  DEALLOCATE ( aux_rho  )
  !
  CALL stop_clock( 'orbital_density' )
  !
  RETURN
  !
ENDSUBROUTINE
