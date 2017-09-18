!
! Copyright (C) 2002-2009 Quantm ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Written and Revised by Carlo Cavazzoni

!=----------------------------------------------------------------------------------=!


SUBROUTINE kcp_runcp_uspp_x &
      ( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec_bgrp, c0_bgrp, cm_bgrp, fromscra, restart )
      !
      !  This subroutine performs a Car-Parrinello or Steepest-Descent step
      !  on the electronic variables, computing forces on electrons
      ! 
      !  on input:
      !  c0_bgrp  wave functions at time t
      !  cm_bgrp  wave functions at time t - dt 
      !
      !  on output:
      !  cm_bgrp  wave functions at time t + dt, not yet othogonalized 
      !
      USE parallel_include
      USE kinds,               ONLY : DP
      USE mp_global,           ONLY : me_bgrp, &
                                      my_bgrp_id, nbgrp, inter_bgrp_comm
      USE mp,                  ONLY : mp_sum
      USE fft_base,            ONLY : dffts, dtgs
      USE fft_parallel,        ONLY : tg_gather
      use wave_base,           only : wave_steepest, wave_verlet
      use control_flags,       only : lwf, tsde
      use uspp,                only : deeq, vkb
      use gvect,  only : gstart
      use electrons_base,      only : nbsp_bgrp, nbspx_bgrp, ispin_bgrp, f_bgrp, nspin, nupdwn_bgrp, iupdwn_bgrp
      use wannier_subroutines, only : ef_potential
      use efield_module,       only : dforce_efield, tefield, dforce_efield2, tefield2
      use gvecw,               only : ngw, ngwx
      USE cp_interfaces,       ONLY : dforce
      USE ldaU_cp,             ONLY : lda_plus_u, vupsi
      USE nksic,               ONLY : do_orbdep, vsic, deeq_sic, vsicpsi
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nfi
      REAL(DP) :: fccc, ccc
      REAL(DP) :: ema0bg(:), dt2bye
      REAL(DP) :: rhos(:,:)
      REAL(DP) :: bec_bgrp(:,:)
      COMPLEX(DP) :: c0_bgrp(:,:), cm_bgrp(:,:)
      LOGICAL, OPTIONAL, INTENT(IN) :: fromscra
      LOGICAL, OPTIONAL, INTENT(IN) :: restart
      !
      !
     real(DP) ::  verl1, verl2, verl3
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: emadt2, emaver, c2, c3, c2tmp, c3tmp, tg_rhos, ftmp, itmp
#endif
#endif
     real(DP),    allocatable :: emadt2(:)
     real(DP),    allocatable :: emaver(:)
     complex(DP), allocatable :: c2(:), c3(:), c2tmp(:), c3tmp(:)
     REAL(DP),    ALLOCATABLE :: tg_rhos(:,:), ftmp(:)
     INTEGER,     ALLOCATABLE :: itmp(:)
     integer :: i, nsiz, incr, idx, idx_in, ierr
     integer :: iwfc, nwfc, is, ii, tg_rhos_siz, c2_siz
     integer :: iflag
     logical :: ttsde
     !
     iflag = 0
     !
     IF ( PRESENT( fromscra ) ) THEN
        !
        IF ( fromscra ) iflag = 1
        ! 
     ENDIF
     ! 
     IF ( PRESENT( restart ) ) THEN
        !
        IF ( restart ) iflag = 2
        !
     ENDIF
     !
     IF ( dtgs%have_task_groups ) THEN
        ! 
        tg_rhos_siz = dtgs%nogrp * dtgs%tg_nnr
        c2_siz      = dtgs%nogrp * ngwx
        !
     ELSE
        !
        tg_rhos_siz = 1
        c2_siz      = ngw 
        !
     ENDIF
     !
     ! ...  set verlet variables 
     !
     verl1 = 2.0d0 * fccc
     verl2 = 1.0d0 - verl1
     verl3 = 1.0d0 * fccc
     !
     ALLOCATE( emadt2( ngw ) )
     ALLOCATE( emaver( ngw ) )
     !
     ccc    = fccc * dt2bye
     emadt2 = dt2bye * ema0bg
     emaver = emadt2 * verl3
     !
     IF ( iflag == 0 ) THEN
        !
        ttsde  = tsde
        !
     ELSEIF( iflag == 1 ) THEN
        !
        ttsde = .TRUE.
        !
     ELSEIF( iflag == 2 ) THEN
        !
        ttsde = .FALSE.
        !
     ENDIF
     !
!============================================================================
! Lingzhu Kong
!     IF( lwf ) THEN
     IF ( .false. ) THEN
        ! 
        call ef_potential( nfi, rhos, bec_bgrp, deeq, vkb, c0_bgrp, cm_bgrp,&
                           emadt2, emaver, verl1, verl2 )
        !
     ELSE
        !
        allocate( c2( c2_siz ), c3( c2_siz ) )
        allocate( tg_rhos( tg_rhos_siz, nspin ) )
        !
        c2      = 0D0
        c3      = 0D0
        ! 
        IF ( dtgs%have_task_groups ) THEN
           !
           !  The potential in rhos is distributed across all processors
           !  We need to redistribute it so that it is completely contained in the
           !  processors of an orbital TASK-GROUP
           !
           DO i = 1, nspin
              !
              CALL tg_gather( dffts, dtgs, rhos(:,i), tg_rhos(:,i) )
              !
           ENDDO
           !
           incr = 2 * dtgs%nogrp
           !
        ELSE
           !
           incr = 2
           !
        ENDIF
        ! 
        DO i = 1, nbsp_bgrp, incr
           !
           IF ( dtgs%have_task_groups ) THEN
              !
              !The input coefficients to dforce cover eigenstates i:i+2*NOGRP-1
              !Thus, in dforce the dummy arguments for c0_bgrp(1,i) and
              !c0_bgrp(1,i+1) hold coefficients for eigenstates i,i+2*NOGRP-2,2
              !and i+1,i+2*NOGRP...for example if NOGRP is 4 then we would have
              !1-3-5-7 and 2-4-6-8
              !
              IF ( tefield .OR. tefield2 ) THEN
                 ! 
                 CALL errore( ' runcp_uspp ', ' electric field with task group not implemented yet ', 1 )
                 !
              ENDIF
              !
              IF ( nspin > 1 .AND. ispin_bgrp(i) /= ispin_bgrp( MIN( nbsp_bgrp, i+incr-1 ) ) ) THEN
                 !   
                 ! when computing force with task group and states with different spin,
                 ! we need to compute spin up and spin down separately because the logics 
                 ! of computing two states with different spin at the same time do not work any longer
                 !
                 ALLOCATE( c2tmp( c2_siz ) )
                 ALLOCATE( c3tmp( c2_siz ) )
                 ALLOCATE( ftmp( nbsp_bgrp ) )
                 ALLOCATE( itmp( nbsp_bgrp ) )
                 !
                 !  spin up
                 !
                 itmp = ispin_bgrp(i)
                 ftmp = f_bgrp(i)
                 c2tmp = 0.0d0
                 c3tmp = 0.0d0
                 !
                 CALL dforce( i, bec_bgrp, vkb, c0_bgrp, c2tmp, c3tmp, tg_rhos, tg_rhos_siz, itmp, ftmp, nbsp_bgrp, nspin )
                 ! 
                 idx_in = 1
                 ! 
                 DO idx = 1, incr, 2
                    ! 
                    IF ( i + idx - 1 <= nbsp_bgrp ) THEN
                       !
                       IF ( ispin_bgrp( i + idx - 1 ) == ispin_bgrp(i) ) THEN
                          !
                          c2( (idx_in-1)*ngw+1 : idx_in*ngw ) = c2tmp( (idx_in-1)*ngw+1 : idx_in*ngw )
                          ! 
                       ENDIF
                       !
                       IF ( ispin_bgrp( i + idx     ) == ispin_bgrp(i) ) THEN
                          !
                          c3( (idx_in-1)*ngw+1 : idx_in*ngw ) = c3tmp( (idx_in-1)*ngw+1 : idx_in*ngw )
                          !
                       ENDIF
                       ! 
                    ENDIF
                    !
                    idx_in = idx_in + 1
                    !
                 ENDDO
                 !
                 !  spin down
                 !
                 itmp = ispin_bgrp( MIN( nbsp_bgrp, i+incr-1 ) )
                 ftmp = f_bgrp( MIN( nbsp_bgrp, i+incr-1 ) )
                 c2tmp = 0.0d0
                 c3tmp = 0.0d0
                 ! 
                 CALL dforce( i, bec_bgrp, vkb, c0_bgrp, c2tmp, c3tmp, tg_rhos, tg_rhos_siz, itmp, ftmp, nbsp_bgrp, nspin )
                 !
                 idx_in = 1
                 ! 
                 DO idx = 1, incr, 2
                    !
                    IF ( i + idx - 1 <= nbsp_bgrp ) THEN
                       !
                       IF ( ispin_bgrp( i + idx - 1 ) == ispin_bgrp( MIN( nbsp_bgrp, i+incr-1 ) ) ) THEN
                          !
                          c2( (idx_in-1)*ngw+1 : idx_in*ngw ) = c2tmp( (idx_in-1)*ngw+1 : idx_in*ngw )
                          !
                       ENDIF
                       !
                       IF ( ispin_bgrp( i + idx     ) == ispin_bgrp( MIN( nbsp_bgrp, i+incr-1 ) ) ) THEN
                          !
                          c3( (idx_in-1)*ngw+1 : idx_in*ngw ) = c3tmp( (idx_in-1)*ngw+1 : idx_in*ngw )
                          ! 
                       ENDIF
                       !
                    ENDIF 
                    !
                    idx_in = idx_in + 1
                    !
                 ENDDO
                 !
                 DEALLOCATE( itmp )
                 DEALLOCATE( ftmp )
                 DEALLOCATE( c3tmp )
                 DEALLOCATE( c2tmp )
                 !
              ELSE
                 !
                 CALL dforce( i, bec_bgrp, vkb, c0_bgrp, c2, c3, tg_rhos, tg_rhos_siz, ispin_bgrp, f_bgrp, nbsp_bgrp, nspin )
                 !
              ENDIF
              !
              IF ( lda_plus_u ) THEN
                 ! 
                 idx_in = 1
                 !
                 DO idx = 1, incr, 2
                    !
                    ii = i+idx-1
                    !
                    IF ( ii <= nbsp_bgrp ) THEN
                       ! 
                       c2( (idx_in-1)*ngw+1 : idx_in*ngw ) = &
                       c2( (idx_in-1)*ngw+1 : idx_in*ngw ) - vupsi(1:ngw,ii)
                       c3( (idx_in-1)*ngw+1 : idx_in*ngw ) = &
                       c3( (idx_in-1)*ngw+1 : idx_in*ngw ) - vupsi(1:ngw,ii+1)
                       !
                    ENDIF
                    !
                    idx_in = idx_in + 1
                    !
                 ENDDO
                 !
              ENDIF
              ! 
           ELSE
              !
              CALL dforce( i, bec_bgrp, vkb, c0_bgrp, c2, c3, rhos, &
                           SIZE(rhos,1), ispin_bgrp, f_bgrp, nbsp_bgrp, nspin )
              ! 
              IF ( lda_plus_u ) THEN
                 !
                 c2(:) = c2(:) - vupsi(:,i)
                 c3(:) = c3(:) - vupsi(:,i+1)
                 ! 
              ENDIF
              !
              IF ( do_orbdep ) THEN
                 !
                 ! faux takes into account spin multiplicity.
                 !
                 CALL nksic_eforce( i, nbsp_bgrp, nbspx_bgrp, vsic, deeq_sic, bec_bgrp, ngw, c0_bgrp(:,i), c0_bgrp(:,i+1), vsicpsi)
                 !
                 c2(:) = c2(:) - vsicpsi(:,1) * f_bgrp(i)
                 c3(:) = c3(:) - vsicpsi(:,2) * f_bgrp(i+1)
                 !
              ENDIF
              !
           ENDIF
           !
           IF ( tefield ) THEN
              !
              CALL dforce_efield ( bec_bgrp, i, c0_bgrp, c2, c3, rhos)
              !
           ENDIF
           !
           IF ( tefield2 ) THEN
              !
              CALL dforce_efield2 ( bec_bgrp, i, c0_bgrp, c2, c3, rhos)
              ! 
           ENDIF
           ! 
           IF ( iflag == 2 ) THEN
              ! 
              DO idx = 1, incr, 2
                 !
                 IF ( i + idx - 1 <= nbsp_bgrp ) THEN
                    !
                    cm_bgrp( :, i+idx-1) = c0_bgrp(:,i+idx-1)
                    cm_bgrp( :, i+idx  ) = c0_bgrp(:,i+idx  )
                    ! 
                 ENDIF
                 !
              ENDDO
              ! 
           ENDIF
           !
           idx_in = 1
           !
           DO idx = 1, incr, 2
              !
              IF ( i + idx - 1 <= nbsp_bgrp ) THEN
                 ! 
                 IF ( tsde ) THEN
                    ! 
                    CALL wave_steepest( cm_bgrp(:, i+idx-1 ), c0_bgrp(:, i+idx-1 ), emaver, c2(:), ngw, idx_in )
                    CALL wave_steepest( cm_bgrp(:, i+idx   ), c0_bgrp(:, i+idx   ), emaver, c3(:), ngw, idx_in )
                    !
                 ELSE
                    !
                    CALL wave_verlet( cm_bgrp(:, i+idx-1 ), c0_bgrp(:, i+idx-1 ), verl1, verl2, emaver, c2(:), ngw, idx_in )
                    CALL wave_verlet( cm_bgrp(:, i+idx   ), c0_bgrp(:, i+idx   ), verl1, verl2, emaver, c3(:), ngw, idx_in )
                    !
                 ENDIF
                 ! 
                 IF ( gstart == 2 ) THEN
                    !
                    cm_bgrp(1,i+idx-1) = CMPLX(real(cm_bgrp(1,i+idx-1)),0.0d0,kind=dp)
                    cm_bgrp(1,i+idx  ) = CMPLX(real(cm_bgrp(1,i+idx  )),0.0d0,kind=dp)
                    !
                 ENDIF
                 !
              ENDIF
              !
              idx_in = idx_in + 1
              !
           ENDDO
           ! 
        ENDDO
        !
        DEALLOCATE( c2 )
        DEALLOCATE( c3 )
        DEALLOCATE( tg_rhos )
        !
     ENDIF
     !
     DEALLOCATE( emadt2 )
     DEALLOCATE( emaver )
     !
ENDSUBROUTINE kcp_runcp_uspp_x
