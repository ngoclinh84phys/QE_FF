!-----------------------------------------------------------------------
SUBROUTINE empty_cp ( nfi, c0, v, tcg )
!-----------------------------------------------------------------------
      !
      ! Performs the minimization on the empty state subspace keeping the
      ! occupied manyfold fixed. A proper orthogonalization of the two
      ! manyfolds is performed. 
      !
      USE kinds,                ONLY : DP
      USE constants,            ONLY : autoev
      USE control_flags,        ONLY : iprsta, tsde, program_name, gamma_only, tortho
      USE io_global,            ONLY : ionode, stdout
      USE cp_main_variables,    ONLY : eigr, ema0bg, collect_lambda, &
                                       rhor, rhog, rhos, eigr, eigrb, irb, bec, bec_emp
      USE descriptors,          ONLY : descla_siz_ , descla_init, nlax_, lambda_node_
      USE cell_base,            ONLY : omega
      USE uspp,                 ONLY : vkb, nkb, okvan
      USE uspp_param,           ONLY : nhm
      USE grid_dimensions,      ONLY : nnrx
      USE electrons_base,       ONLY : nbsp, nbspx, ispin, nspin, f, nudx, iupdwn, nupdwn
      USE electrons_module,     ONLY : iupdwn_emp, nupdwn_emp, n_emp, ei_emp,  &
                                       max_emp, ethr_emp, etot_emp, eodd_emp
      USE ions_base,            ONLY : nat, nsp
      USE gvecw,                ONLY : ngw
      USE orthogonalize_base,   ONLY : calphi, updatc
      USE reciprocal_vectors,   ONLY : gzero, gstart
      USE wave_base,            ONLY : wave_steepest, wave_verlet, frice
      USE cvan,                 ONLY : nvb
      USE cp_electronic_mass,   ONLY : emass
      USE time_step,            ONLY : delt
      USE check_stop,           ONLY : check_stop_now
      USE cp_interfaces,        ONLY : writeempty, readempty, gram_empty, ortho, &
                                       wave_rand_init, wave_atom_init, elec_fakekine, &
                                       crot, dforce, nlsm1, grabec, &
                                       bec_csv, readempty_twin, writeempty_twin
      USE mp,                   ONLY : mp_comm_split, mp_comm_free, mp_sum
      USE mp_global,            ONLY : intra_image_comm, me_image
      USE nksic,                ONLY : do_orbdep, do_pz, do_wxd, vsicpsi, wtot, sizwtot, &
                                       odd_alpha, valpsi, nkscalfact, odd_alpha_emp
      USE nksic,                ONLY : do_spinsym, pink_emp, allocate_nksic_empty
      USE hfmod,                ONLY : do_hf, vxxpsi
      USE control_flags,        ONLY : tatomicwfc, trane, ndr, ndw
      USE electrons_module,     ONLY : wfc_centers_emp, wfc_spreads_emp, icompute_spread
      USE core,                 ONLY : nlcc_any, rhoc
      USE input_parameters,     ONLY : odd_nkscalfact_empty,  &
                                       restart_from_wannier_cp, wannier_empty_only, &
                                       fixed_band, print_wfc_anion, wo_odd_in_empty_run, &
                                       odd_nkscalfact, index_empty_to_save
      USE wavefunctions_module, ONLY : c0fixed_emp
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: nfi
      COMPLEX(DP)            :: c0(:,:)
      REAL(DP)               :: v(:,:)
      logical, optional, intent(IN) :: tcg
      !
      INTEGER  :: i, iss, j, in, in_emp, iter, iter_ortho
      INTEGER  :: n_occs, n_emps, n_empx, nudx_emp, issw, n
      INTEGER  :: nlax_emp, nlam_emp
      LOGICAL  :: exst, do_wxd_, tcg_
      !
      REAL(DP) :: fccc, ccc, csv, dt2bye, bigr
      REAL(DP) :: verl1, verl2, verl3
      REAL(DP) :: dek, ekinc, ekinc_old, detothf
      !
      REAL(DP),    ALLOCATABLE :: emadt2(:)
      REAL(DP),    ALLOCATABLE :: emaver(:)
      COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
      COMPLEX(DP), ALLOCATABLE :: c0_emp(:,:), cm_emp(:,:), phi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: becsum_emp(:,:,:)
      type(twin_matrix) :: bephi_emp! !modified:giovanni
      type(twin_matrix) :: becp_emp !modified:giovanni
      type(twin_matrix) :: bec_occ !(:,:) !modified:giovanni
      type(twin_matrix), dimension(:),  ALLOCATABLE :: lambda_emp !(:,:,:) !, 
      REAL(DP),    ALLOCATABLE :: f_emp(:)
      REAL(DP),    ALLOCATABLE :: lambda_rep(:,:)
      COMPLEX(DP), ALLOCATABLE :: lambda_rep_c(:,:)
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      REAL(DP),    ALLOCATABLE :: fsic_emp(:)
      REAL(DP),    ALLOCATABLE :: vsic_emp(:,:)
      REAL(DP),    ALLOCATABLE :: wxd_emp(:,:)
      REAL(DP),    ALLOCATABLE :: deeq_sic_emp(:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: vxxpsi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: exx_emp(:)
      REAL(DP),    ALLOCATABLE :: old_odd_alpha(:)
      !
      INTEGER, SAVE :: np_emp(2), me_emp(2), emp_comm, color
      INTEGER, SAVE :: desc_emp( descla_siz_ , 2 )
      LOGICAL, SAVE :: first = .true.
      LOGICAL :: lgam !added:giovanni
      LOGICAL :: done_extra !added:giovanni
      COMPLEX(DP), PARAMETER :: c_zero=CMPLX(0.d0,0.d0)
      INTEGER :: sizvsic_emp
      INTEGER :: ndr_loc, ndw_loc
      !
      LOGICAL :: odd_nkscalfact_old
      INTEGER :: nbnd_, ib, start_is
      COMPLEX(DP), ALLOCATABLE :: c0_anion(:,:)
      !
      if(nbgrp > 1) &
        call errore(' cp_empty ', ' parallelization over bands not yet implemented ', 1 )
      !
      if(present(tcg)) THEN
         !
         tcg_=tcg
         !
      ELSE
         !
         tcg_=.false.
         !
      ENDIF
      !
      ! ...  quick exit if empty states have not to be computed
      !
      f_emp      = 2.0d0 / DBLE(nspin)
      !
      ! ... setup main variables for cp empty run
      !
      call jsahjfah 
      !   
      call prefor( eigr, vkb )
      !
      call nlsm1 ( n_occs, 1, nvb, eigr, c0, bec_occ, 1)
      !
      ! here is initialize wfcs
      !
      exst = readempty( c0_emp, n_empx, ndr_loc )
      !
      if (.not. exst ) then
         !
         write(stdout, * ) 'Linh: oopp restart from minimizing orbital does not work for emptystate'
         write(stdout, * ) 'Linh: initialize random states and orthogonalize to filled ones'
         !
         ! ...  initial random states orthogonal to filled ones
         !
         if ( .not.do_spinsym .or. nspin == 1 ) then
            !
            call wave_rand_init( c0_emp, n_emps, 1 )
            !
         else 
            !  
            if ( nupdwn_emp(1) < nupdwn_emp(2) ) &
               cell errore('empty_cp','unexpec emp nupdwn(1) < nupdwn(2)',10)
            !
            write(stdout, "(24x, 'spin symmetry applied to init wave')" )
            !
            call wave_rand_init( c0_emp, nupdwn_emp(1) , 1 )
            !
            do i = 1, min(nupdwn_emp(1),nupdwn_emp(2))
               !
               j = i+iupdwn_emp(2)-1
               c0_emp(:,j) = c0_emp(:,i)
               !
            enddo
            !
         endif 
         !
         if (gzero) c0_emp( 1, : ) = (0.0d0, 0.0d0)
         !
         call nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp, 1, lgam )
         !
         do iss = 1, nspin
            !
            in_emp = iupdwn_emp(iss)
            !
            issw   = iupdwn(iss)
            !
            if (nupdwn(iss)>0.and.nupdwn_emp(iss)>0) then
               !
               call gram_schmidt_empty ( .false., eigr, vkb, bec_emp, bec_occ, nkb, &
                                          c0_emp( :, in_emp: ), c0( :, issw: ), ngw,&
                                          nupdwn_emp(iss), nupdwn(iss), in_emp, issw )
               !
            endif
            !
         enddo
         !
      else
         !
         write(stdout, * ) 'Linh: the code restarts not random wfc'
         !
      endif
      !
      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp, 1, lgam )
      !
      cm_emp = c0_emp
      !  
      if (tcg_) then ! compute empty states with conjugate gradient
         ! 
         !call runcg_uspp_emp(c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
         !                    n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
         !                    max_emp, wxd_emp, vsic_emp, sizvsic_emp, pink_emp, nnrx, becsum_emp, &
         !                    deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, v, &
         !                    nfi, .true., .true., eigr, bec, irb, eigrb, &
         !                    rhor, rhog, rhos, rhoc, ema0bg, desc_emp)   
         !
      else ! compute empty states with damped dynamics
         !
         ! ...  set verlet variables
         !
         if ( tsde ) then
            !
            fccc = 1.0d0
            ! 
         else
            ! 
            fccc = 1.0d0/(1.0d0+frice)
            !
         endif
         !
         verl1 = 2.0d0 * fccc
         verl2 = 1.0d0 - verl1
         verl3 = 1.0d0 * fccc
         !
         allocate( c2( ngw ) )
         allocate( c3( ngw ) )
         allocate( emadt2( ngw ) )
         allocate( emaver( ngw ) )
         !  
         dt2bye = delt * delt / emass
         !
         ccc    = fccc   * dt2bye
         emadt2 = dt2bye * ema0bg
         emaver = emadt2 * verl3
         ! 
         ekinc_old = 0.0
         ekinc     = 0.0
         !     
         done_extra=.false.
         !
         if ( ionode ) then
            ! 
            write ( stdout, "(/,3X,'Empty states damping/stepest-decent minimization starting ...', &
                            & /,3x,'nfi         dekinc         ekinc' )")
            !   
         endif
         !
         ITERATIONS: do iter = 1, max_emp
            !
            ! ... standard terms
            !
            do i = 1, n_emps, 2
               !
               call dforce( i, bec_emp, vkb, c0_emp, c2, c3, vloc, SIZE(vloc,1), ispin_emp, f_emp, n_emps, nspin )
               !
               if ( tsde ) then
                  !
                  call wave_steepest( cm_emp(:, i  ), c0_emp(:, i  ), emaver, c2 )
                  ! 
                  call wave_steepest( cm_emp(:, i+1), c0_emp(:, i+1), emaver, c3 )
                  !
               else
                  ! 
                  call wave_verlet( cm_emp(:, i  ), c0_emp(:, i  ), verl1, verl2, emaver, c2 )
                  call wave_verlet( cm_emp(:, i+1), c0_emp(:, i+1), verl1, verl2, emaver, c3 )
                  !  
               endif
               !
               if ( gstart == 2) then
                  ! 
                  cm_emp(1,  i)=CMPLX(DBLE(cm_emp(1,  i)),0.d0)
                  cm_emp(1,i+1)=CMPLX(DBLE(cm_emp(1,i+1)),0.d0)
                  ! 
               endif
               !
            enddo
            ! 
            ! ortho cm_emp with c0
            !     
            do iss = 1, nspin
               !
               in     = iupdwn(iss)
               in_emp = iupdwn_emp(iss)
               !
               issw   = iupdwn( iss )
               !
               if (nupdwn(iss)>0 .and. nupdwn_emp(iss)>0) then
                  !
                  call gram_schmidt_empty( .true. , eigr, vkb, bec_emp, bec_occ, nkb, &
                                   cm_emp( :, in_emp: ), c0( :, issw: ), ngw, &
                                   nupdwn_emp(iss), nupdwn(iss), in_emp, issw )
                  !
               endif
               !
            enddo
            !
            ! ... calphi calculates phi
            ! ... the electron mass rises with g**2
            !
            call calphi_bgrp_empty ( c0_emp, ngw, bec_emp, nkb, vkb, phi_emp, n_emps, ema0bg )
            !
            if ( tortho ) then
               !
               call ortho_empty( eigr(1:ngw,1:nat), cm_emp(1:ngw,1:n_emps), phi_emp(1:ngw,1:n_emps), &
                                   ngw, lambda_emp(1:nspin), desc_emp(1:descla_siz_,1:nspin), &
                                   bigr, iter_ortho, ccc, bephi_emp, becp_emp, n_emps, nspin, &
                                   nupdwn_emp, iupdwn_emp )
               !
            else
               !
               call gram_empty( vkb, bec_emp, nkb, cm_emp, ngw, n_emps )
               !
            endif
            !
            do iss = 1, nspin
               !
               CALL updatc_empty( ccc, n_emps, lambda_emp(:,:,iss), n_empx, phi_emp, ngw, &
                                  bephi_emp, nkb, becp_emp, bec_emp, &
                                   cm_emp, nupdwn_emp(iss), iupdwn_emp(iss), desc_emp(:,iss) )
               !
            enddo
            !
            call nlsm1 ( n_emps, 1, nsp, eigr, cm_emp, bec_emp, 1 )
            !
            call elec_fakekine( ekinc, ema0bg, emass, c0_emp, cm_emp, ngw, n_emps, 1, delt )
            !
            call dswap( 2*SIZE( c0_emp ), c0_emp, 1, cm_emp, 1 )
            !
            dek = ekinc - ekinc_old
            !
            write(stdout, 112) iter, ekinc, dek
            !    
            ! ...   check for exit
            !
            if ( check_stop_now() ) then
               !
               exit ITERATIONS
               !
            endif
            !  
            ! ...  check for convergence
            !     
            if ( ( ekinc <  ethr_emp ) .AND. ( iter > 3 ) ) then
               !
               if (done_extra) then
                  !
                  if ( ionode ) write( stdout, 113)
                  !  
                  exit ITERATIONS
                  !  
               else
                  ! 
                  done_extra=.true.
                  !
               endif
               !
            endif
            ! 
            ekinc_old = ekinc
            !
         enddo ITERATIONS
         ! 
         deallocate( c2 )
         deallocate( c3 )
         deallocate( emadt2 )
         deallocate( emaver )
         !  
      endif !if clause to choose between cg and damped dynamics
      !
      ! ...  Compute eigenvalues and bring wave functions on Kohn-Sham orbitals
      !
      ALLOCATE( lambda_rep( nudx_emp, nudx_emp ) )
      !
      DO iss = 1, nspin
         !
         i = iupdwn_emp(iss)
         n = nupdwn_emp(iss)
         !
         CALL collect_lambda( lambda_rep, lambda_emp(:,:, iss), desc_emp( :, iss ) )
         CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep, nudx_emp, ei_emp(:,iss) )
         !   
         ei_emp( 1:n, iss ) = ei_emp( 1:n, iss ) / f_emp( i : i + n - 1 )
         !
      ENDDO
      !
      DEALLOCATE( lambda_rep)
      ! 
      DEALLOCATE( ispin_emp )
      DEALLOCATE( f_emp )
      DEALLOCATE( c0_emp )
      DEALLOCATE( cm_emp )
      DEALLOCATE( phi_emp )
      !
      DEALLOCATE( bec_emp )
      DEALLOCATE( bec_occ )
      DEALLOCATE( bephi_emp )
      DEALLOCATE( becp_emp )
      !
      DO iss=1,nspin
         CALL deallocate_twin(lambda_emp(iss))
      ENDDO
      !  
112   FORMAT(I5,2X,2D14.6)
113   FORMAT(/,3X,'Empty states: convergence achieved')
      !
      return 
      !
endsubroutine empty_cp_twin_x
