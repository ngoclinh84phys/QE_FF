SUBROUTINE empty_cp ( nfi, eigr, eigrb, rhor, rhog, rhos, rhoc, &
                      vpot, c0, cm, bec_occ, phi )
      !
      ! Performs the minimization of the effective emptystates energy keeping the
      ! occupied manyfold fixed. A proper orthogonalization of the two
      ! manyfolds is performed. 
      !
      use kinds,          only: DP
      use gvecw,          only: ngw
      use io_global,      only: stdout
      use ions_base,      only: na, nat, nax, nsp, rcmax
      use fft_base,       only: dffts, dfftp
      use gvect,          only: ngm, gstart
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt,&
                                nbspx, nbsp, ispin
      use control_flags,  only: ndw, ndr
      use mp_global,      only: me_image, my_image_id, nbgrp, intra_bgrp_comm
      use descriptors,    only: la_descriptor, descla_init
      USE mp,             only: mp_rank, mp_comm_free, mp_comm_split, mp_sum
      use uspp,           only: vkb, nkb 
      use uspp_param,     only: nvb, nhm
      use smallbox_gvec,  only: ngb
      use cp_interfaces,  only: nlsm1, collect_lambda, crot
      use emptystates_electrons_module, only: iupdwn_emp, nupdwn_emp, & 
                                n_emp, ei_emp, max_emp 
      use nksic,          only : do_orbdep, pink_emp, allocate_nksic_empty
      use cp_main_variables, only: irb, ema0bg
      use cell_base, only : omega 
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN)  :: nfi 
      complex(dp) :: eigr(ngw,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp) :: rhor(dfftp%nnr,nspin)
      real(dp) :: vpot(dfftp%nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp) :: rhos(dffts%nnr,nspin)
      real(dp) :: rhoc(dfftp%nnr)
      complex(dp) :: c0( ngw, nbspx )
      complex(dp) :: cm( ngw, nbspx )
      complex(dp) :: phi( ngw, nbspx )
      real(dp)    :: bec_occ(nkb,nbspx)
      !
      ! local variables
      !
      INTEGER  :: i, n, in, iss, j, in_emp, iter 
      INTEGER  :: n_occs  
      INTEGER  :: n_emps, n_empx, nudx_emp
      INTEGER  :: issw  
      INTEGER  :: nlam_emp
      !
      COMPLEX(DP), ALLOCATABLE :: c0_emp(:,:), cm_emp(:,:), phi_emp(:,:)
      REAL(DP),    ALLOCATABLE :: becsum_emp(:,:,:)
      REAL(DP),    ALLOCATABLE :: bec_emp(:,:) 
      REAL(DP),    ALLOCATABLE :: f_emp(:)
      REAL(DP),    ALLOCATABLE :: lambda_rep(:,:)
      REAL(DP),    ALLOCATABLE :: lambda_emp(:,:,:)
      INTEGER,     ALLOCATABLE :: ispin_emp(:)
      !
      INTEGER :: sizvsic_emp
      REAL(DP) :: etot_emp, eodd_emp
      REAL(DP),    ALLOCATABLE :: fsic_emp(:)
      REAL(DP),    ALLOCATABLE :: vsic_emp(:,:)
      REAL(DP),    ALLOCATABLE :: wxd_emp(:,:)
      REAL(DP),    ALLOCATABLE :: deeq_sic_emp(:,:,:,:)
      !
      INTEGER, SAVE :: np_emp(2), me_emp(2), emp_comm, color
      TYPE(la_descriptor), ALLOCATABLE :: descla_emp(:)
      LOGICAL, SAVE :: first = .true.
      INTEGER :: ndr_loc, ndw_loc
      !
      !
      if (nbgrp > 1) &
         call errore(' cp_empty ', ' parallelization over bands not yet implemented ', 1 )
      !
      ! ...  quick exit if empty states have not to be computed
      !
      if ( n_emp < 1 ) return
      !  
      f_emp = 2.0d0 / DBLE(nspin)
      !
      ! ... setup main variables for cp empty run
      !
      ! restart directories
      !
      IF ( first ) THEN
         ndr_loc = ndr
         ndw_loc = ndw
      ELSE
         ndr_loc = ndw
         ndw_loc = ndw
      ENDIF
      !
      !  Here set the group of processors for empty states
      !
      ALLOCATE( descla_emp( nspin ) )
      !
      IF ( .NOT. first ) THEN
         CALL mp_comm_free( emp_comm )
      ENDIF
      !
      me_image = mp_rank( intra_bgrp_comm ) 
      np_emp = 1
      !
      IF ( me_image < np_emp(1) * np_emp(2) ) THEN
         color = 1
      ELSE
         color = 0
      ENDIF
      !
      CALL mp_comm_split( intra_bgrp_comm, color, me_image, emp_comm )
      !  
      if( me_image <  np_emp(1) * np_emp(2) ) then
          me_emp(1) = me_image / np_emp(1)
          me_emp(2) = MOD( me_image, np_emp(1) )
      else
          me_emp(1) = me_image
          me_emp(2) = me_image
      endif
      !
      first = .FALSE.
      !
      !  Done with the group
      !
      ! n_occs    == nbsp
      ! n_emps    => nbsp   (corresponds to)
      ! n_empx    => nbspx  
      ! nudx_emp  => nudx
      !
      !
      n_occs = nupdwn( 1 )
      if ( nspin == 2 ) n_occs = n_occs + nupdwn( 2 )
      !
      n_emps = nupdwn_emp( 1 )
      if ( nspin == 2 ) n_emps = n_emps + nupdwn_emp( 2 )
      !
      nudx_emp = nupdwn_emp( 1 )
      IF( nspin == 2 ) nudx_emp = MAX( nudx_emp, nupdwn_emp( 2 ) )
      !       
      n_empx = nupdwn_emp( 1 )
      IF( nspin == 2 ) n_empx = n_empx + nupdwn_emp( 2 )
      n_empx = n_empx + MOD( n_empx, 2)
      !
      DO iss = 1, nspin
         CALL descla_init( descla_emp( iss ), nupdwn_emp( iss ), nudx_emp, np_emp, me_emp, emp_comm, -1, color )
      ENDDO
      !
      nlam_emp = 1
      IF ( SIZE( descla_emp ) < 2 ) THEN
         IF ( descla_emp(1)%active_node > 0 ) &
            nlam_emp = descla_emp(1)%nrcx
      ELSE
         IF ( ( descla_emp(1)%active_node > 0 ) .OR. ( descla_emp(2)%active_node > 0 ) ) &
            nlam_emp = MAX( descla_emp(1)%nrcx, descla_emp(2)%nrcx )
      ENDIF
      !
      ALLOCATE( lambda_emp( nlam_emp, nlam_emp, nspin ) )
      !
      ALLOCATE( c0_emp( ngw, n_empx * nspin ) )
      ALLOCATE( cm_emp( ngw, n_empx * nspin ) )
      ALLOCATE( phi_emp( ngw, n_empx * nspin ) )
      ALLOCATE( bec_emp( nkb, n_emps ) )
      ALLOCATE( f_emp( n_empx * nspin ) )
      ALLOCATE( ispin_emp( n_empx * nspin ) )
      ALLOCATE( becsum_emp(nhm*(nhm+1)/2,nat,nspin))
      !  
      ! initializes
      ! 
      c0_emp     = (0.0_dp,0.0_dp)
      cm_emp     = (0.0_dp,0.0_dp) 
      phi_emp    = (0.0_dp,0.0_dp)
      bec_emp    = 0.0_dp
      becsum_emp = 0.0_dp 
      !
      lambda_emp(:,:,:) = 0.0_dp
      !
      f_emp     = 2.0d0 / DBLE(nspin)
      !
      ispin_emp = 0
      ispin_emp( 1:nupdwn_emp( 1 ) ) = 1
      IF( nspin == 2 ) ispin_emp( iupdwn_emp(2) : ) = 2
      !
      IF ( do_orbdep ) THEN
         !
         ALLOCATE( fsic_emp( n_empx ) )
         ! n_empx_odd=n_empx
         ALLOCATE( vsic_emp(dfftp%nnr, n_empx) )
         ALLOCATE( wxd_emp (dfftp%nnr, 2) )
         ALLOCATE( deeq_sic_emp (nhm,nhm,nat,n_empx) )
         ! 
         CALL allocate_nksic_empty(n_empx)
         sizvsic_emp=dfftp%nnr
         !
         fsic_emp = 0.0d0
         vsic_emp = 0.0d0
         wxd_emp  = 0.0d0
         ! 
      ELSE
         !
         ALLOCATE( fsic_emp( n_empx ) )
         ! n_empx_odd=1
         ALLOCATE( vsic_emp(1, n_empx) )
         ALLOCATE( wxd_emp (1, 2) )
         ALLOCATE( deeq_sic_emp (nhm,nhm,nat,n_empx) )
         !
         call allocate_nksic_empty(n_empx)
         sizvsic_emp=1
         !
         fsic_emp = 0.0d0
         vsic_emp = 0.0d0
         wxd_emp  = 0.0d0
         !
      ENDIF
      !
      ! here is initialization for wfcs
      !
      !exst = readempty( c0_emp, n_empx, ndr_loc )
      !
      c0_emp( :, : ) = (1.0d0, 0.0d0)
      ! 
      !if (.not. exst ) then
         !
      !   write(stdout, * ) 'Linh: oopp restart from minimizing orbital does not work for emptystate'
      !   write(stdout, * ) 'Linh: initialize random states and orthogonalize to filled ones'
         !
         ! ...  initial random states orthogonal to filled ones
         !
      !   if ( .not.do_spinsym .or. nspin == 1 ) then
            !
      !      call wave_rand_init( c0_emp, n_emps, 1 )
            !
      !   else 
            !  
      !      if ( nupdwn_emp(1) < nupdwn_emp(2) ) &
      !         call errore('empty_cp','unexpec emp nupdwn(1) < nupdwn(2)',10)
            !
      !      write(stdout, "(24x, 'spin symmetry applied to init wave')" )
            !
      !      call wave_rand_init( c0_emp, nupdwn_emp(1) , 1 )
            !
      !      do i = 1, min(nupdwn_emp(1),nupdwn_emp(2))
               !
      !         j = i+iupdwn_emp(2)-1
      !         c0_emp(:,j) = c0_emp(:,i)
               !
      !      enddo
            !
      !   endif 
         !
      !   if (gzero) c0_emp( 1, : ) = (0.0d0, 0.0d0)
         !
      !   call nlsm1 ( n_emps, 1, nvb, eigr, c0_emp, bec_emp, 1, lgam )
         !
      !   do iss = 1, nspin
            !
      !      in_emp = iupdwn_emp(iss)
            !
      !      issw   = iupdwn(iss)
            !
      !      if (nupdwn(iss)>0.and.nupdwn_emp(iss)>0) then
               !
      !         call gram_schmidt_empty ( .false., eigr, vkb, bec_emp, bec_occ, nkb, &
      !                                    c0_emp( :, in_emp: ), c0( :, issw: ), ngw,&
      !                                    nupdwn_emp(iss), nupdwn(iss), in_emp, issw )
               !
      !      endif
            !
      !   enddo
         !
      !else
         !
      !  write(stdout, * ) 'Linh: the code restarts from the given wfcs'
         !
      !endif
      !
      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp )
      !
      DO iss = 1, nspin
         !
         in     = iupdwn(iss)
         in_emp = iupdwn_emp(iss)
         !
         issw   = iupdwn( iss )
         !
         CALL gram_empty( .false. , eigr, vkb, bec_emp( :, in_emp: ), bec_occ( :, in: ), nkb, &
                           c0_emp( :, in_emp: ), c0( :, issw: ), ngw, nupdwn_emp(iss), nupdwn(iss) )
         !
      ENDDO
      !
      CALL nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp )
      cm_emp = c0_emp
      !  
      call emptystates_cg_sub (c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
                          n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
                          max_emp, wxd_emp, vsic_emp, sizvsic_emp, pink_emp, dfftp%nnr, becsum_emp, &
                          deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, vpot, &
                          nfi, .true., .true., eigr, c0, bec_occ, irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ema0bg, descla_emp, nlam_emp)   
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
         CALL collect_lambda( lambda_rep, lambda_emp(:,:, iss), descla_emp( iss ) )
         CALL crot( cm_emp, c0_emp, ngw, n, i, i, lambda_rep, nudx_emp, ei_emp(:,iss) )
         !   
         ei_emp( 1:n, iss ) = ei_emp( 1:n, iss ) / f_emp( i : i + n - 1 )
         !
      ENDDO
      !
      DEALLOCATE( descla_emp )
      DEALLOCATE( ispin_emp )
      DEALLOCATE( f_emp )
      DEALLOCATE( c0_emp )
      DEALLOCATE( cm_emp )
      DEALLOCATE( phi_emp )
      !
      DEALLOCATE( bec_emp )
      !
      DEALLOCATE( lambda_rep)
      DEALLOCATE (lambda_emp)
      DEALLOCATE (becsum_emp)
      !  
      DEALLOCATE( fsic_emp )
      DEALLOCATE( vsic_emp )
      DEALLOCATE( wxd_emp )
      DEALLOCATE( deeq_sic_emp )
      !
      return 
      !
endsubroutine empty_cp
