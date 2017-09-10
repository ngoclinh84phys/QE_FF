
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine emptystates_cg_sub( c0_emp, cm_emp, bec_emp, f_emp, fsic_emp, n_empx,&
                          n_emps, ispin_emp, iupdwn_emp, nupdwn_emp, phi_emp, lambda_emp, &
                          maxiter_emp, wxd_emp, vsic_emp, sizvsic_emp, pink_emp, nnrx, rhovan_emp, &
                          deeq_sic_emp, nudx_emp, eodd_emp, etot_emp, &
                          filledstates_potential, nfi, tfirst, tlast, eigr, c0, bec, irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ema0bg, desc_emp, nlam_emp)
      !
      use kinds,          only: dp
      !
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt,&
                                nbspx, nbsp, ispin
      use gvect,          only: ngm, gstart
      use gvecs,          only: ngms
      use smallbox_gvec,  only: ngb
      use gvecw,          only: ngw
      use ions_base,      only: na, nat, nax, nsp
      use cell_base,      only: omega, alat, tpiba2
      use local_pseudo,   only: vps, rhops
      use io_global,                only : stdout, ionode, ionode_id
      use mp_global,                only : intra_bgrp_comm
      use constants,                only : pi, au_gpa
      use uspp,                     only : nhsa=> nkb, nhsavb=> nkbus, &
                                           betae => vkb, rhovan => becsum, &
                                           deeq, qq, nlcc_any
      use uspp_param,               only : nh, nvb, ish, nhm
      use cg_module,                only : ene_ok, niter_cg_restart, &
                                           conv_thr, passop, enever,itercg,c0old
      use mp,                       only : mp_sum 
      use cp_electronic_mass,       only : emass_cutoff
      use orthogonalize_base,       only : calphi_bgrp
      use cp_interfaces,            only : rhoofr_generalized, dforce, prefor
      use cp_interfaces,            only : distribute_lambda, nlsm1
      use descriptors,              only : la_descriptor
      use fft_base,                 only : dffts, dfftp
      use nksic,                    only : do_orbdep, do_innerloop, do_innerloop_cg, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, innerloop_init_n, innerloop_cg_ratio, &
                                           vsicpsi, vsic, wtot, fsic, fion_sic, deeq_sic, f_cutoff, &
                                           pink, do_wxd, sizwtot, do_bare_eigs, innerloop_until, &
                                           valpsi, odd_alpha, eodd
      use emptystates_electrons_module, only : wfc_spreads_emp, wfc_centers_emp
      use electrons_module,             only : icompute_spread
      !
      implicit none
      !
      integer     :: nfi, nlam_emp
      logical     :: tfirst , tlast
      integer     :: sizvsic_emp
      integer     :: n_emps, n_empx, iupdwn_emp(nspin), nupdwn_emp(nspin), maxiter_emp, nnrx, &
                     nudx_emp, ispin_emp(n_empx) 
      complex(dp) :: eigr(ngw,nat)
      complex(dp) :: c0(ngw,nbspx)
      real(dp)    :: bec(nhsa,nbspx) 
      real(dp)    :: bec_emp(nhsa, n_emps)
      real(dp)    :: dbec_emp(nhsa,nbspx,3,3)
      real(dp)    :: lambda_emp(nlam_emp,nlam_emp,nspin)
      integer     :: irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp)    :: rhor(dfftp%nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp)    :: rhos(dffts%nnr,nspin)
      real(dp)    :: rhoc(dfftp%nnr)
      real(dp)    :: ema0bg(ngw)
      real(dp)    :: f_emp(n_empx), fsic_emp(n_empx), wxd_emp(sizvsic_emp,2), vsic_emp(sizvsic_emp, n_empx), &
                     pink_emp(n_empx), rhovan_emp(nhm*(nhm+1)/2, nat, nspin), &
                     deeq_sic_emp(nhm,nhm,nat,n_empx), eodd_emp, etot_emp, & 
                     filledstates_potential(dffts%nnr,nspin)
      complex(dp) :: c0_emp(ngw, n_empx), cm_emp(ngw, n_empx), phi_emp(ngw, n_empx)
      TYPE(la_descriptor)::  desc_emp(nspin) 
      !
      ! local variables
      ! 
      integer  :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk, ip, isp
      integer  :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot , comm
      real(dp) :: enb, enbi, x
      real(dp) :: entmp, sta
      real(dp) :: gamma 
      complex(dp), allocatable :: c2(:)
      complex(dp), allocatable :: c3(:)
      complex(dp), allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:)
      real(dp), allocatable :: s_minus1(:,:)!factors for inverting US S matrix
      real(dp), allocatable :: k_minus1(:,:)!factors for inverting US preconditioning matrix
      real(DP), allocatable :: lambda_repl(:,:) ! replicated copy of lambda
      real(DP), allocatable :: lambda_dist(:,:) ! replicated copy of lambda
      !
      real(dp)    :: sca, dumm(1)
      logical     :: newscheme, firstiter
      integer     :: maxiter3
      !
      real(DP), allocatable :: bec0(:,:), becm(:,:) 
      real(DP), allocatable :: fmat_(:,:)!average kinetic energy for preconditioning
      real(DP) :: esse, essenew !factors in c.g.
      logical     :: ltresh!flag for convergence on energy
      real(DP)    :: passo!step to minimum
      real(DP)    :: etotnew, etotold!energies
      real(DP)    :: spasso!sign of small step
      logical     :: restartcg!if .true. restart again the CG algorithm, performing a SD step
      integer     :: numok!counter on converged iterations
      integer     :: iter3
      real(DP)    :: passof,passov !step to minimum: effective, estimated
      real(DP)    :: ene0,ene1,dene0,enesti !energy terms for linear minimization along hi
      !
      real(DP),    allocatable :: faux(:) ! takes into account spin multiplicity
      complex(DP), allocatable :: hitmp(:,:)
      integer     :: ninner,nbnd1,nbnd2,itercgeff
      real(DP)    :: dtmp, temp
      real(dp)    :: tmppasso, ene_save(100), ene_save2(100)
      !
      logical     :: switch=.false., ortho_switch=.false., okvan, steepest=.false.
      complex(DP) :: phase
      integer     :: ierr, northo_flavor
      real(DP)    :: deltae, sic_coeff1, sic_coeff2 !coefficients which may change according to the flavour of SIC
      integer     :: iunit_manifold_overlap, iunit_spreads
      character(len=10) :: tcpu_cg_here
      real(DP):: ekin_emp, enl_emp, dekin_emp(6), denl_emp(3,3), epot_emp
      real(DP), allocatable :: rhor_emp(:,:), rhos_emp(:,:), rhoc_emp(:)
      real(DP), allocatable :: drhor_emp(:,:,:,:) 
      complex(DP), allocatable :: drhog_emp(:,:,:,:)
      complex(DP), allocatable :: rhog_emp(:,:)
      real(DP), allocatable    :: faux_emp(:)
      integer                  :: in_emp, issw
      COMPLEX(DP), PARAMETER   :: c_zero=CMPLX(0.d0,0.d0)
      !
      ! vars for numerical derivatives
      !
      REAl(DP):: etot_emp_tmp1, etot_emp_tmp2, etot_emp_tmp
      !
      call do_allocation_initialization()
      !
      ! Initializing clock for minimization
      !
      call start_clock('runcg_uspp_empty')
      ! 
      if ( tfirst .and. ionode ) &
         write(stdout,"(/,a,/)") 'PERFORMING CONJUGATE GRADIENT MINIMIZATION FOR EMPTY STATES'
      !         
      ! set tpa mass preconditioning
      !
      call emass_precond_tpa ( ema0bg, tpiba2, emass_cutoff )
      ! 
      call prefor(eigr, betae) 
      !
      call nlsm1( n_emps, 1, nsp, eigr, c0_emp, bec_emp )
      !
      ! recompute phi (the augmented wave function) from the new c0_empty
      !
      CALL calphi_bgrp( c0_emp, SIZE(c0_emp,1), bec_emp, nhsa, betae, phi_emp, n_emps)
      !
      ! calculates the factors for S and K inversion in US case -- they are important for preconditioning
      ! see paper by Hasnip and Pickard, Computer Physics Communications 174 (2006) 24–29
      !
      if (nvb.gt.0) then
         !
         allocate( s_minus1(nhsavb,nhsavb))
         allocate( k_minus1(nhsavb,nhsavb))
         call set_x_minus1(betae,s_minus1,dumm,.false.)
         call set_x_minus1(betae,k_minus1,ema0bg,.true.)
         ! 
      else
         !
         allocate( s_minus1(1,1))
         allocate( k_minus1(1,1)) 
         !
      endif
      !
      ! set index on number of converged iterations
      !
      numok = 0
      !
      allocate( hpsi(ngw, n_emps) )
      allocate( hpsi0(ngw, n_emps) )
      allocate( gi(ngw, n_emps), hi(ngw, n_emps) )
      allocate( hitmp(ngw, n_emps) )
      ! 
      hitmp(:,:) = CMPLX(0.d0,0.d0)
      gi(:,:)=CMPLX(0.d0,0.d0)
      hi(:,:)=CMPLX(0.d0,0.d0)
      !
      ene_ok=.false.
      !
      !=======================================================================
      !                 begin of the main loop
      !=======================================================================
      !
      OUTER_LOOP: &
      !
      do while ( itercg .lt. maxiter_emp .and. (.not.ltresh) )
        !
        call start_clock( "outer_loop" )
        ! 
        ENERGY_CHECK: &
        if (.not. ene_ok ) then
           !
           !call nlsm1( n_emps, 1, nsp, eigr, c0_emp, bec_emp) 
           ! 
           call rhoofr_generalized (n_empx, n_emps, ispin_emp, f_emp, nfi, c0_emp, irb, eigrb, &
                             bec_emp, dbec_emp, rhovan_emp, rhor_emp, drhor_emp, rhog_emp, drhog_emp, & 
                             rhos_emp, enl_emp, denl_emp, ekin_emp, dekin_emp)
           !
           etot_emp = enl_emp + ekin_emp                                 
           !
           call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
           !
           etot_emp = etot_emp + epot_emp
           ! 
           etotnew = etot_emp
           !
        else
           !
           etot_emp = enever
           !
           etotnew = etot_emp
           ene_ok = .false.
           !
        endif ENERGY_CHECK
        !
        call print_out_observables()
        !
        ! here we store the etot in ene0, to keep track of the energy of the initial point
        !
        call check_convergence_cg()        
        !
        call compute_hpsi()
        !
        ! HPSI IS ORTHOGONALIZED TO c0
        !
        ! comp. <beta|hpsi>
        !  
        !if (switch.or.(.not.do_orbdep).or.(do_orbdep.and.wo_odd_in_empty_run) ) then
        if (switch.or.(.not.do_orbdep)) then
           call pcdaga2_generalized(c0_emp, phi_emp, hpsi, n_emps, ispin_emp )
           !
        else
           !
           call pc2nc_nksic(c0_emp, hpsi, n_emps, ispin_emp)
           !
        endif
        !
        CALL nlsm1 ( n_emps, 1, nsp, eigr, hpsi, becm )
        !
        do iss=1,nspin
           !
           in_emp = iupdwn_emp( iss )
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true. , eigr, betae, becm, bec, nhsa, &
                            hpsi( :, in_emp: ), c0( :, issw: ), ngw, &
                            nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo
        ! 
        call nlsm1( n_emps, 1, nsp, eigr, hpsi, becm )  
        !
        ! TWO VECTORS INITIALIZED TO HPSI
        ! 
        hpsi0(:,:) = hpsi(:,:) 
        gi(:,:)    = hpsi(:,:)
        !
	! COMPUTES ULTRASOFT-PRECONDITIONED HPSI,
        ! non kinetic-preconditioned, 
        ! is the subsequent reorthogonalization necessary 
        ! in the norm conserving case???: giovanni
        call xminus1_generalized(n_emps, hpsi, betae, dumm, becm, s_minus1,.false.)
        !
        call orthogonalize(c0_emp, hpsi, becm, bec_emp)
        !
        call xminus1_generalized(n_emps, gi, betae, ema0bg, becm, k_minus1, .true.)
        !
        call orthogonalize(c0_emp, gi, becm, bec_emp)
        !    
        !  calculates gamma
        !
        gamma = 0.0_dp
        !
        do i=1, n_emps
           !
           do ig=1,ngw
              ! 
              gamma=gamma+2.d0*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
              ! 
           enddo
           !
           if (gstart==2) then 
              !
              gamma=gamma-DBLE(CONJG(gi(1,i))*hpsi(1,i))
              !
           endif
           !
        enddo
        !        
        call mp_sum( gamma, intra_bgrp_comm )
        !
        if (nvb.gt.0) then
           !
           do i=1, n_emps
              ! 
              do is=1,nvb
                 !
                 do iv=1,nh(is)
                    !
                    do jv=1,nh(is)
                       !
                       do ia=1,na(is)
                          !
                          inl=ish(is)+(iv-1)*na(is)+ia
                          !
                          jnl=ish(is)+(jv-1)*na(is)+ia
                          ! 
                          gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*bec0(jnl,i)
                          !
                       enddo
                       !
                    enddo
                    !
                 enddo
                 !
              enddo
              !
           enddo
           ! 
        endif
        !  
        ! case of first iteration
        ! 
        if (steepest) then
           !
           ! steepest descent
           !  
           gamma=0.d0
           !
        endif
        !
        if (itercg==1 .or. mod(itercg, niter_cg_restart)==0 .or. restartcg) then
           !
           restartcg=.false.
           !
           passof=passop
           !
           ! hi is the search direction
           !
           hi(:,:)=gi(:,:) 
           ! 
           esse=gamma
           ! 
        else
           !
           ! find direction hi for general case 
           ! calculates gamma for general case, not using Polak Ribiere
           !
           if (.not.steepest) then
              !
              essenew=gamma
              !
              gamma=gamma/esse
              ! 
              esse=essenew
              !
           else
              !    
              esse=0.d0
              !
              essenew=0.d0
              !
           endif
           !   
           hi(:,:) = gi(:,:) + gamma * hi(:,:)
           !   
        endif
        !
        ! note that hi is saved on gi, because we need it before projection on conduction states
        !     
        ! ... find minimum along direction hi:
        !
        ! project hi on conduction sub-space
        !
        call orthogonalize(c0_emp, hi, bec0, bec_emp)
        !    
        ! do quadratic minimization
        !             
        ! calculate derivative with respect to lambda along direction hi
        !
        dene0=0.d0
        ! 
        do i=1, n_emps
           !  
           do ig=1,ngw
              !
              dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
              !
           enddo
           !
           if (gstart==2) then
              ! 
              dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
              !   
           endif
           !
        enddo
        !
        ! We need the following because n for spin 2 is double that for spin 1!
        ! need to be check, Linh
        !
        dene0 = dene0 * 2.d0/nspin
        !
        call mp_sum( dene0, intra_bgrp_comm )
        !
        ! if the derivative is positive, search along opposite direction
        ! 
        if (dene0.gt.0.d0) then
           !
           spasso=-1.D0
           ! 
        else
           !
           spasso=1.d0
           !
        endif
        !
        ! calculates wave-functions on a point on direction hi
        !
        cm_emp(:,:) = c0_emp(:,:) + spasso * passof * hi(:,:)
        !
        if (gstart==2) cm_emp(1,:) = 0.5d0*(cm_emp(1,:) + CONJG(cm_emp(1,:)))
        ! 
        ! orthonormalize
        !
        call orthogonalize_wfc_only(cm_emp, becm)
        !
        call rhoofr_generalized (n_empx, n_emps, ispin_emp, f_emp, nfi, cm_emp, irb, eigrb, &
                             becm, dbec_emp, rhovan_emp, rhor_emp, drhor_emp, rhog_emp, drhog_emp, &
                             rhos_emp, enl_emp, denl_emp, ekin_emp, dekin_emp)
        !
        etot_emp =  enl_emp + ekin_emp 
        !
        call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp)
        !
        etot_emp =  etot_emp + epot_emp
        !
        ene1=etot_emp
        !    
        ! find the minimum
        !
        call minparabola(ene0, spasso*dene0, ene1, passof, passo, enesti)
        !
        write(stdout,"(6f20.12)") ene0,dene0,ene1,passo,gamma,esse
        ! 
        ! set new step
        !
        passov=passof
        passof=2.0*passo
        !      
        ! calculates wave-functions at minimum
        !
        cm_emp(:,:) = c0_emp(:,:) + spasso * passo * hi(:,:)
        !
        IF (gstart==2) cm_emp(1,:) = 0.5d0*(cm_emp(1,:)+CONJG(cm_emp(1,:)))
        !
        call orthogonalize_wfc_only(cm_emp, becm)
        !
        call rhoofr_generalized (n_empx, n_emps, ispin_emp, f_emp, nfi, cm_emp, irb, eigrb, &
                             becm, dbec_emp, rhovan_emp, rhor_emp, drhor_emp, rhog_emp, drhog_emp, &
                             rhos_emp, enl_emp, denl_emp, ekin_emp, dekin_emp)
        !
        etot_emp =  enl_emp + ekin_emp
        !
        call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp)
        !
        etot_emp =  etot_emp + epot_emp
        !
        enever=etot_emp
        !
        ! check with  what supposed
        !
        write(stdout,*) 'ene0, dene0, ene1, enesti,enever, passo, passov, passof'
        write(stdout,"(7f18.12)") ene0, dene0, ene1, enesti,enever, passo, passov, passof
        !
        ! if the energy has diminished with respect to ene0 and ene1 , everything ok
        !
        if (((enever.lt.ene0) .and. (enever.lt.ene1))) then
           !
           c0_emp(:,:) = cm_emp(:,:)
           !
           bec_emp(:,:)= becm(:,:)
           !
           ene_ok=.true.
           !
           ! if  ene1 << energy < ene0; go to  ene1
           !
        elseif((enever.ge.ene1) .and. (enever.lt.ene0)) then
           !
           if (ionode) then
              ! 
              write(stdout,"(2x,a,i5,f20.12)") 'cg_sub: missed minimum, case 1, iteration',itercg, passof
              ! 
           endif
           ! 
           c0_emp(:,:) = c0_emp(:,:) + spasso*passov*hi(:,:)
           ! 
           restartcg = .true.
           !
           call orthogonalize_wfc_only(c0_emp, bec_emp)
           !
           ene_ok=.false.
           !
           ! if  ene1 << ene0 <= energy; go to  ene1
           !
        elseif((enever.ge.ene0).and.(ene0.gt.ene1)) then
           !   
           if (ionode) then
              write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 2, iteration',itercg
           endif
           ! 
           c0_emp(:,:) = c0_emp(:,:) + spasso * passov * hi(:,:)
           !
           restartcg = .true.
           !
           call orthogonalize_wfc_only(c0_emp, bec_emp)
           !
           ene_ok=.false.
           !  
           ! if ene > ene0, en1 do a steepest descent step
           !
        elseif((enever.ge.ene0).and.(ene0.le.ene1)) then
           !
           if(ionode) then
             write(stdout,"(2x,a,i5)") 'cg_sub: missed minimum, case 3, iteration, doing steepest descent',itercg
           endif
           !
           iter3=0
           !  
           do while(enever.ge.ene0 .and. iter3.lt.maxiter3)
             ! 
             iter3 = iter3 + 1
             !  
             passov = passov*0.5d0
             !
             cm_emp(:,:) = c0_emp(:,:) + spasso*passov*hi(:,:)
             !
             ! change the searching direction
             !
             spasso=spasso*(-1.d0)
             !
             call orthogonalize_wfc_only(cm_emp, becm)
             !   
             itercgeff = itercgeff+1
             call rhoofr_generalized (n_empx, n_emps, ispin_emp, f_emp, nfi, cm_emp, irb, eigrb, &
                             becm, dbec_emp, rhovan_emp, rhor_emp, drhor_emp, rhog_emp, drhog_emp, &
                             rhos_emp, enl_emp, denl_emp, ekin_emp, dekin_emp)
             !
             etot_emp = enl_emp + ekin_emp
             !
             call v_times_rho(filledstates_potential, nspin, rhos_emp, epot_emp) 
             ! 
             etot_emp = etot_emp + epot_emp 
             !
             enever=etot_emp
             !
           enddo
           !
           if (ionode) write(stdout,"(2x,a,i5)") 'iter3 = ',iter3
           !
           if (iter3 == maxiter3 .and. enever.gt.ene0) then
              ! 
              write(stdout,"(2x,a)") 'missed minimum: iter3 = maxiter3'
              write(stdout,*) enever, ene0
              !
           elseif (enever.le.ene0) then
              !
              c0_emp(:,:)=cm_emp(:,:)
              ! 
              bec_emp(:,:)=becm(:,:)
              !
           endif
           !
           restartcg=.true.
           ene_ok=.false.
           !
           if (iter3 == maxiter3) then
              !
              passof=passop
              !
           endif
           !
        endif
        !  
        if (.not.ene_ok) call nlsm1 ( n_emps, 1, nsp, eigr, c0_emp, bec_emp )
        !
        ! calculates phi for pc_daga
        !
        call calphi_bgrp( c0_emp, SIZE(c0_emp,1), bec_emp, nhsa, betae, phi_emp, n_emps )
        ! 
        !=======================================================================
        !                 end of the outer loop
        !=======================================================================
        !
        itercg=itercg+1
        !
        itercgeff=itercgeff+1
        !
        call stop_clock( "outer_loop" )
        !
      enddo OUTER_LOOP
      !
      !=======================================================================
      !                 end of the main loop
      !=======================================================================
      !
      ! faux takes into account spin multiplicity.
      !
      faux(:) = f_emp(:) * DBLE( nspin ) / 2.0d0
      !
      do i = 1, n_emps, 2
         !
         call start_clock( 'dforce2' )
         !
         call dforce(i, bec_emp, betae, c0_emp, c2, c3, filledstates_potential, dffts%nnr, ispin_emp, faux, n_emps, nspin)
         !
         call start_clock( 'dforce2' )
         ! 
         do ig=1, ngw
            !
            gi(ig, i)=c2(ig)
            !
            if(i+1 <= n_emps) gi(ig,i+1)=c3(ig)
            ! 
         enddo
         !
         if (gstart==2) then
            ! 
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0, kind=DP)
            !
            if(i+1 <= n_emps) gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0, kind=DP)
            !    
         endif
         !
      enddo
      !
      allocate(lambda_repl(nudx_emp, nudx_emp))
      !
      hitmp(:,:) = c0_emp(:,:)
      !
      do is = 1, nspin
         !
         nss = nupdwn_emp(is)
         istart = iupdwn_emp(is)
         lambda_repl = 0.d0
         !
         do i = 1, nss
            ! 
            do j = i, nss
               !
               ii = i + istart - 1
               jj = j + istart - 1
               !
               do ig = 1, ngw
                  !
                  lambda_repl( i, j ) = lambda_repl( i, j ) - &
                     2.d0 * DBLE( CONJG( c0_emp( ig, ii ) ) * gi( ig, jj) )
                  !
               enddo
               !
               if ( gstart == 2 ) then
                  !  
                  lambda_repl( i, j ) = lambda_repl( i, j ) + &
                      DBLE( CONJG( c0_emp( 1, ii ) ) * gi( 1, jj ) )
                  !
               endif
               !
               lambda_repl( j, i ) = lambda_repl( i, j )
               !
            enddo
            !
         enddo
         !
         call mp_sum( lambda_repl, intra_bgrp_comm )
         call distribute_lambda( lambda_repl, lambda_emp(:,:, is),  desc_emp(is) )
         !
      enddo
      !
      deallocate( lambda_repl )
      !
      call do_deallocation()
      !
      return
      !
      contains
      !
     subroutine do_allocation_initialization()
         !  
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! INITIALIZATION PART (variables, allocation of arrays, minimization parameters)
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !!! Initialize some basic variables
         !
         okvan=(nvb>0)
         !
         deltae = 2.d0*conv_thr
         etotnew=0.d0
         etotold=0.d0
         !
         northo_flavor=1
         !
         IF(northo_flavor==2) THEN
            !
            sic_coeff1=1.d0
            sic_coeff2=1.d0
            !
         ELSE
            !
            sic_coeff1=0.5d0
            sic_coeff2=0.5d0
            !
         ENDIF
         !
         allocate (faux(n_empx))
         allocate (faux_emp(n_empx))
         !
         faux_emp=0.d0
         !
         allocate (c2(ngw),c3(ngw))
         !
         allocate(rhor_emp(dfftp%nnr,nspin), rhos_emp(dffts%nnr,nspin), rhog_emp(ngm,nspin))
         !
         allocate( drhog_emp( 1, 1, 1, 1 ) )
         allocate( drhor_emp( 1, 1, 1, 1 ) )
         !
         if (nlcc_any) allocate(rhoc_emp(dfftp%nnr))
         !  
         ! Allocation of twin-type variables
         !
         allocate (bec0(nhsa, n_emps) )
         allocate (becm(nhsa, n_emps) )
         !
         ! initializing variables for checking iterations and convergence
         !
         newscheme=.false.
         firstiter=.true.
         maxiter3=12
         !
         if(do_orbdep) maxiter3=10
         !
         ninner=0
         !
         ltresh    = .false.
         itercg    = 1
         etotold   = 1.d8
         restartcg = .true.
         passof = passop
         ene_ok = .false.
         !
         itercgeff = 1
         !
     end subroutine do_allocation_initialization
     
     subroutine do_deallocation()
        !
        deallocate(hpsi0,hpsi,gi,hi)
        deallocate(hitmp, STAT=ierr)
        !        
        deallocate(s_minus1, k_minus1)
        !
        call stop_clock('runcg_uspp_empty')
        !
        deallocate(bec0, becm)
        !
        deallocate(c2,c3)
        deallocate(faux)
        deallocate(faux_emp)
        deallocate(rhor_emp,rhos_emp,rhog_emp)
        deallocate(drhor_emp,drhog_emp)
        !
        if (nlcc_any) deallocate(rhoc_emp)
        !
     end subroutine do_deallocation
     ! 
     subroutine do_innerloop_subroutine()
        !
	!if (do_innerloop_empty .and. innerloop_until>=itercgeff) then
           !
        !   call start_clock( "inner_loop" )
        !   !
        !   eodd_emp= sum(pink_emp(:))
        !   etot_emp= etot_emp - eodd_emp
        !   etotnew = etotnew  - eodd_emp
        !   ninner  = 0
           !
        !   if (.not.do_innerloop_cg) then
              ! 
        !      write(stdout,*)  "WARNING, do_innerloop_cg should be .true."
              ! 
        !   else
              !
        !      call nksic_rot_emin_cg_general(itercg,innerloop_init_n,ninner,etot_emp,deltae*innerloop_cg_ratio,lgam, &
        !                             n_emps, n_empx, nudx_emp, iupdwn_emp, nupdwn_emp, ispin_emp, c0_emp, rhovan_emp, bec_emp, rhor, rhoc, &
        !                             vsic_emp, pink_emp, deeq_sic_emp, wtot, fsic_emp, sizwtot, .false.,  wfc_centers_emp, wfc_spreads_emp, .true.) 
              !
        !   endif
           !
        !   eodd_emp= sum(pink_emp(:)) 
        !   etot_emp= etot_emp + eodd_emp 
        !   etotnew = etotnew  + eodd_emp
           ! 
        !   call stop_clock( "inner_loop" )
           !
        !endif
        !
     endsubroutine do_innerloop_subroutine
     !   
     subroutine print_out_observables()
        ! 
        call print_clock('CP')
        !
        if ( ionode ) then
           !
           if (itercg>2) then
              write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etotemp (Ha) =",F22.14," delta_Eemp=",E22.14)')&
              itercg, itercgeff, etotnew, deltae
           else
              write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4,"   Etotemp (Ha) =",F22.14)')&
              itercg, itercgeff, etotnew
           endif
           !
           write(stdout,'(5x,  "Ekin (Ha) = ",F22.14 , " Enl (Ha) = ",F22.14, " Eloc (Ha) =" , F22.14)' )&
              ekin_emp, enl_emp, epot_emp 
           !
           !if ( do_orbdep .and. (.not. wo_odd_in_empty_run)  ) then 
           !   write(stdout,'(1x,  "Fake EODD (Ha) = ",F22.14) ') eodd_emp
           !endif
           !
        endif
        !
        if ( ionode .and. mod(itercg,10) == 0 ) write(stdout,"()" )
        !
        !if ( ionode .and. mod(itercg, iprint_spreads)==0) then
        if ( .false.) then
           !
           if(nspin==1) then
             !
             write( iunit_spreads, '(400f20.14)') wfc_spreads_emp(:,1,2)               
             !
           elseif(nspin==2) then
             !
             write( iunit_spreads, '(2(400f20.14)(3x))') wfc_spreads_emp(:,1,2), wfc_spreads_emp(:,2,2)
             !
           endif
           !
        endif
        !
     end subroutine print_out_observables
     !
     subroutine check_convergence_cg()
        !
        deltae=abs(etotnew-etotold)
        !
        if( deltae < conv_thr ) then
           numok=numok+1
        else 
           numok=0
        endif
        !
        if( numok >= 4 ) ltresh=.true.
        !
        if(ltresh.or.itercg==maxiter_emp-1) icompute_spread=.true.
        !
        etotold=etotnew
        ene0=etot_emp
        !
     end subroutine check_convergence_cg
     !
     subroutine compute_hpsi()
        ! 
        ! faux takes into account spin multiplicity.
        !
        faux(:)=0.d0
        faux(1:n_emps) = f_emp(1:n_emps) * DBLE( nspin ) / 2.0d0
        !
        do i=1, n_emps ,2
           ! 
           call dforce( i, bec_emp, betae, c0_emp, c2, c3, filledstates_potential, dffts%nnr, ispin_emp, faux, n_emps, nspin) 
           !
           ! ODD terms
           !
           !if ( do_orbdep .and. (.not. wo_odd_in_empty_run) ) then
              !   
           !   CALL nksic_eforce( i, n_emps, n_empx, vsic_emp, deeq_sic_emp, bec_emp, ngw, &
           !                            c0_emp(:,i), c0_emp(:,i+1), vsicpsi, lgam )
           !   !
           !   c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
           !   !
           !   if( i+1 <= n_emps )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
              !
           !endif
           !
           hpsi(1:ngw, i)=c2(1:ngw)
           !
           if (i+1 <= n_emps ) then
              hpsi(1:ngw, i+1)=c3(1:ngw)
           endif
           !
           if (gstart == 2) then
              ! 
              hpsi(1, i)=CMPLX(DBLE(hpsi(1, i)), 0.d0)
              !
              if (i+1 <= n_emps) then
                 !
                 hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
                 ! 
              endif
              !
           endif
           !
        enddo
        !  
     end subroutine compute_hpsi

     subroutine orthogonalize(wfc0, wfc, becwfc, bec0)
        ! 
        real(DP) :: becwfc(:,:), bec0(:,:)
        complex(DP) :: wfc(:,:), wfc0(:,:)
        !
        if (switch.or.(.not.do_orbdep)) then !.or.(do_orbdep .and.wo_odd_in_empty_run)) then
           !
           call pc2_generalized(wfc0, bec0, wfc, becwfc, n_emps, &
                        nupdwn_emp, iupdwn_emp, ispin_emp)
           !
        else
           !
           if (.not.okvan) then
              !
              call pc2nc_nksic(wfc0, wfc, n_empx, ispin_emp )
              !
           else
              ! 
              call pc2_generalized(wfc0, bec0, wfc, becwfc, n_emps, &
                        nupdwn_emp, iupdwn_emp, ispin_emp)
              !
           endif
           !
        endif
        ! 
        CALL nlsm1 ( n_emps, 1, nsp, eigr, wfc, becwfc )
        !
        do iss=1,nspin
           !
           in_emp = iupdwn_emp( iss )
           issw   = iupdwn( iss )
           !
           CALL gram_empty(.true., eigr, betae, becwfc, bec, nhsa, &
                            wfc( :, in_emp: ), c0( :, issw: ), &
                            ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)
           !
        enddo
        !  
     end subroutine orthogonalize
     !
     subroutine orthogonalize_wfc_only(wfc,becwfc)
        !
        real(DP) :: becwfc(:,:)
        complex(DP) :: wfc(:,:)
        !
        call nlsm1 ( n_emps, 1, nsp, eigr, wfc, becwfc )
        !
        do iss=1,nspin
           !
           issw   = iupdwn(iss)
           in_emp = iupdwn_emp(iss)
           !
           CALL gram_empty(.false., eigr, betae, becwfc, bec, nhsa, &
                            wfc( :, in_emp: ), c0( :, issw: ), &
                            ngw, nupdwn_emp(iss), nupdwn(iss), in_emp, issw)           
           !
        enddo
        ! 
        call nlsm1 ( n_emps, 1, nsp, eigr, wfc, becwfc ) 
        !
     end subroutine orthogonalize_wfc_only
     !
     subroutine v_times_rho(v, nspin, rhos_emp, epot_emp)
        !
        use kinds, only: DP
        use mp,                   only : mp_sum
        use mp_global,            only : intra_bgrp_comm 
        use cell_base,            only : omega
        use fft_base,             only : dffts 
        !
        implicit none
        ! 
        integer, intent(in)   :: nspin
        real(DP), intent(in)  :: v(dffts%nnr,nspin), rhos_emp(dffts%nnr,nspin)
        real(DP), intent(out) :: epot_emp
        !
        ! local vars
        !   
        integer  :: i
        real(DP) :: etemp, fact, rhosum(2)
        ! 
        etemp=0.d0
        rhosum=0.d0
        fact=omega/DBLE(dffts%nr1*dffts%nr2*dffts%nr3)
        ! 
        do i=1,nspin
           !
           etemp = etemp + sum(v(:,i) * rhos_emp(:,i))
           ! 
           rhosum(i) =  sum(rhos_emp(:,i))
           !
        enddo
        !  
        call mp_sum(etemp, intra_bgrp_comm )
        !  
        call mp_sum(rhosum, intra_bgrp_comm)
        !
        epot_emp = etemp*fact
        !
        rhosum = rhosum*fact
        !
        return
        !
     end subroutine v_times_rho
     !                     
END SUBROUTINE emptystates_cg_sub
