!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! THIS ROUTINE HAS BEEN ADOPTED FROM CG_SUB.F90 TO KOOPMANS CODE
!
!=======================================================================
!
subroutine kcp_runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, fion, ema0bg, becdr, &
      lambdap, lambda, nlam, vpot, c0, cm, phi, dbec,l_cprestart  )
      !
      use kinds,                    only : dp
      use control_flags,            only : tpre, iverbosity, tfor, tprnfor
      ! 
      !---ensemble-DFT
      !
      use energies,                 only : eht, epseu, exc, etot, eself, enl, ekin,&
                                           atot, entropy, egrand
      use electrons_base,           only : f, nspin, nel, iupdwn, nupdwn, nudx, nelt,&
                                           nbspx, nbsp, ispin
      use ensemble_dft,             only : tens, ef, z0t, c0diag,  &
                                           becdiag, fmat0, e0,  id_matrix_init
      !
      !---
      !
      use gvect,                    only : ngm, gstart
      use gvecs,                    only : ngms
      use smallbox_gvec,            only : ngb
      use gvecw,                    only : ngw
      use ions_base,                only : na, nat, nax, nsp, rcmax
      use cell_base,                only : omega, alat, tpiba2
      use local_pseudo,             only : vps, rhops
      use io_global,                only : stdout, ionode, ionode_id
      use mp_global,                only : intra_bgrp_comm, np_ortho, me_ortho, ortho_comm
      use dener
      use constants,                only : pi, au_gpa
      use io_files,                 only : tmp_dir, prefix
      use uspp,                     only : nhsa=> nkb, nhsavb=> nkbus, &
                                           betae => vkb, rhovan => becsum, &
                                           deeq, qq, nlcc_any
      use uspp_param,               only : nh, nvb, ish
      use cg_module,                only : ene_ok,  maxiter,niter_cg_restart, &
                                           conv_thr, passop, enever, itercg,c0old
      use ions_positions,           only : tau0
      use efield_module,            only : tefield, evalue, ctable, qmat, detq, ipolp, &
                                           berry_energy, ctabin, gqq, gqqm, df, pberryel, &
                                           tefield2, evalue2, ctable2, qmat2, detq2, ipolp2, &
                                           berry_energy2, ctabin2, gqq2, gqqm2, pberryel2
      use mp,                       only : mp_sum, mp_bcast
      use cp_electronic_mass,       only : emass_cutoff
      use orthogonalize_base,       only : calphi_bgrp
      use cp_interfaces,            only : rhoofr, dforce, compute_stress, vofrho, nlfl_bgrp, prefor
      use cp_interfaces,            only : nlsm2_bgrp, calbec, caldbec_bgrp, nlfq_bgrp
      use cp_interfaces,            only : collect_lambda, distribute_lambda
      use cp_main_variables,        only : descla, drhor, drhog
      use descriptors,              only : la_descriptor, ldim_cyclic
      use mp_global,                only : me_image, my_image_id, nbgrp
      use fft_base,                 only : dffts, dfftp
      !
      ! --- Koopmans 
      ! 
      use nksic,                    only : do_orbdep, do_innerloop, do_innerloop_cg, innerloop_cg_nsd, &
                                           innerloop_cg_nreset, innerloop_init_n, innerloop_cg_ratio, &
                                           vsicpsi, vsic, wtot, fsic, fion_sic, deeq_sic, f_cutoff, &
                                           pink, do_wxd, sizwtot, do_bare_eigs, innerloop_until, &
                                           valpsi, odd_alpha, eodd
      use kcp_electrons_module,     only : wfc_spreads, wfc_centers, icompute_spread
      !
      implicit none
      !
      logical, intent(in) :: l_cprestart !if true prepares a CG->CP restart
      character(len=80)   :: uname
      character(len=6), external :: int_to_char
      integer :: nfi, nlam
      logical :: tfirst , tlast
      complex(dp) :: eigr(ngw,nat)
      real(dp) :: bec(nhsa,nbspx)
      real(dp) :: becdr(nhsa,nbspx,3)
      integer irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp) :: rhor(dfftp%nnr,nspin)
      real(dp) :: vpot(dfftp%nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp) :: rhos(dffts%nnr,nspin)
      real(dp) :: rhoc(dfftp%nnr)
      complex(dp) :: ei1(-dfftp%nr1:dfftp%nr1,nat)
      complex(dp) :: ei2(-dfftp%nr2:dfftp%nr2,nat)
      complex(dp) :: ei3(-dfftp%nr3:dfftp%nr3,nat)
      complex(dp) :: sfac( ngms, nsp )
      real(dp) :: fion(3,nat)
      real(dp) :: ema0bg(ngw)
      real(dp) :: lambdap(nlam,nlam,nspin)
      real(dp) :: lambda(nlam,nlam,nspin)
      complex(dp) :: c0( ngw, nbspx )
      complex(dp) :: cm( ngw, nbspx )
      complex(dp) :: phi( ngw, nbspx )
      real(dp)    :: dbec(nhsa,nbspx,3,3)
      !
      ! internal variables
      !
      integer :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk, ip, nrlx
      integer :: inl, jnl, niter, istart, nss, nrl, me_rot, np_rot , comm
      real(dp):: enb, enbi, x
      complex(dp), allocatable :: c2(:)
      complex(dp), allocatable :: c3(:)
      real(dp) :: gamma, entmp, sta
      complex(dp),allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:)
      real(DP), allocatable :: s_minus1(:,:)!factors for inverting US S matrix
      real(DP), allocatable :: k_minus1(:,:)!factors for inverting US preconditioning matrix 
      real(DP), allocatable :: lambda_repl(:,:) ! replicated copy of lambda
      real(DP), allocatable :: lambda_dist(:,:) ! replicated copy of lambda
      real(dp) :: sca, dumm(1)
      logical  :: newscheme, firstiter
      integer  :: maxiter3
      !
      real(DP), allocatable :: bec0(:,:), becm(:,:), becdrdiag(:,:,:)
      real(DP), allocatable :: ave_ene(:)!average kinetic energy for preconditioning
      real(DP), allocatable :: fmat_(:,:)!average kinetic energy for preconditioning
      real(DP), allocatable :: faux(:)!takes into account spin multiplicity 
      logical :: pre_state!if .true. does preconditioning state by state
      real(DP):: esse,essenew !factors in c.g.
      logical :: ltresh!flag for convergence on energy
      real(DP):: passo!step to minimum
      real(DP):: etotnew,etotold!energies
      real(DP):: eoddnew,eoddold !odd energies | KOOPMANS FUNCTIONALS
      real(DP):: spasso!sign of small step
      logical :: restartcg!if .true. restart again the CG algorithm, performing a SD step
      integer :: numok!counter on converged iterations
      integer :: iter3
      real(DP):: passof,passov !step to minimum: effective, estimated
      real(DP):: ene0,ene1,dene0,enesti !energy terms for linear minimization along hi
      integer :: i_max
      real(DP):: max_sca
      real(DP):: deltae
      logical :: switch=.false.
      integer :: ninner, itercgeff ! KOOPMANS FUNCTIONALS
      ! 
      nrlx = maxval(descla(:)%nrlx)
      !
      allocate(bec0(nhsa,nbspx),becm(nhsa,nbspx), becdrdiag(nhsa,nbspx,3))
      allocate(ave_ene(nbspx))
      allocate(c2(ngw),c3(ngw))
      allocate(faux(nbspx))
      ! 
      call start_clock('runcg_uspp')
      !
      if(nbgrp > 1) &
        call errore(' runcg_uspp ', ' parallelization over bands not yet implemented ', 1 )
      !
      ! default parameter
      !
      newscheme=.false.
      firstiter=.true.
      pre_state=.false.!normally is disabled
      maxiter3=250
      !
      if(ionode) then
         uname = TRIM( tmp_dir ) // trim(prefix) // '.' &
                 // trim(int_to_char( my_image_id )) // '_' // trim(int_to_char( me_image))
         !open(37,file='convergence.dat',status='unknown')!for debug and tuning purposes
         open(37,file=uname,status='unknown')!for debug and tuning purposes
      endif
      ! 
      if( tfirst .and. ionode) &
         write(stdout,*) 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'
      ! 
      !set tpa preconditioning
      ! 
      call  emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
      !      
      call prefor(eigr,betae) 
      !
      ltresh    = .false.
      itercg    = 1
      itercgeff = 1
      etotold   = 1.d8
      restartcg = .true.
      passof    = passop
      ene_ok    = .false.
      !
      ! orthonormalize c0
      !
      call orthogonalize_wfc_only(c0, bec)
      !call calbec(1,nsp,eigr,c0,bec)
      !call gram_bgrp( betae, bec, nhsa, c0, ngw )
      ! 
      ! calculates phi for pcdaga
      !
      call calphi_bgrp( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp )
      !
      ! calculates the factors for S and K inversion in US case
      ! 
      if (nvb.gt.0) then
         !
         allocate( s_minus1(nhsavb,nhsavb))
         allocate( k_minus1(nhsavb,nhsavb))
         call  set_x_minus1(betae,s_minus1,dumm,.false.)
         call  set_x_minus1(betae,k_minus1,ema0bg,.true.)
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
      ! initialize z0t
      !
      call id_matrix_init( descla, nspin )
      !
      allocate(hpsi(ngw,nbspx),hpsi0(ngw,nbspx),gi(ngw,nbspx),hi(ngw,nbspx))
      !
      do while ( itercg .lt. maxiter .and. (.not.ltresh) )
         !
         ENERGY_CHECK: if(.not. ene_ok) then
            !
            call calbec(1,nsp,eigr,c0,bec)
            !
            if (.not.tens) then
               !
               call rhoofr(nfi,c0(:,:),irb,eigrb,bec,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
               !
            else
               !
               if (newscheme.or.firstiter) then 
                  !
                  call inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                                        rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                                        sfac,c0,bec,dbec,firstiter,vpot)
                  firstiter=.false.
                  !
               endif
               !
               ! calculation of the rotated quantities
               !
               call rotate( nrlx, z0t, c0(:,:), bec, c0diag, becdiag )
               ! 
               ! calculation of rho corresponding to the rotated wavefunctions
               ! 
               call rhoofr(nfi,c0diag,irb,eigrb,becdiag,dbec, &
                           rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
               ! 
            endif
            ! 
            ! when cycle is restarted go to diagonal representation
            !
            if (mod(itercg,niter_cg_restart)==1 .and. itercg >=2) then
               ! 
               call rotate( nrlx, z0t, c0(:,:), bec, c0diag, becdiag )
               c0(:,:)=c0diag(:,:)
               bec(:,:)=becdiag(:,:)
               call id_matrix_init( descla, nspin )
               !
            endif
            ! 
            ! calculates the potential
            !
            ! put core charge (if present) in rhoc(r)
            !
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
            !
            !---ensemble-DFT
            !
            vpot = rhor
            !
            call vofrho(nfi,vpot,drhor,rhog,drhog,rhos,rhoc,tfirst,tlast,&
                        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            !
            if (.not.tens) then
               etotnew=etot
            else
               etotnew=etot+entropy
            endif
            !
            if (tefield) then!just in this case calculates elfield stuff at zeo field-->to be bettered
               ! 
               call berry_energy( enb, enbi, bec, c0(:,:), fion )
               etot=etot+enb+enbi
               !
            endif
            !
            if (tefield2) then!just in this case calculates elfield stuff at zeo field-->to be bettered
               !
               call berry_energy2( enb, enbi, bec, c0(:,:), fion )
               etot=etot+enb+enbi
               !
            endif
            !
            if (do_orbdep) then
               !
               fsic(:) = f(:)
               !
               call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                     ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, & 
                                     vsic, do_wxd, pink, nudx, wfc_centers, &
                                     wfc_spreads, icompute_spread, .false. )
               ! 
               eodd=sum(pink(1:nbsp))
               ! 
               etot=etot + eodd
               !
            endif
            !
         else
            !
            etot=enever
            if (.not.tens) then 
               etotnew=etot
            else
               etotnew=etot+entropy
            endif
            ene_ok=.false.
            ! 
         endif ENERGY_CHECK
         !
         ! KOOPMANS FUNCTIONAL
         !  
         if ( do_orbdep ) then
            !
            call do_innerloop_subroutine()
            ! 
         endif
         !
         call print_out_observables()
         !
! Linh put all below to rotuine check_convergence_cg()
         !
         deltae=abs(etotnew-etotold) 
         !
         if (do_orbdep) deltae=deltae+abs(eoddnew-eoddold)
         !   
         if (deltae.lt.conv_thr) then
            !  
            numok=numok+1
            !
         else
            ! 
            numok=0
            !
         endif
         !
         if(numok.ge.4) then
           ltresh=.true.
         endif
         !
         etotold=etotnew
         eoddold=eoddnew
         ene0=etot
         if(tens .and. newscheme) then
            ene0=ene0+entropy
         endif
         !
!Linh end check_convergence_cg()
         !
         ! updated
         !
         call newd(vpot,irb,eigrb,rhovan,fion)
         !  
         call prefor(eigr,betae)!ATTENZIONE
         !
         call compute_hpsi() 
         !
         if(pre_state) call ave_kin(c0,SIZE(c0,1),nbsp,ave_ene)
         !
         if(switch.or.(.not.do_orbdep)) then
           !  
           call pcdaga2(c0,phi,hpsi)
           !
         else 
           ! 
           call pc2nc_nksic(c0,hpsi, nbsp, ispin)
           ! 
         endif
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
         ! 
         call calbec(1,nsp,eigr,hpsi,becm)
         call xminus1(hpsi,betae,dumm,becm,s_minus1,.false.)
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !look if the following orthogonalization are really needed
         call orthogonalize(c0,hpsi,becm,bec)
         !call calbec(1,nsp,eigr,hpsi,becm)
         !call pc2(c0,bec,hpsi,becm)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !
         !  COMPUTES ULTRASOFT+KINETIC-preconditioned GI
         !  
         if (.not.pre_state) then
            !
            call xminus1(gi,betae,ema0bg,becm,k_minus1,.true.)
            !
         else
            !
            call xminus1_state(gi,betae,ema0bg,becm,k_minus1,.true.,ave_ene)
            !
         endif
         ! 
         call orthogonalize(c0,gi,becm,bec)
         !call calbec(1,nsp,eigr,gi,becm)
         !call pc2(c0,bec,gi,becm)
         ! 
         if (tens) call calcmt( nrlx, f, z0t, fmat0 )
         !  
         call calbec(1,nsp,eigr,hpsi,bec0) 
         ! 
         !calculates gamma
         !
         gamma=0.d0
         ! 
         if (.not.tens) then
            ! 
            do i=1,nbsp
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
               do i=1,nbsp
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
                              jnl=ish(is)+(jv-1)*na(is)+ia
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
         else
            !
            do iss=1,nspin
               !
               nss=nupdwn(iss)
               istart=iupdwn(iss)
               me_rot = descla( iss )%mype
               np_rot = descla( iss )%npc * descla( iss )%npr
               ! 
               allocate( fmat_ ( nrlx, nudx ) )
               !
               do ip = 1, np_rot
                  ! 
                  if ( me_rot == ( ip - 1 ) ) then
                     fmat_ = fmat0(:,:,iss)
                  endif
                  ! 
                  nrl = ldim_cyclic( nss, np_rot, ip - 1 )
                  call mp_bcast( fmat_ , ip - 1 , intra_bgrp_comm )
                  !
                  do i=1,nss
                     !
                     jj = ip
                     !
                     do j=1,nrl
                        !
                        do ig=1,ngw
                           ! 
                           gamma=gamma+2.d0*DBLE(CONJG(gi(ig,i+istart-1))*hpsi(ig,jj+istart-1))*fmat_(j,i)
                           ! 
                        enddo
                        !
                        if (gstart==2) then
                           !
                           gamma=gamma-DBLE(CONJG(gi(1,i+istart-1))*hpsi(1,jj+istart-1))*fmat_(j,i)
                           !
                        endif
                        !
                        jj = jj + np_rot
                        ! 
                     enddo
                     !
                  enddo
                  !
               enddo
               !
               deallocate( fmat_ )
               !
            enddo
            !
            if (nvb.gt.0) then
               ! 
               do iss=1,nspin
                  !
                  nss=nupdwn(iss)
                  istart=iupdwn(iss)
                  me_rot = descla( iss )%mype
                  np_rot = descla( iss )%npc * descla( iss )%npr
                  !
                  allocate( fmat_ ( nrlx, nudx ) )
                  !
                  do ip = 1, np_rot
                     ! 
                     if ( me_rot == ( ip - 1 ) ) then
                        fmat_ = fmat0(:,:,iss)
                     endif
                     !
                     nrl = ldim_cyclic( nss, np_rot, ip - 1 )
                     call mp_bcast( fmat_ , ip - 1 , intra_bgrp_comm )
                     !
                     do i=1,nss
                        !
                        jj = ip 
                        !
                        do j=1,nrl
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
                                       jnl=ish(is)+(jv-1)*na(is)+ia
                                       gamma=gamma+ qq(iv,jv,is)*becm(inl,i+istart-1)*bec0(jnl,jj+istart-1)*fmat_(j,i)
                                    enddo
                                    !
                                 enddo
                                 !
                              enddo
                              ! 
                           enddo
                           !
                           jj = jj + np_rot
                           !
                        enddo
                        !
                     enddo
                     !
                  enddo
                  ! 
                  deallocate( fmat_ )
                  ! 
               enddo
               !
            endif
            !
            call mp_sum( gamma, intra_bgrp_comm )
            ! 
         endif
         !
         ! case of first iteration
         ! 
         if (itercg==1.or.(mod(itercg,niter_cg_restart).eq.1).or.restartcg) then
            !
            restartcg=.false.
            passof=passop
            hi=gi !hi is the search direction
            esse=gamma
            !
         else
            !
            !find direction hi for general case 
            !calculates gamma for general case, not using Polak Ribiere
            ! 
            essenew=gamma
            gamma=gamma/esse
            esse=essenew
            !
            hi=gi+gamma*hi
            !
         endif
         ! 
         !note that hi, is saved on gi, because we need it before projection on conduction states
         !
         !find minimum along direction hi:
         !
         !project hi on conduction sub-space
         !
         call orthogonalize(c0,hi,bec0,bec)
         !call calbec(1,nsp,eigr,hi,bec0)
         !call pc2(c0,bec,hi,bec0)
         !
         !do quadratic minimization
         !             
         !calculate derivative with respect to  lambda along direction hi
         ! 
         dene0=0.0
         if (.not.tens) then
            !  
            do i=1,nbsp
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
         else
            ! 
            !in the ensemble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
            !calculation of the kinetic energy x=xmin      
            call calcmt( nrlx, f, z0t, fmat0 )
            do iss = 1, nspin
               !
               nss    = nupdwn(iss)
               istart = iupdwn(iss)
               me_rot = descla( iss )%mype
               np_rot = descla( iss )%npc * descla( iss )%npr
               ! 
               allocate( fmat_ ( nrlx, nudx ) )
               !
               do ip = 1, np_rot
                  ! 
                  if ( me_rot == ( ip - 1 ) ) then
                     !
                     fmat_ = fmat0(:,:,iss)
                     ! 
                  endif
                  !
                  nrl = ldim_cyclic( nss, np_rot, ip - 1 )
                  !
                  call mp_bcast( fmat_ , ip - 1 , intra_bgrp_comm )
                  !
                  do i=1,nss
                     !
                     jj = ip
                     !
                     do j=1,nrl
                        !
                        do ig=1,ngw
                           ! 
                           dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i+istart-1))*hpsi0(ig,jj+istart-1))*fmat_(j,i)
                           dene0=dene0-2.d0*DBLE(CONJG(hpsi0(ig,i+istart-1))*hi(ig,jj+istart-1))*fmat_(j,i)
                           !
                        enddo
                        !  
                        if (gstart==2) then
                           !  
                           dene0=dene0+DBLE(CONJG(hi(1,i+istart-1))*hpsi0(1,jj+istart-1))*fmat_(j,i)
                           dene0=dene0+DBLE(CONJG(hpsi0(1,i+istart-1))*hi(1,jj+istart-1))*fmat_(j,i)
                           !
                        endif
                        !
                        jj = jj + np_rot
                        !    
                     enddo
                     !
                  enddo
                  !
               enddo
               !
               deallocate( fmat_ )
               !
            enddo
            !
         endif
         !
         call mp_sum( dene0, intra_bgrp_comm )
         !
         !if the derivative is positive, search along opposite direction
         !
         if (dene0.gt.0.d0) then
            spasso=-1.D0
         else
            spasso=1.d0
         endif
         !
         !calculates wave-functions on a point on direction hi
         !
         cm=c0+spasso*passof*hi
         !
         !orthonormalize
         call orthogonalize_wfc_only(cm, becm)
         !
   !      call calbec(1,nsp,eigr,cm,becm)
    !     call gram_bgrp( betae, becm, nhsa, cm, ngw )
         ! 
         !calculate energy
         ! 
         if (.not.tens) then
            ! 
            call rhoofr(nfi,cm(:,:),irb,eigrb,becm,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
            !call rhoofr_generalized( nbspx, nbsp, ispin, f ,nfi,cm(:,:),irb,eigrb,becm,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)

            !
         else
            !
            if (newscheme) then 
               !
               call inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                        rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,dbec,.false., vpot  )  
               !
            endif
            !
            ! calculation of the rotated quantities
            ! 
            call rotate( nrlx, z0t, cm(:,:), becm, c0diag, becdiag )
            !
            ! calculation of rho corresponding to the rotated wavefunctions
            !
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag,dbec,rhovan,rhor, &
                        drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
            !
         endif
         !
         ! calculate potential
         !
         ! put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(irb,eigrb,rhoc)
         !
         vpot = rhor
         !
         call vofrho(nfi,vpot,drhor,rhog,drhog,rhos,rhoc,tfirst,tlast, &
                     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
         !
         if ( tefield ) then!to be bettered
            !  
            call berry_energy( enb, enbi, becm, cm(:,:), fion )
            !
            etot=etot+enb+enbi
            ! 
         endif
         ! 
         if ( tefield2 ) then!to be bettered
            ! 
            call berry_energy2( enb, enbi, becm, cm(:,:), fion )
            !
            etot=etot+enb+enbi
            !
         endif
         !
         if (do_orbdep) then
            !
            !
            call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic, & 
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, sizwtot, &
                                  vsic, do_wxd, pink, nudx, &
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            !
            eodd=sum(pink(1:nbsp))
            !
            etot = etot + eodd
            !
         endif
         ! 
         ene1=etot
         !
         if (tens.and.newscheme) then
            !
            ene1=ene1+entropy
            !
         endif
         !              
         !find the minimum
         ! 
         call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)
         ! 
         if( iverbosity > 0 ) write(6,*) ene0,dene0,ene1,passo, gamma, esse
         ! 
         !set new step
         !
         passov=passof
         passof=2.d0*passo
         !      
         !calculates wave-functions at minimum
         ! 
         cm=c0+spasso*passo*hi
         !
         if (gstart==2) then
            ! 
            cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
            !
         endif
         ! 
         call orthogonalize_wfc_only(cm, becm)
   !      call calbec(1,nsp,eigr,cm,becm)
    !     call gram_bgrp( betae, becm, nhsa, cm, ngw )
         ! 
         !test on energy: check the energy has really diminished
         !
         !call calbec(1,nsp,eigr,cm,becm)
         if (.not.tens) then
            !  
            call rhoofr(nfi,cm(:,:),irb,eigrb,becm,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
            !
         else
            !
            if (newscheme)  then
               !
               call inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,dbec,.false., vpot  )
               !
            endif
            !
            ! calculation of the rotated quantities
            !
            call rotate( nrlx, z0t, cm(:,:), becm, c0diag, becdiag )
            !
            ! calculation of rho corresponding to the rotated wavefunctions
            ! 
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag,dbec,rhovan,rhor, &
                     drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
            !
         endif
         !
         ! calculates the potential
         !
         ! put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(irb,eigrb,rhoc)
         !
         vpot = rhor
         !
         call vofrho(nfi,vpot,drhor,rhog,drhog,rhos,rhoc,tfirst,tlast,&
                     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
         !  
         if ( tefield )  then!to be bettered
            !
            call berry_energy( enb, enbi, becm, cm(:,:), fion )
            etot=etot+enb+enbi
            ! 
         endif
         !  
         if (tefield2 )  then!to be bettered
            !
            call berry_energy2( enb, enbi, becm, cm(:,:), fion )
            etot=etot+enb+enbi
            !
         endif
         ! 
         if (do_orbdep) then
            !
            call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, &
                                  sizwtot, vsic, do_wxd, pink, nudx, &
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            !
            eodd = sum(pink(1:nbsp))
            etot = etot + eodd
            !
         endif
         !  
         enever=etot
         ! 
         if (tens.and.newscheme) then
            ! 
            enever=enever+entropy
            ! 
         endif
         !
         if (tens.and.newscheme) then
            !
            if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0,ene1,enesti,enever
            if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
            ! 
         else
            !
            if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
            if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
            !
         endif
         !
         !check with  what supposed
         !
         if (ionode) then
            !
            if (iverbosity > 1) then
               !
               write(stdout,*) 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
               write(stdout,*) 'cg_sub: minmum   :'  , enever,passo,passov
               !
            endif
            !
         endif
         !
         !if the energy has diminished with respect to ene0 and ene1, everything ok
         !  
         if (((enever.lt.ene0) .and. (enever.lt.ene1)).or.(tefield.or.tefield2)) then
            ! 
            c0(:,:) =cm(:,:)
            bec(:,:)=becm(:,:)
            ene_ok  =.true.
            !
         elseif ((enever.ge.ene1) .and. (enever.lt.ene0)) then
            !
            if (ionode) then
               ! 
               write(stdout,*) 'cg_sub: missed minimum, case 1, iteration',itercg
               !
            endif
            !
            c0=c0+spasso*passov*hi
            restartcg=.true.
            !
            call orthogonalize_wfc_only(c0, bec) 
            ! 
            ene_ok=.false.
            !
            !if  ene1 << energy <  ene0; go to  ene1
            !
         elseif ((enever.ge.ene0).and.(ene0.gt.ene1)) then
            !  
            if (ionode) then
               !
               write(stdout,*) 'cg_sub: missed minimum, case 2, iteration',itercg
               !
            endif  
            !
            c0=c0+spasso*passov*hi
            restartcg=.true.!ATTENZIONE
            !         
            call orthogonalize_wfc_only(c0,bec) 
            !call calbec(1,nsp,eigr,c0,bec)
            !call gram_bgrp( betae, bec, nhsa, c0, ngw )
            ene_ok=.false.
            ! 
            !if ene > ene0,en1 do a steepest descent step
         elseif ((enever.ge.ene0).and.(ene0.le.ene1)) then
            ! 
            if (ionode) then
               !
               write(stdout,*) 'cg_sub: missed minimum, case 3, iteration',itercg
               !
            endif
            !
            iter3=0
            do while (enever.gt.ene0 .and. iter3.lt.maxiter3)
               !
               iter3=iter3+1
               passov=passov*0.5d0
               cm=c0+spasso*passov*hi
               !
               ! change the searching direction
               ! 
               spasso=spasso*(-1.d0)
               ! 
               call orthogonalize_wfc_only(cm,becm)
               !
               itercgeff = itercgeff + 1
               !call calbec(1,nsp,eigr,cm,becm)
               !call gram_bgrp( betae, bec, nhsa, cm, ngw )
               !call calbec(1,nsp,eigr,cm,becm)
               ! 
               if (.not.tens) then
                  !
                  call rhoofr(nfi,cm(:,:),irb,eigrb,becm,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
                  !
               else
                  !
                  if (newscheme)  then
                     ! 
                     call inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,dbec,.false., vpot  )
                     !  
                  endif
                  !
                  ! calculation of the rotated quantities
                  !
                  call rotate( nrlx, z0t, cm(:,:), becm, c0diag, becdiag )
                  ! 
                  ! calculation of rho corresponding to the rotated wavefunctions
                  !
                  call rhoofr(nfi,c0diag,irb,eigrb,becdiag,dbec,rhovan,rhor, &
                         drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
                  !
               endif
               !  
               ! calculates the potential
               !
               ! put core charge (if present) in rhoc(r)
               !
               if (nlcc_any) call set_cc(irb,eigrb,rhoc)
               !
               vpot = rhor
               !
               call vofrho(nfi,vpot,drhor,rhog,drhog,rhos,rhoc,tfirst,tlast, &
                           ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
               !
               if (tefield) then !to be bettered
                  !
                  call berry_energy( enb, enbi, becm, cm(:,:), fion )
                  etot=etot+enb+enbi
                  !
               endif
               !
               if (tefield2)  then !to be bettered
                  ! 
                  call berry_energy2( enb, enbi, becm, cm(:,:), fion )
                  etot=etot+enb+enbi
                  !
               endif
               !
               if (do_orbdep) then
                  !
                  call nksic_potential( nbsp, nbspx, cm, fsic, becm, rhovan, deeq_sic, &
                                        ispin, iupdwn, nupdwn, rhor, rhoc, wtot, &
                                        sizwtot, vsic, do_wxd, pink, nudx, &
                                        wfc_centers, wfc_spreads, &
                                        icompute_spread, .false. )
                  !
                  eodd = sum(pink(1:nbsp))
                  !
                  etot = etot + eodd
                  !
               endif
               ! 
               enever=etot
               !
               if (tens.and.newscheme) then
                  !
                  enever=enever+entropy
                  !
               endif
               !
            enddo
            !
            if (iter3 == maxiter3) write(stdout,*) 'missed minimun: iter3 = maxiter3'
            !
            c0(:,:)=cm(:,:)
            restartcg=.true.
            ene_ok=.false.
            ! 
         endif
         !        
         if (tens.and.newscheme) enever=enever-entropy
         !  
         if (.not. ene_ok) call calbec (1,nsp,eigr,c0,bec)
         !
         ! calculates phi for pc_daga
         !
         call calphi_bgrp( c0, SIZE(c0,1), bec, nhsa, betae, phi, nbsp )
         !    
         !=======================================================================
         !
         !                 start of the inner loop
         !                 (Uij degrees of freedom)
         !
         !=======================================================================
         if (tens.and. .not.newscheme) then
            !
            call inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                    rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,dbec,firstiter, vpot  )
            ! 
            ! the following sets up the new energy
            !
            enever=etot
            !
         endif
         ! 
         !=======================================================================
         !                 end of the inner loop
         !=======================================================================
         !
         itercg=itercg+1
         !
         itercgeff=itercgeff+1 ! including the damping steps in miss case 3 
         !
      enddo!on conjugate gradient iterations
      !
      ! calculates atomic forces and lambda
      !
      ! if pressure is need the following is written because of caldbec
      !
      if (tpre) then
         !   
         call calbec(1,nsp,eigr,c0,bec)
         !
         if (.not.tens) then
            !
            call caldbec_bgrp( eigr, c0, dbec, descla )
            call rhoofr(nfi,c0(:,:),irb,eigrb,bec,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
            !   
         else
            !
            ! calculation of the rotated quantities
            !
            call rotate( nrlx, z0t, c0(:,:), bec, c0diag, becdiag )
            ! 
            ! calculation of rho corresponding to the rotated wavefunctions
            !
            call caldbec_bgrp( eigr, c0diag, dbec, descla )
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag,dbec,rhovan,rhor, &
                        drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6)
            !
         endif
         !
         ! calculates the potential
         !
         ! put core charge (if present) in rhoc(r)
         !
         if (nlcc_any) call set_cc(irb,eigrb,rhoc)
         !
         !---ensemble-DFT
         ! 
         vpot = rhor
         !  
         call vofrho(nfi,vpot,drhor,rhog,drhog,rhos,rhoc,tfirst,tlast,&
                     ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
         !
         if (do_orbdep) then
            !
            call nksic_potential( nbsp, nbspx, c0, fsic, bec, rhovan, deeq_sic, &
                                  ispin, iupdwn, nupdwn, rhor, rhoc, wtot, &
                                  sizwtot, vsic, do_wxd, pink, nudx, &
                                  wfc_centers, wfc_spreads, &
                                  icompute_spread, .false.)
            !
            eodd = sum(pink(1:nbsp))
            etot = etot + eodd
            !
         endif
         ! 
      endif
      !
      call calcmt( nrlx, f, z0t, fmat0 )
      !
      call newd(vpot,irb,eigrb,rhovan,fion)
      !
      if (.not.tens) then
         !
         if (tfor .or. tprnfor) call nlfq_bgrp( c0, eigr, bec, becdr, fion ) ! call nlfq(c0,eigr,bec,becdr,fion)
         !
      else
         !
         if (tfor .or. tprnfor) call nlfq_bgrp( c0diag, eigr, becdiag, becdrdiag, fion ) ! call nlfq(c0diag,eigr,becdiag,becdrdiag,fion)
         ! 
      endif
      ! 
      call prefor(eigr,betae)
      ! 
      ! faux takes into account spin multiplicity.
      !
      faux(1:nbsp) = max(f_cutoff, f(1:nbsp))
      ! 
      do i=1,nbsp,2
         !
         call dforce(i,bec,betae,c0,c2,c3,rhos,dffts%nnr,ispin,faux,nbsp,nspin)
         ! 
         if (tefield.and.(evalue .ne. 0.d0)) then
            !
            call dforceb &
               (c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            !
            do ig=1,ngw
               !
               c2(ig)=c2(ig)+evalue*df(ig)
               !
            enddo
            !
            call dforceb &
               (c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            !
            do ig=1,ngw
               !
               c3(ig)=c3(ig)+evalue*df(ig)
               !
            enddo
            !
         endif
         ! 
         if (tefield2.and.(evalue2 .ne. 0.d0)) then
            !  
            call dforceb &
               (c0, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            !
            do ig=1,ngw
               ! 
               c2(ig)=c2(ig)+evalue2*df(ig)
               !
            enddo
            ! 
            call dforceb &
               (c0, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            !
            do ig=1,ngw
               ! 
               c3(ig)=c3(ig)+evalue2*df(ig)
               ! 
            enddo
            !
         endif
         !
         if ( do_orbdep ) then
            !
            ! faux takes into account spin multiplicity.
            !
            CALL nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi )
            !
            c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
            !
            if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
            !
         endif
         ! 
         do ig=1,ngw
            !
            gi(ig,i)=c2(ig)
            !
            if (i+1 <= nbsp) then
               !
               gi(ig,i+1)=c3(ig)
               !
            endif
            !
         enddo
         !
         if (gstart==2) then
            ! 
            gi(1, i)=CMPLX(DBLE(gi(1,i)),0.d0,kind=DP)
            !
            if (i+1 <= nbsp) then
               !
               gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0,kind=DP)
               !
            endif
            !
          endif
          !
       enddo
       !
       allocate( lambda_repl( nudx, nudx ) )
       !
       do is = 1, nspin
          !
          nss = nupdwn(is)
          istart = iupdwn(is)
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
                       2.d0 * DBLE( CONJG( c0( ig, ii ) ) * gi( ig, jj) )
                   ! 
                enddo
                !
                if (gstart == 2) then
                   ! 
                   lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( c0( 1, ii ) ) * gi( 1, jj ) )
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
          !
          call distribute_lambda( lambda_repl, lambda( :, :, is ), descla( is ) )
          !
       enddo
       !
       if (l_cprestart .and. .not.tens .and. nspin==1 .and. nvb<1 .and. .not.do_orbdep) then
          ! 
          !if required project c0 on previous manifold of occupied states                                                                                    
          !
          !NOT IMPLEMENTED YET FOR ENSEMBLE DFT AND NSPIN==2
          !NOT IMPLEMENTED FOR US PSEUDOPOTENTIALS
          !
          lambda_repl=0.d0
          !  
          do i = 1, nss
             !
             do j = 1, nss
                !
                ii = i + istart - 1
                jj = j + istart - 1
                !  
                do ig = 1, ngw
                   ! 
                   lambda_repl( i, j ) = lambda_repl( i, j ) + &
                          2.d0 * DBLE( CONJG( c0old( ig, ii ) ) * c0( ig, jj) )
                   ! 
                enddo
                ! 
                if (gstart == 2) then
                   !
                   lambda_repl( i, j ) = lambda_repl( i, j ) - &
                          DBLE( CONJG( c0old( 1, ii ) ) * c0( 1, jj ) )
                   !
                endif
                !
             enddo
             !
          enddo
          !
          call mp_sum( lambda_repl, intra_bgrp_comm )
          ! 
          cm(:,:)=c0(:,:)
          c0=(0.d0,0.d0)
          !
          do i=1,nss
             ! 
             do j=1,nss
                !
                c0(1:ngw,i)=c0(1:ngw,i)+lambda_repl(i,j)*cm(1:ngw,j)
                !
             enddo
             ! 
          enddo
          ! 
          call calbec (1,nsp,eigr,c0,bec)
          call gram_bgrp( betae, bec, nhsa, c0, ngw )
          call calbec(1,nsp,eigr,c0,bec)
          !          
          do i=1,nbsp,2
             !
             call dforce(i,bec,betae,c0,c2,c3,rhos,dffts%nnr,ispin,f,nbsp,nspin)
             ! 
             do ig=1,ngw
                !
                gi(ig,  i)=c2(ig)
                !
                if (i+1 <= nbsp) then
                   !
                   gi(ig,i+1)=c3(ig)
                   !
                endif
                !
             enddo
             ! 
             if (gstart==2) then
                !  
                gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0,kind=DP)
                !
                if (i+1 <= nbsp) then
                   !
                   gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0,kind=DP)
                   !
                endif
                !
             endif
             ! 
          enddo
          ! 
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
                          2.d0 * DBLE( CONJG( c0( ig, ii ) ) * gi( ig, jj) )
                enddo
                !
                if ( gstart == 2 ) then
                   !
                   lambda_repl( i, j ) = lambda_repl( i, j ) + &
                          DBLE( CONJG( c0( 1, ii ) ) * gi( 1, jj ) )
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
          !
          call distribute_lambda( lambda_repl, lambda( :, :, 1 ), descla( 1 ) )
          ! 
          cm(:,:)=c0(:,:)
          call calbec (1,nsp,eigr,cm,becm)
          !
       endif
       !
       deallocate ( lambda_repl )
       ! 
       if ( tens ) then
          !
          ! in the ensemble case matrix labda must be multiplied with f
          ! 
          allocate ( lambda_dist( nlam, nlam ) )
          ! 
          do iss = 1, nspin
             !
             nss    = nupdwn( iss )
             !
             lambdap(:,:,iss) = 0.0d0
             !
             call cyc2blk_redist( nss, fmat0(1,1,iss), nrlx, SIZE(fmat0,2), &
                                  lambda_dist, nlam, nlam, descla(iss) )
             !
             ! Perform lambdap = lambda * fmat0
             !
             call sqr_mm_cannon( 'N', 'N', nss, 1.0d0, lambda(1,1,iss), nlam, lambda_dist, nlam, &
                                  0.0d0, lambdap(1,1,iss), nlam, descla(iss) )
             !
             lambda_dist      = lambda(:,:,iss)
             lambda(:,:,iss)  = lambdap(:,:,iss)
             lambdap(:,:,iss) = lambda_dist
             !
          enddo
          !
          deallocate ( lambda_dist )
          !
          call nlsm2_bgrp( ngw, nhsa, eigr, c0, becdr, nbspx, nbsp )
          !
       endif
       !
       call nlfl_bgrp( bec, becdr, lambda, descla, fion ) 
       !
       ! bforceion adds the force term due to electronic berry phase
       ! only in US-case
       ! 
       if (tefield.and.(evalue .ne. 0.d0)) then
          ! 
          call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)
          ! 
       endif
       !
       if (tefield2.and.(evalue2 .ne. 0.d0) ) then
          !
          call bforceion(fion,tfor.or.tprnfor,ipolp2, qmat2,bec,becdr,gqq2,evalue2)
          !  
       endif
       !
       deallocate(hpsi0,hpsi,gi,hi)
       deallocate(s_minus1,k_minus1)
       !
       if (ionode) close(37)!for debug and tuning purposes
       !
       call stop_clock('runcg_uspp')
       ! 
       deallocate(bec0,becm,becdrdiag)
       deallocate(ave_ene)
       deallocate(c2,c3)
       deallocate(faux)
       !
       return
       !
  contains
    !
    subroutine compute_hpsi()
       ! 
       ! faux takes into account spin multiplicity.
       !
       faux(1:nbspx)=0.d0
       faux(1:nbsp) = max(f_cutoff,f(1:nbsp)) 
       !   
       do i=1,nbsp,2
          ! 
          call dforce( i, bec, betae, c0,c2,c3,rhos, dffts%nnr, ispin,faux,nbsp,nspin)
          !
          if (tefield .and. (evalue.ne.0.d0)) then
             ! 
             call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
             c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw) 
             call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
             c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
             !
          endif
          ! 
          if (tefield2 .and. (evalue2.ne.0.d0)) then
             !
             call dforceb(c0, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
             c2(1:ngw)=c2(1:ngw)+evalue2*df(1:ngw) 
             call dforceb(c0, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
             c3(1:ngw)=c3(1:ngw)+evalue2*df(1:ngw)
             !
          endif
          !
          if ( do_orbdep ) then
             !
             ! faux takes into account spin multiplicity.
             !
             call nksic_eforce( i, nbsp, nbspx, vsic, deeq_sic, bec, ngw, c0(:,i), c0(:,i+1), vsicpsi )
             !
             c2(:) = c2(:) - vsicpsi(:,1) * faux(i)
             !
             if( i+1 <= nbsp )   c3(:) = c3(:) - vsicpsi(:,2) * faux(i+1)
             !
          endif
          !
          hpsi(1:ngw,i)=c2(1:ngw)
          !
          if (i+1 <= nbsp) then
             !
             hpsi(1:ngw,i+1)=c3(1:ngw)
             !
          endif
          !
          if (gstart==2) then
             ! 
             hpsi(1,i)=CMPLX(DBLE(hpsi(1,i)), 0.d0,kind=DP)
             ! 
             if (i+1 <= nbsp) then
                ! 
                hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0,kind=DP)
                !
             endif
             !
          endif  
          !
       enddo
       !
    endsubroutine compute_hpsi
    !
    subroutine orthogonalize(wfc0, wfc, becwfc, bec0)
       ! 
       real(dp) :: becwfc(:,:), bec0(:,:)
       complex(dp) :: wfc(:,:), wfc0(:,:)
       ! 
       if (switch.or.(.not.do_orbdep)) then
          !
          !if (fixed_state) then
            !
            !wfc(:, fixed_band) = cmplx(0.d0, 0.d0)
            ! 
          !endif
          !
          call calbec(1,nsp,eigr,wfc,becwfc)
          !  
          call pc2(wfc0,bec0,wfc,becwfc)
          !
       else
          !
          !if (fixed_state) then
            !
            !call pc3nc_fixed(wfc0, wfc, lgam)
            !
          !endif
          !   
          call calbec(1,nsp,eigr,wfc,becwfc)
          ! 
          if (.not.(nvb>0)) then
             !
             call pc2nc_nksic(wfc0,wfc, nbsp, ispin)
             !
          else
             ! 
             call pc2(wfc0,bec0,wfc,becwfc)
             !
          endif
          !
       endif
       !  
    endsubroutine orthogonalize
    !
    subroutine orthogonalize_wfc_only(wfc,becwfc)
        !
        real(DP)    :: becwfc(:,:)
        complex(DP) :: wfc(:,:)
        !
        call calbec(1,nsp,eigr,wfc,becwfc)
        !
           !if (fixed_state) then
             !
           !  call gram_swap(betae,becwfc,nhsa,wfc,ngw,nbsp, fixed_band)
             ! 
           !else
             ! 
        call gram_bgrp( betae, becwfc, nhsa, wfc, ngw ) 
             ! 
        !   endif
            !
        !ENDIF
   endsubroutine orthogonalize_wfc_only
   !
   subroutine do_innerloop_subroutine()
        !
        if (do_innerloop) then ! .and. innerloop_until>=itercgeff) then
           !
           call start_clock( "inner_loop" )
           !
           eodd    = sum(pink(1:nbsp))
           etot    = etot - eodd
           etotnew = etotnew - eodd
           ninner  = 0
           ! 
           !if (.not.do_innerloop_cg) then
           !   call nksic_rot_emin(itercg,ninner,etot,Omattot, lgam)
           !else
           !call
           !nksic_rot_emin_cg(itercg,innerloop_init_n,ninner,etot,Omattot,deltae*innerloop_cg_ratio,lgam)
           call nksic_rot_emin_cg_general(itercg, innerloop_init_n, ninner, etot, deltae*innerloop_cg_ratio, &
                                      nbsp, nbspx, nudx, iupdwn, nupdwn, ispin, c0, rhovan, bec, rhor, rhoc, &
                         vsic, pink, deeq_sic, wtot, fsic, sizwtot, do_wxd, wfc_centers, wfc_spreads, .false.)
           !
           ! endif
           !  
           eodd    = sum(pink(1:nbsp))
           etot    = etot + eodd
           etotnew = etotnew + eodd
           eoddnew = eodd
           !
           call stop_clock( "inner_loop" )
           !
        endif
        !
        eodd = sum(pink(1:nbsp))
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
              !
              write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4," Etot (Ha) =",F22.14," delta_E=",E22.14)')&
              itercg, itercgeff, etotnew, deltae
              ! 
           else
              !  
              write(stdout,'(5x,"iteration =",I4,"  eff iteration =",I4," Etot (Ha) =",F22.14)')&
              itercg, itercgeff, etotnew
              !  
           endif
           !
        endif
        !
        !if(ionode) write(stdout, *) 'itercg, etotnew,pberryel,pberryel2'!for
        !debug and tuning purposes
        !if(ionode) write(stdout, *) itercg, etotnew,pberryel,pberryel2!for
        !debug and tuning purposes
        !
   endsubroutine print_out_observables
   !    
endsubroutine kcp_runcg_uspp 
