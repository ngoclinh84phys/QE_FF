!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine nksic_potential( nbsp, nx, c, f_diag, bec, becsum, &
                            deeq_sic, ispin, iupdwn, nupdwn,  &
                            rhor, rhoc, wtot, sizwtot, vsic, do_wxd_, pink, nudx, &
                            wfc_centers, wfc_spreads, icompute_spread, is_empty)
!-----------------------------------------------------------------------
      !
      ! ....calculate orbital dependent potentials 
      !     Perdew-Zunger (PZ),
      !     Koopmans' integral definition (KI)
      !     and KIPZ methods
      !
      use kinds,                      only: dp
      use io_global,                  only: stdout
      use gvect,                      only: ngm
      use gvecs,                      only: ngms
      use gvecw,                      only: ngw
      use fft_base,                   only: dfftp
      USE electrons_base,             only: nspin
      use funct,                      only: dft_is_gradient, dft_is_meta
      use nksic,                      only: do_nki, do_pz, do_nkpz, do_nkipz, do_nk, &
                                            fion_sic, valpsi, odd_alpha, nkscalfact, odd_nkscalfact
      use ions_base,                  only: nat
      use control_flags,              only: gamma_only
      use uspp,                       only: nkb, nlcc_any
      use uspp_param,                 only: nhm
      !
      implicit none
      !
      ! in/out vars
      !
      integer,     intent(in)  :: nbsp, nx, nudx, sizwtot
      integer,     intent(in)  :: ispin(nx)
      integer,     intent(in)  :: iupdwn(nspin), nupdwn(nspin)
      real(dp),    intent(in)  :: bec(nkb,nbsp)
      real(dp),    intent(in)  :: becsum( nhm*(nhm+1)/2, nat, nspin)
      real(dp),    intent(in)  :: f_diag(nx)
      real(dp)                 :: rhor(dfftp%nnr,nspin)
      real(dp),    intent(in)  :: rhoc(dfftp%nnr)
      real(dp),    intent(out) :: vsic(dfftp%nnr,nx), wtot(sizwtot,2)
      real(dp),    intent(out) :: deeq_sic(nhm,nhm,nat,nx)
      real(dp),    intent(out) :: pink(nx)
      complex(dp), intent(in)  :: c(ngw,nx)
      logical,     intent(in)  :: do_wxd_
      logical  :: icompute_spread
      real(DP) :: wfc_centers(4,nudx,nspin)
      real(DP) :: wfc_spreads(nudx,nspin,2)
      logical  :: is_empty
      !
      ! local variables
      !
      integer  :: i,j,jj,ibnd,isp,ir
      real(dp) :: focc,pinkpz,shart
      real(dp), allocatable    :: vsicpz(:)
      real(dp), allocatable    :: rhor_nocc(:,:)
      !
      complex(dp), allocatable :: rhobarg(:,:)
      real(dp), allocatable    :: grhobar(:,:,:)
      !
      real(dp), allocatable    :: rhobar(:,:)
      real(dp), allocatable    :: rhoref(:,:)
      real(dp), allocatable    :: orb_rhor(:,:)
      real(dp), allocatable    :: wxdsic(:,:)
      real(dp), allocatable    :: wrefsic(:)
      real(dp), allocatable    :: vxc_aux(:,:)
      real(dp) :: etxc_aux
      ! 
      logical  :: is_empty_
      !
      ! main body
      !
      call start_clock( 'nksic_potential' )
      !
      is_empty_=is_empty
      !
      ! initialization, should put below to subrotuine 
      !
      if (dft_is_gradient()) then
         !
         allocate(rhobarg(ngm,2))
         allocate(grhobar(dfftp%nnr,3,2))
         !
      else
         !
         allocate(rhobarg(1,1))
         allocate(grhobar(1,1,1))
         !  
      endif
      !
      allocate( orb_rhor(dfftp%nnr,2) )
      ! 
      if (do_nki .or. do_nkipz .or. do_nk .or. do_nkpz) then
         ! 
         allocate( rhobar(dfftp%nnr,2) )
         allocate( rhoref(dfftp%nnr,2) )
         !
         if ( do_nk .or. do_nkpz ) then
            ! 
            allocate( wxdsic(dfftp%nnr,2) )
            allocate( wrefsic(dfftp%nnr) )
            wxdsic=0.0_dp
            wrefsic=0.0_dp
            !
         endif
         ! 
         if ( do_nki .or. do_nkipz ) then
            ! 
            allocate( wxdsic(dfftp%nnr, 2) )
            wxdsic=0.0_dp
            ! 
         endif
         !
         if ( do_nkpz .or. do_nkipz) then
            ! 
            allocate(vsicpz(dfftp%nnr))
            vsicpz=0.0_dp
            !
         endif
         !
         if (do_nki .or. do_nk .or. do_nkipz) then
            !
            allocate(vxc_aux(dfftp%nnr, 2))
            vxc_aux=0.0_dp
            etxc_aux=0.0_dp
            !
         endif
         !  
      endif
       !
      if (nlcc_any) then
         !
         allocate(rhor_nocc(dfftp%nnr,nspin))
         rhor_nocc(:,:) = rhor(:,:)
         !
         ! add core charge
         !
         call add_cc_rspace(rhoc, rhor)
         !
      endif
      !
      if ( do_nki .or. do_nkipz .or. do_nk .or. do_nkpz ) then
         ! 
         wtot=0.0_dp
         !
      endif
      !
      !
      pink=0.0_dp
      vsic=0.0_dp
      !
      ! done initialization
      !
      ! loop over bands (2 ffts at the same time)
      !
      do j=1, nbsp, 2
         !
         ! compute orbital densities
         ! n odd => c(:,n+1) is already set to zero
         !
         call nksic_get_orbitalrho( bec, ispin, nbsp, &
                                    c(:,j), c(:,j+1), orb_rhor, j, j+1 )
         !
         ! compute centers and spreads of nksic or pz
         ! minimizing orbitals
         !
         !if (icompute_spread) then
            !
         !   call compute_nksic_centers(dfftp%nnr, nx, nudx, nbsp, nspin, iupdwn, &
         !            nupdwn, ispin, orb_rhor, wfc_centers, wfc_spreads, j, j+1)
         !   !
         !endif
         !
         shart=0.d0
         !
         ! compute orbital potentials
         !
         inner_loop: do jj=1,2
           !
           i=j+jj-1
           !
           ! this condition is important when n is odd
           !
           if ( i > nbsp ) exit inner_loop
           !
           ibnd=i
           !
           if ( nspin==2 ) then
              !
              if ( i >= iupdwn(2) ) ibnd=i-iupdwn(2)+1
              !
           endif
           !
           ! note: iupdwn(2) is set to zero if nspin = 1
           !
           focc=f_diag(i)*DBLE(nspin)/2.0d0
           !
           if (do_nki .or. do_nkipz .or. do_nk .or. do_nkpz ) then
              !
              ! define rhoref and rhobar
              !
              call nksic_get_rhoref( i, dfftp%nnr, ispin(i), nspin,  &
                                     focc, rhor, orb_rhor(:,jj), &
                                     rhoref, rhobar, rhobarg, grhobar )
           endif
           !
           ! compute pz potentials and energy
           !
           if ( do_pz ) then
              !
              call nksic_correction_pz ( focc, ispin(i), orb_rhor(:,jj), &
                                         vsic(:,i), pink(i), ibnd, shart )
              !
         !     wfc_spreads(ibnd, ispin(i), 2) = shart
              !
           endif
           !
           ! compute nki pieces to build the potentials and the energy
           !
           if ( do_nki .or. do_nkipz) then
              !
              call nksic_correction_nki( focc, ispin(i), orb_rhor(:,jj), &
                                         rhor, rhoref, rhobar, rhobarg, grhobar, &
                                         vsic(:,i), wxdsic, do_wxd_, pink(i), ibnd, shart, vxc_aux, etxc_aux, is_empty_)
              !
              ! here information is accumulated over states
              ! (wtot is added in the next loop)
              !
              wtot(:,1:2) = wtot(:,1:2) + wxdsic(:,1:2)
              !
              ! the sic potential is partly updated here to save some memory
              !
              vsic(:,i) = vsic(:,i) - wxdsic(:, ispin(i) )
              !
              !wfc_spreads(ibnd, ispin(i), 2) = shart
              !
           endif
           !
           if ( do_nkipz ) then
              !
              call nksic_correction_nkipz( focc, ispin(i), orb_rhor(:,jj), vsicpz, &
                                           pinkpz, ibnd, shart, is_empty_ )
              !
              vsic(:,i) = vsic(:,i) + vsicpz(:)
              !
              pink(i) = pink(i) + pinkpz
              !
              !wfc_spreads(ibnd, ispin(i), 2) = shart
              !
           endif
           !
           ! take care of spin symmetry, Linh ???? 
           !
           if (.not.is_empty_) then
              !
              pink(i) = 2.d0 * pink(i)/nspin
              ! 
           else
              !
              pink(i) = 2.d0 * pink(i)/nspin
              ! 
           endif
           !
           if ( do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
              !
              if ( nspin== 1 ) then
                 !
                 wtot(:,1) = wtot(:,1) + wxdsic(:,2)
                 ! 
                 wtot(:,2) = wtot(:,2) + wxdsic(:,1)
                 !
              endif
              !
           endif
           !
         enddo inner_loop
         !
      enddo
      !
      ! Switch off the icompute_spread flag if present
      !
      if (icompute_spread) then
         !
         icompute_spread=.false.
         !
      endif
      !
      ! now wtot is completely built and can be added to vsic
      !
      if ( do_nk .or. do_nkpz .or. do_nki .or. do_nkipz ) then
         !
         do i = 1, nbsp
            !
            vsic(:,i) = vsic(:,i) + wtot(:, ispin(i))
            !
         enddo
         !
      endif
      !
      ! computing orbital dependent alpha
      !
      if ( odd_nkscalfact ) then
         !
         do j=1,nbsp,2
            !
            inner_loop_odd_alpha: do jj=1,2
               !
               i=j+jj-1
               !  
               if ( i > nbsp ) exit inner_loop_odd_alpha
               !
      !         vsic(1:nnrx,i) = vsic(1:nnrx,i)*odd_alpha(i)/nkscalfact
               !
      !         valpsi(i,:) = valpsi(i,:) * pink(i)/nkscalfact
               ! 
      !         pink(i) = pink(i)*odd_alpha(i)/nkscalfact
               !
            enddo inner_loop_odd_alpha
            !
         enddo
         !
      endif
      !
      ! USPP:
      ! compute corrections to the D coefficients of the pseudopots
      ! due to vsic(r, i) in the case of orbital dependent functionals.
      ! The corresponding contributions to the forces are computed.
      !
      ! IMPORTANT: the following call makes use of newd.
      !            It must be done before we call newd for the
      !            total potentials, because deeq is overwritten at every call
      !
      fion_sic(:,:) = 0.0d0
      !
      if ( nhm > 0 ) then
         !
         deeq_sic(:,:,:,:) = 0.0d0
         !
         do i = 1, nbsp
            !
       !     CALL nksic_newd( i, nnrx, ispin(i), nspin, vsic(:,i), nat, nhm, &
       !                      becsum, fion_sic, deeq_sic(:,:,:,i) ) 
            !this is for ultrasoft! watch out! warning:giovanni 
            !this has to be modified in order to run ultrasoft
            !
         enddo
         !
      endif
      !
      !
      if (nlcc_any) then
         !
         rhor(:,:)=rhor_nocc(:,:)
         deallocate(rhor_nocc)
         !
      endif
      !
      deallocate(rhobarg)
      deallocate(grhobar)
      deallocate(orb_rhor)
      if(allocated(rhobar)) deallocate(rhobar)
      if(allocated(rhoref)) deallocate(rhoref)
      if(allocated(wxdsic)) deallocate(wxdsic)
      if(allocated(wrefsic)) deallocate(wrefsic)
      if(allocated(vsicpz)) deallocate(vsicpz)
      if(allocated(vxc_aux)) deallocate(vxc_aux)
      !
      !
      call stop_clock('nksic_potential')
      ! 
      return
      !
endsubroutine nksic_potential

subroutine nksic_get_orbitalrho( bec, ispin, nbsp, &
                                 c1, c2, orb_rhor, i1, i2 )
      !
      ! Computes orbital densities on the real (not smooth) grid
      !
      use kinds,                      only: dp
      use fft_interfaces,             only: fwfft, invfft
      use fft_base,                   only: dffts, dfftp
      use cp_interfaces,              only: calrhovan
      use cell_base,                  only: omega
      use gvecw,                      only: ngw
      use gvect,                      only: ngm,  nl, nlm
      use gvecs,                      only: ngms, nls, nlsm
      use cp_main_variables,          only: irb, eigrb
      use uspp_param,                 only: nhm
      use electrons_base,             only: nspin
      use ions_base,                  only: nat
      use uspp,                       only: okvan, nkb

      use io_global,          only: stdout
      use mp_global,          only: intra_bgrp_comm
      use mp,                 only: mp_sum

      !
      implicit none
      !
      ! input/output vars
      !
      integer,     intent(in) :: i1,i2
      integer,     intent(in) :: nbsp, ispin(nbsp)
      real(dp),    intent(in) :: bec(nkb, nbsp)
      real(dp),    intent(out):: orb_rhor(dfftp%nnr,2)
      complex(dp), intent(in) :: c1(ngw),c2(ngw)
      !
      ! local vars
      !
      character(20) :: subname='nksic_get_orbitalrho'
      integer       :: ir, ig, ierr, nnrx, nnrsx
      real(dp)      :: sa1
      complex(dp)   :: fm, fp, ci
      complex(dp), allocatable :: psis(:), psi(:)
      complex(dp), allocatable :: orb_rhog(:,:)
      real(dp),    allocatable :: orb_rhos(:)
      real(dp),    allocatable :: rhovan(:,:,:)
      real(dp),    allocatable :: rhovanaux(:,:,:)
      !
      ci = ( 0.0d0, 1.0d0 )
      !
      nnrx  = dfftp%nnr
      nnrsx = dffts%nnr 
      !
      if ( okvan ) then
         !
         allocate(rhovan(nhm*(nhm+1)/2,nat,nspin), stat=ierr )
         if ( ierr/=0 ) call errore(subname,'allocating rhovan',abs(ierr))
         allocate(rhovanaux(nhm*(nhm+1)/2,nat,nspin), stat=ierr)
         if ( ierr/=0 ) call errore(subname,'allocating rhovanaux',abs(ierr))
         !
      endif
      !
      allocate(psi(nnrx),stat=ierr)
      if ( ierr/=0 ) call errore(subname,'allocating psi',abs(ierr))
      !
      allocate(orb_rhog(ngm,2),stat=ierr)
      if ( ierr/=0 ) call errore(subname,'allocating orb_rhog',abs(ierr))
      !
      sa1 = 1.0d0 / omega
      !
      ! check whether it is necessary to
      ! deal with the smooth and dense grids separately
      !
      if ( nnrsx == nnrx ) then
         !
         ! This case should be the case when using NCPP
         !
         call c2psi( psi, nnrx, c1, c2, ngw, 2 )
         !
         call invfft('Dense', psi, dfftp )
         !
         ! computing the orbital charge in real space on the full grid
         !
         do ir = 1, nnrx
            !
            orb_rhor(ir,1) = sa1 * ( DBLE(psi(ir)) )**2
            orb_rhor(ir,2) = sa1 * ( AIMAG(psi(ir)) )**2
            !
         enddo
         !
      else
         !
         ! this is the general case,
         ! normally used with USPP
         !
         allocate( psis(nnrsx), stat=ierr )
         if ( ierr/=0 ) call errore(subname,'allocating psis',abs(ierr))
         allocate( orb_rhos(2), stat=ierr )
         if ( ierr/=0 ) call errore(subname,'allocating orb_rhos',abs(ierr))
         !
         call c2psi( psis, nnrsx, c1, c2, ngw, 2 )
         !
         call invfft('Wave',psis, dffts )
         !
         ! computing the orbital charge
         ! in real space on the smooth grid
         !
         do ir = 1, nnrsx
            !
            orb_rhos(1) = sa1 * ( DBLE(psis(ir)) )**2
            orb_rhos(2) = sa1 * ( AIMAG(psis(ir)) )**2
            !
            psis( ir )  = CMPLX( orb_rhos(1), orb_rhos(2),kind=DP)
            ! 
         enddo
         !
         ! orbital charges are taken to the G space
         !
         CALL fwfft('Smooth',psis, dffts )
         !
         do ig = 1, ngms 
            !
            fp=psis(nls(ig))+psis(nlsm(ig))
            fm=psis(nls(ig))-psis(nlsm(ig))
            orb_rhog(ig,1)=0.5d0*CMPLX(DBLE(fp),AIMAG(fm), kind=DP)
            orb_rhog(ig,2)=0.5d0*CMPLX(AIMAG(fp),-DBLE(fm), kind=DP)
            !
         enddo
         !
         psi (:) = (0.d0, 0.d0)
         do ig=1,ngms
            !
            psi(nlm(ig)) = CONJG(orb_rhog(ig,1)) +ci*CONJG(orb_rhog(ig,2))
            psi(nl (ig)) = orb_rhog(ig,1)        +ci*orb_rhog(ig,2)
            !
         enddo
         !
         call invfft('Dense',psi,dfftp)
         !
         do ir=1,nnrx
            !
            orb_rhor(ir,1) = DBLE(psi(ir))
            orb_rhor(ir,2) = AIMAG(psi(ir))
            !
         enddo
         !
         deallocate( psis )
         deallocate( orb_rhos )
         !
      endif
      !
      ! add Vanderbilt contribution to orbital density
      !
      if( okvan ) then
        !
      !  rhovan(:,:,:) = 0.0d0
        !
       ! if ( nspin == 2 ) then
            !
        !    if ( i1 <= nbsp ) then
         !       call calrhovan(rhovanaux,bec,i1)
         !       rhovan(:,:,1)=rhovanaux(:,:,ispin(i1))
         !   endif
            !
          !  if ( i2 <= nbsp ) then
          !      call calrhovan(rhovanaux,bec,i2)
          !      rhovan(:,:,2)=rhovanaux(:,:,ispin(i2))
          !  endif
            !
          !  call rhov(irb,eigrb,rhovan,orb_rhog,orb_rhor, .true.)
        !else
            !
         !   if ( i1 <= nbsp ) then
          !      call calrhovan(rhovanaux,bec,i1)
           !     rhovan(:,:,1)=rhovanaux(:,:,ispin(i1))*0.5d0 ! 1/2 factor since rhovanaux is counted twice in the case nspin=2
            !    !
             !   call rhov(irb,eigrb,rhovan,orb_rhog(:,1),orb_rhor(:,1), .true.)
              !  !
            !endif
            !
            !if ( i2 <= nbsp ) then
            !    call calrhovan(rhovanaux,bec,i2)
            !    rhovan(:,:,1)=rhovanaux(:,:,ispin(i2))*0.5d0 ! 1/2 factor since rhovanaux is counted twice in the case nspin=2
            !    !
            !    call rhov(irb,eigrb,rhovan,orb_rhog(:,2),orb_rhor(:,2), .true.)
                !
           ! endif
            !
       ! endif
        !
      endif
      !
      deallocate(psi)
      deallocate(orb_rhog)
      !
      if ( okvan ) then
          deallocate(rhovan)
          deallocate(rhovanaux)
      endif
      !
      return
      !
end subroutine nksic_get_orbitalrho
!
!
!
subroutine nksic_get_rhoref( i, nnrx, ispin, nspin, f, &
              rhor, orb_rhor, rhoref_, rhobar_,rhobarg, grhobar_)
      !
      ! Computes rhoref and rhobar
      !
      use io_global,          only: stdout
      use kinds,                      only : dp
      use gvect,                      only : ngm
      use funct,                      only : dft_is_gradient
      use fft_interfaces,             only : fwfft, invfft
      use cp_interfaces,              only : fillgrad
      use fft_base,                   only : dffts, dfftp
      use gvecs,                      only : ngms, nls, nlsm
      use nksic,                      only : fref, rhobarfact
      !
      implicit none
      !
      ! input/output vars
      !
      integer,       intent(in)  :: i, nnrx
      integer,       intent(in)  :: ispin, nspin
      real(dp),      intent(in)  :: f
      real(dp),      intent(in)  :: rhor(nnrx,nspin)
      real(dp),      intent(in)  :: orb_rhor(nnrx)
      real(dp),      intent(out) :: rhoref_(nnrx,2)
      real(dp),      intent(out) :: rhobar_(nnrx,2)
      complex(dp)                :: rhobarg(ngm,2)
      real(dp),      intent(out) :: grhobar_(nnrx,3,2)
      !
      ! internal variables
      !
      integer      :: ig
      complex(dp)  :: fp, fm
      complex(dp),   allocatable :: psi(:)
      !
      ! main body
      !
      call start_clock( 'nksic_get_rhoref' )
      !
      ! define rhobar_i = rho - f_i * rho_i
      !
      if ( nspin == 1 ) then
         !
         rhobar_(:,1) = rhor(:,1) * 0.5_dp
         rhobar_(:,2) = rhor(:,1) * 0.5_dp
         !
      else
         !
         rhobar_(:,1:2) = rhor(:,1:2)
         !
      endif
      !
      rhobar_(:,ispin) = rhobar_(:,ispin) - f * orb_rhor(:)
      !
      ! probably obsolete
      !
      if ( rhobarfact < 1.0d0 ) then
         !
         rhobar_ = rhobar_ * rhobarfact
         !
      endif
      !
      ! define rhoref = rho + (f_ref -f_i) rho_i = rhobar_i + f_ref * rho_i
      ! build rhoref from scratch
      !
      rhoref_(:,1:2)   = rhobar_(:,1:2)
      rhoref_(:,ispin) = rhoref_(:,ispin) + fref * orb_rhor(:)
      !
      ! compute the gradient of rhobar, if needed
      !
      if ( dft_is_gradient() ) then
         !
         allocate( psi(nnrx) )
         !
         psi(:) = CMPLX ( rhobar_(:,1), rhobar_(:,2), kind=DP )
         !
         call fwfft('Dense', psi, dfftp )
         !
         do ig=1,ngm
            ! 
            fp = psi( nls(ig) ) +psi( nlsm(ig) )
            fm = psi( nls(ig) ) -psi( nlsm(ig) )
            !
            rhobarg(ig,1) = 0.5d0 *CMPLX( DBLE(fp),AIMAG(fm), kind=DP)
            rhobarg(ig,2) = 0.5d0 *CMPLX(AIMAG(fp),-DBLE(fm), kind=DP)
            ! 
         enddo
         !
         call fillgrad( 2, rhobarg, grhobar_)
         !
         deallocate( psi )
         !
      endif
      !
      call stop_clock( 'nksic_get_rhoref' )
      !
      return
      !
endsubroutine nksic_get_rhoref
!
!
subroutine add_cc_rspace( rhoc, rhor )
      !
      ! add core correction to the charge density for exch-corr calculation
      ! this subroutine performs the addition in r-space only
      ! see also in nlcc.f90 for g-space 
      !
      use kinds,              only: dp
      use electrons_base,     only: nspin
      use control_flags,      only: iverbosity
      use io_global,          only: stdout
      use mp_global,          only: intra_bgrp_comm
      use cell_base,          only: omega
      use mp,                 only: mp_sum
      use gvect,              only: gstart, ngm, nl
      use fft_interfaces,     only: fwfft
      use fft_base,           only: dfftp
      !
      implicit none
      !
      real(dp), intent(in)   :: rhoc( dfftp%nnr )
      real(dp), intent(inout):: rhor( dfftp%nnr, nspin )
      !
      integer :: ig, ir, iss, isup, isdw
      real(dp):: rsum
      !
      if ( iverbosity > 1 ) then
         !   
         rsum = sum( rhoc ) * omega / dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         call mp_sum( rsum, intra_bgrp_comm )
         write( stdout, 10 ) rsum
10       format( 3X, 'Core Charge = ', D14.6 )
         ! 
      endif
      !
      ! In r-space:
      !
      if ( nspin .eq. 1 ) then
         !
         iss=1
         call daxpy(dfftp%nnr,1.d0,rhoc,1,rhor(1,iss),1)
         !  
      else
         ! 
         isup=1
         isdw=2
         call daxpy(dfftp%nnr,0.5d0,rhoc,1,rhor(1,isup),1)
         call daxpy(dfftp%nnr,0.5d0,rhoc,1,rhor(1,isdw),1)
         !
      endif
      !
      return
      !
endsubroutine add_cc_rspace
!
!
!
subroutine nksic_correction_pz( f, ispin, orb_rhor, &
                                vsic, pink, ibnd, shart)
      !
      ! ... calculate the non-Koopmans potential from the orbital density
      !
      use kinds,                only : dp
      use constants,            only : e2, fpi, hartree_si, electronvolt_si
      use cell_base,            only : tpiba2,omega
      use nksic,                only : nkscalfact, hartree_only_sic
      use gvect,                only : ngm, nl, nlm, gstart, gg 
      use eecp_mod,             only : do_comp
      use fft_interfaces,       only : fwfft, invfft
      use cp_interfaces,        only : fillgrad
      use fft_base,             only : dffts, dfftp
      use funct,                only : dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_bgrp_comm
      !  
      !use control_flags,        only : gamma_only, do_wf_cmplx
      !use control_flags,        only : hartree_only_sic
      !
      implicit none
      ! 
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(dfftp%nnr) 
      real(dp),    intent(out) :: vsic(dfftp%nnr)
      real(dp),    intent(out) :: pink, shart
      !
      ! internal variable
      !  
      integer       :: ig, nnrx
      real(dp)      :: etxc, ehele, fact
      !
      real(dp),    allocatable :: rhoelef(:,:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: vxc(:,:)
      ! 
      real(dp) :: dexc_dummy(3,3)
      !
      !==================
      ! main body
      !==================
      !
      nnrx = dfftp%nnr
      !
      vsic=0.0_dp
      pink=0.0_dp
      !
      ! this make sure that pz method does not work for
      ! empty states
      !   
      if ( f < 1.0d-6 ) return 
      !
      CALL start_clock( 'nk_corr' )
      CALL start_clock( 'nk_corr_h' )
      !
      fact=omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      allocate(rhogaux(ngm,2))
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      allocate(vhaux(nnrx))
      allocate(vxc(nnrx,2)) 
      !
      ! Compute self-hartree contributions
      !
      rhogaux=0.0_dp
      !
      ! rhoelef contains occupations
      !
      vhaux(:) = f*orb_rhor(:)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
         rhogaux(ig,ispin) = vhaux( nl(ig) )
      enddo
      !
      ! compute hartree-like potential
      !
      if ( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
         vtmp(ig) = rhogaux(ig,ispin) * fpi/( tpiba2*gg(ig) )
      enddo
      !
      ! compute periodic corrections
      !
      !if( do_comp ) then
          !
          !call calc_compensation_potential( vcorr, rhogaux(:,ispin), .true.)
          !vtmp(:) = vtmp(:) + vcorr(:)
          !
      !endif
      !
      vhaux=0.0_dp
      !  
      do ig=1,ngm
         !
         vhaux(nl(ig)) = vtmp(ig)
         vhaux(nlm(ig)) = CONJG(vtmp(ig))
         !
      enddo
      !
      call invfft('Dense',vhaux,dfftp)
      !
      ! init vsic
      !
      vsic(1:nnrx) =  -DBLE( vhaux(1:nnrx) )
      ehele        =   0.5_dp * sum( DBLE( vhaux(1:nnrx) ) &
                              * orb_rhor(1:nnrx) )
      !
      ! set ehele as measure of spread
      !
      shart=abs(ehele)*fact*hartree_si/electronvolt_si
      call mp_sum(shart, intra_bgrp_comm )
      !
      ! below is to make ehele quadratic in f (check this)
      ehele=ehele*f 
      !
      ! partial cleanup
      !
      deallocate( vtmp )
      deallocate( vcorr )
      deallocate( vhaux )
      !
      CALL stop_clock( 'nk_corr_h' )
      !
      ! Compute xc-contributions
      !
      if (.not.hartree_only_sic) then
         !
         if ( dft_is_gradient()) then
            !
            allocate(grhoraux(nnrx,3,2))
            !
            ! note: rhogaux contains the occupation
            !
            grhoraux=0.0_dp
            call fillgrad( 1, rhogaux(:,ispin:ispin), grhoraux(:,:,ispin:ispin))
            !
         else
            ! 
            allocate(grhoraux(1,1,1))
            !
            grhoraux=0.0_dp
            !
         endif
         !
         vxc=0.0_dp
         etxc=0.0_dp
         !
         vxc(:,ispin)=f*orb_rhor(:)
         !  
         call exch_corr_cp(nnrx, 2, grhoraux, vxc, etxc) 
         !
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc should be dexc(3,3), is lgam a variable here?
            ! 
            call gradh( 2, grhoraux, rhogaux, vxc, dexc_dummy)
            !
         endif
         !
         vsic(1:nnrx) =  vsic(1:nnrx) -vxc(1:nnrx,ispin)
         !
      else
         !
         etxc=0.0_dp
         !
      endif
      !
      ! energy correction terms
      !
      pink = fact * ( -etxc -ehele )
      !
      call mp_sum(pink, intra_bgrp_comm )
      !
      pink = pink * nkscalfact
      vsic = vsic * nkscalfact
      !
      deallocate( grhoraux )
      deallocate( rhogaux )
      deallocate( vxc )
      !
      call stop_clock( 'nk_corr' )
      !
      return
      !
end subroutine nksic_correction_pz
!
!
!
subroutine nksic_correction_nki( f, ispin, orb_rhor, rhor, &
                         rhoref, rhobar, rhobarg, grhobar, &
                         vsic, wxdsic, do_wxd_, pink, ibnd, shart, vxc_aux, etxc_aux, is_empty )
      !
      ! ... calculate the non-Koopmans (integrated, NKI)
      !     potential from the orbital density
      !
      !     note that fref=1.0 when performing NKI (i.e. it has a diff
      !     meaning)
      !     then  rho_ref = rho - rho_i + n_i
      !           rho_bar = rho - rho_i
      !
      use kinds,                only : dp
      use constants,            only : e2, fpi, hartree_si, electronvolt_si
      use cell_base,            only : tpiba2, omega
      use nksic,                only : fref, rhobarfact, nkscalfact 
      use fft_interfaces,       only : fwfft, invfft
      use fft_base,             only : dffts, dfftp
      use cp_interfaces,        only : fillgrad
      use electrons_base,       only : nspin
      use gvect,                only : ngm, nl, nlm, gstart, gg
      use eecp_mod,             only : do_comp
      use funct,                only : dmxc_spin, dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_bgrp_comm
      !
      implicit none
      ! 
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(dfftp%nnr)
      real(dp),    intent(in)  :: rhor(dfftp%nnr,nspin)
      real(dp),    intent(in)  :: rhoref(dfftp%nnr,2)
      real(dp),    intent(in)  :: rhobar(dfftp%nnr,2)
      complex(dp), intent(in)  :: rhobarg(ngm,2)
      real(dp),    intent(in)  :: grhobar(dfftp%nnr,3,2)
      logical,     intent(in)  :: do_wxd_
      real(dp),    intent(out) :: pink, shart
      real(dp),    intent(out) :: vxc_aux(dfftp%nnr,2)
      real(dp),    intent(out) :: etxc_aux
      real(dp),    intent(out) :: wxdsic(dfftp%nnr,2)
      real(dp),    intent(out) :: vsic(dfftp%nnr)
      logical, optional, intent(in) :: is_empty
      !
      ! local variables
      !
      integer       :: ig, nnrx
      real(dp)      :: fact, ehele, etmp
      real(dp)      :: etxcref, etxc0, w2cst
      !
      real(dp),    allocatable :: rhoele(:,:)
      real(dp),    allocatable :: rhoraux(:,:)
      real(dp),    allocatable :: vxc0(:,:)
      real(dp),    allocatable :: vxcref(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp),    allocatable :: orb_grhor(:,:,:)
      complex(dp), allocatable :: orb_rhog(:,:)
      ! 
      logical :: is_empty_
      real(dp):: icoeff
      real(dp):: dexc_dummy(3,3)
      !
      !==================
      ! main body
      !==================
      !
      nnrx = dfftp%nnr
      ! 
      icoeff=2.d0
      !
      IF(present(is_empty)) THEN
         !
         is_empty_=is_empty
         !
      ELSE
         !
         is_empty_=.false.
         !
      ENDIF
      ! 
      CALL start_clock( 'nki_corr' )
      !
      fact=omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3) 
      !
      allocate(rhoele(nnrx,2))
      allocate(rhogaux(ngm,2))
      allocate(vtmp(ngm))
      allocate(orb_rhog(ngm,1))
      allocate(vcorr(ngm))
      allocate(vhaux(nnrx))
      !
      rhoele=0.0d0
      rhoele(:,ispin) = orb_rhor(:)
      !
      vsic=0.0_dp
      wxdsic=0.0_dp
      pink=0.0_dp
      !
      ! Compute self-hartree contributions
      !
      orb_rhog=0.0_dp
      !
      ! rhoele has no occupation
      !
      vhaux(:) = rhoele(:,ispin)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
         !
         orb_rhog(ig,1) = vhaux( nl(ig) )
         !
      enddo
      !
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      !
      do ig=gstart,ngm
         !  
         vtmp(ig) = orb_rhog(ig,1) * fpi/( tpiba2*gg(ig) )
         !  
      enddo
      !
      ! compute periodic corrections
      !
      !if( do_comp ) then
          !
      !    call calc_compensation_potential( vcorr, orb_rhog(:,1),.true.)
      !    vtmp(:) = vtmp(:) + vcorr(:)
          !
      !endif
      !
      vhaux=0.0_dp
      do ig=1,ngm
         !
         vhaux(nl(ig)) = vtmp(ig)
         vhaux(nlm(ig)) = CONJG(vtmp(ig))
         !  
      enddo
      !
      call invfft('Dense',vhaux,dfftp)
      !
      ! init here vsic to save some memory
      !
      ! this is just the self-hartree potential
      !
      vsic(1:nnrx) = (1.0_dp-f) * DBLE( vhaux(1:nnrx) )
      !
      ! self-hartree contrib to pink
      ! and w2cst for vsic
      !
      ehele = icoeff * DBLE ( DOT_PRODUCT( vtmp(1:ngm), orb_rhog(1:ngm,1)))
      !
      if ( gstart == 2 ) then
         !
         ehele = ehele +(1.d0-icoeff)*DBLE ( CONJG( vtmp(1) )* orb_rhog(1,1) )
         !
      endif
      !
      shart=abs(ehele)*omega*0.5d0*hartree_si/electronvolt_si
      !
      call mp_sum(shart, intra_bgrp_comm)
      !
      ! self-hartree energy to be added to the vsic potential
      ! the scalar Hatree term of both empty and occupied states is 
      ! in the same form: -E_H[n_i]
      !
      w2cst = 0.0_dp
      !
      w2cst = -0.5_dp * ehele * omega
      !
      call mp_sum(w2cst, intra_bgrp_comm)
      !
      vsic  = vsic + w2cst
      !
      ! the f * (1-f) term is added here
      !
      IF(.not.is_empty_) THEN
         !
         ehele = 0.5_dp * f * (1.0_dp-f) * ehele * omega / fact
         !
      ELSE !this is for the fake functional for empty states
         !
         ehele = 0.5_dp * ehele * omega / fact
         !
      ENDIF
      !
      deallocate(vtmp)
      deallocate(vcorr)
      deallocate(vhaux)
      !
      !   add self-xc contributions
      !
      if ( dft_is_gradient() ) then
         !
         allocate(grhoraux(nnrx,3,2))
         allocate(orb_grhor(nnrx,3,1))
         !
         ! compute the gradient of n_i(r)
         call fillgrad( 1, orb_rhog, orb_grhor(:,:,1:1))
         !
      else
         !
         allocate(grhoraux(1,1,1))
         grhoraux=0.0_dp
         !
      endif
      !
      allocate(vxc0(nnrx,2))
      allocate(vxcref(nnrx,2))
      !
      ! this term is computed for ibnd, ispin == 1 and stored
      ! or if rhobarfact < 1
      !
      if ( ( ibnd == 1 .and. ispin == 1) .OR. rhobarfact < 1.0_dp ) then 
         !
         etxc_aux=0.0_dp
         vxc_aux=0.0_dp
         !
         ! some memory can be same in the nspin-2 case,
         ! considering that rhobar + f*rhoele is identical to rho
         ! when rhobarfact == 1
         !
         if ( dft_is_gradient() ) then
            !
            grhoraux(:,:,1:2)   = grhobar(:,:,1:2)
            grhoraux(:,:,ispin) = grhobar(:,:,ispin) &
                                + f * orb_grhor(:,:,1)
            !
            rhogaux(:,1:2) = rhobarg(:,1:2)
            rhogaux(:,ispin) = rhobarg(:,ispin) + f * orb_rhog(:,1)
            !
         endif
         !
         allocate( rhoraux(nnrx, 2) )
         !
         rhoraux = rhobar + f*rhoele
         vxc_aux=rhoraux
         !
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxc_aux, etxc_aux)
         !
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc
            !  should be dexc(3,3), is lgam a variable here?
            !     
            call gradh( 2, grhoraux, rhogaux, vxc_aux, dexc_dummy)
            !
         endif
         !
         deallocate( rhoraux )
         !
      endif
      !
      etxcref=0.0_dp
      vxcref=0.0_dp
      !
      if ( f == 1.0_dp ) then
         !
         vxcref=vxc_aux
         etxcref=etxc_aux
         !
      else
         !
         if ( dft_is_gradient() ) then
            !
            grhoraux(:,:,1:2)   = grhobar(:,:,1:2)
            grhoraux(:,:,ispin) = grhobar(:,:,ispin) &
                                + fref * orb_grhor(:,:,1)
            !
            rhogaux(:,1:2) = rhobarg(:,1:2)
            rhogaux(:,ispin) = rhobarg(:,ispin) + fref * orb_rhog(:,1)
            !
         endif
         !
         vxcref=rhoref
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxcref, etxcref)
         !
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc
            !  should be dexc(3,3), is lgam a variable here?
            ! 
            call gradh(2, grhoraux, rhogaux, vxcref, dexc_dummy)
            !
         endif
         !
      endif
      !
      !rhoraux = rhobar
      !
      etxc0=0.0_dp
      vxc0=0.0_dp
      !
      vxc0=rhobar
      CALL exch_corr_cp(nnrx, 2, grhobar, vxc0, etxc0)
      !
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !  Need a dummy dexc here, need to cross-check gradh! dexc
         !  should be dexc(3,3), is lgam a variable here?
         !
         call gradh(2, grhobar, rhobarg, vxc0, dexc_dummy)
         ! 
      endif
      !
      ! update potential (including other constant terms)
      ! and define pink
      !
      if (.not.is_empty_) then
         !
         etmp  = sum( vxcref(1:nnrx,ispin) * rhoele(1:nnrx,ispin) )
         w2cst = ( etxcref-etxc0 ) - etmp
         w2cst = w2cst * fact
         !
         call mp_sum(w2cst, intra_bgrp_comm)
         !
         pink = (1.0_dp - f) * etxc0 - etxc_aux + f * etxcref + ehele
         !
      else
         !
         etmp  = sum( vxcref(1:nnrx,ispin) * rhoele(1:nnrx,ispin) )
         w2cst = ( etxcref - etxc0 ) -etmp
         w2cst = w2cst * fact
         !
         call mp_sum(w2cst, intra_bgrp_comm)
         !
         etmp  = sum( vxc_aux(1:nnrx,ispin) * rhoele(1:nnrx,ispin) )
         !
         pink = etxcref - etxc0 - etmp + ehele
         !
      endif
      !
      pink = pink*fact
      !
      call mp_sum(pink, intra_bgrp_comm)
      !
      vsic(1:nnrx) = vsic(1:nnrx) &
                   + vxcref(1:nnrx,ispin) - vxc_aux(1:nnrx,ispin) + w2cst
      !
      !   calculate wxd
      !
      wxdsic(:,:) = 0.0d0
      !
      if ( do_wxd_ ) then
         !
         wxdsic(:,1:2)= (1.0_dp-f)*vxc0(:,1:2) - vxc_aux(:,1:2) + f*vxcref(:,1:2)
         !
      endif
      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      pink = pink * nkscalfact
      vsic = vsic * nkscalfact
      !
      if ( do_wxd_ ) then
         !
         wxdsic = wxdsic * nkscalfact
         !
      else
         !
         wxdsic = 0.d0
         !
      endif
      !
      deallocate(vxc0)
      deallocate(vxcref)
      deallocate(rhoele)
      !
      deallocate(grhoraux)
      deallocate(rhogaux)
      !
      if(allocated(orb_grhor)) deallocate(orb_grhor)
      !
      call stop_clock( 'nki_corr' )
      ! 
      return
      !
endsubroutine nksic_correction_nki
!
!
!
subroutine nksic_correction_nkipz( f, ispin, orb_rhor, &
                                   vsic, pink, ibnd, shart, is_empty)
      !
      ! calculate the non-Koopmans potential from the orbital density
      !
      use kinds,                only : dp
      use constants,            only : e2, fpi, hartree_si, electronvolt_si
      use cell_base,            only : tpiba2, omega
      use nksic,                only : nkscalfact
      use fft_interfaces,       only : fwfft, invfft
      use fft_base,             only : dffts, dfftp
      use cp_interfaces,        only : fillgrad
      use electrons_base,       only : nspin
      use gvect,                only : ngm, nl, nlm, gstart, gg
      use eecp_mod,             only : do_comp
      use funct,                only : dft_is_gradient
      use mp,                   only : mp_sum
      use mp_global,            only : intra_bgrp_comm
      !
      implicit none
      !
      integer,     intent(in)  :: ispin, ibnd
      real(dp),    intent(in)  :: f, orb_rhor(dfftp%nnr)
      real(dp),    intent(out) :: vsic(dfftp%nnr)
      real(dp),    intent(out) :: pink, shart
      logical, optional, intent(in) :: is_empty
      !
      !character(19) :: subname='nksic_correction_pz'
      !
      integer       :: ig, nnrx
      real(dp)      :: ehele, fact, w2cst, etmp, etxc_
      !
      real(dp),    allocatable :: rhoele(:,:)
      real(dp),    allocatable :: vxc_(:,:)
      complex(dp), allocatable :: rhogaux(:,:)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: vtmp(:)
      real(dp),    allocatable :: grhoraux(:,:,:)
      real(dp) :: icoeff
      real(dp) :: dexc_dummy(3,3)
      logical  :: is_empty_
      !
      !==================
      ! main body
      !==================
      !
      nnrx = dfftp%nnr
      !
      icoeff=2.d0
      !
      IF(present(is_empty)) THEN
         !
         is_empty_ = is_empty
         !
      ELSE
         !
         is_empty_ = .false.
         !
      ENDIF
      !
      vsic=0.0_dp
      pink=0.0_dp
      !
      CALL start_clock( 'nkipz_corr' )
      !
      fact=omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
      !
      allocate(rhogaux(ngm,2))
      allocate(vtmp(ngm))
      allocate(vcorr(ngm))
      allocate(vxc_(nnrx,2))
      allocate(vhaux(nnrx))
      !
      ! Compute self-hartree contributions
      !
      rhogaux=0.0_dp
      !
      ! vhaux does not contain occupations
      !
      vhaux(:) = orb_rhor(:)
      !
      call fwfft('Dense',vhaux,dfftp )
      !
      do ig=1,ngm
          rhogaux(ig,ispin) = vhaux( nl(ig) )
      enddo
      !
      ! compute hartree-like potential
      !
      if( gstart == 2 ) vtmp(1)=(0.d0,0.d0)
      do ig=gstart,ngm
          vtmp(ig) = rhogaux(ig,ispin) * fpi/( tpiba2*gg(ig) )
      enddo
      !
      ! compute periodic corrections
      !
      !if( do_comp ) then
          !
      !    call calc_compensation_potential( vcorr, rhogaux(:,ispin),.true. )
      !    vtmp(:) = vtmp(:) + vcorr(:)
          !
      !endif
      !
      vhaux=0.0_dp
      do ig=1,ngm
         !
         vhaux(nl(ig)) = vtmp(ig)
         vhaux(nlm(ig)) = CONJG(vtmp(ig))
         !  
      enddo
      !
      call invfft('Dense',vhaux,dfftp)
      !
      ! init vsic
      !
      vsic(1:nnrx) =  -DBLE( vhaux(1:nnrx) )
      !
      ehele = icoeff * DBLE ( DOT_PRODUCT( vtmp(1:ngm), rhogaux(1:ngm,ispin)))
      if ( gstart == 2 ) ehele = ehele + (1.d0-icoeff)*DBLE ( CONJG( vtmp(1) ) * rhogaux(1,ispin) )
      !
      w2cst = 0.0_dp
      ! 
      w2cst = 0.5_dp * ehele * omega
      !
      call mp_sum(w2cst, intra_bgrp_comm)
      !
      vsic  = vsic + w2cst
      !
      ehele = 0.5d0 * ehele * omega / fact
      !
      shart=abs(ehele)*fact*hartree_si/electronvolt_si
      !
      call mp_sum(shart, intra_bgrp_comm)
      !
      ! partial cleanup
      !
      deallocate( vtmp )
      deallocate( vcorr )
      deallocate( vhaux )
      !
      ! Compute xc-contributions
      !
      if ( dft_is_gradient() ) then
         ! 
         allocate(grhoraux(nnrx,3,2))
         !
         ! note: rhogaux does not contain the occupation
         !
         grhoraux=0.0_dp
         call fillgrad( 1, rhogaux(:,ispin:ispin), grhoraux(:,:,ispin:ispin) )
         ! 
      else
         !
         allocate(grhoraux(1,1,1))
         !
         grhoraux=0.0_dp
         !
      endif
      !
      vxc_=0.0_dp
      etxc_=0.0_dp
      !
      vxc_(:,ispin)=orb_rhor(:)
      call exch_corr_cp(nnrx, 2, grhoraux, vxc_, etxc_)
      !
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !  Need a dummy dexc here, need to cross-check gradh! dexc
         !  should be dexc(3,3), is lgam a variable here?
         !
         call gradh( 2, grhoraux, rhogaux, vxc_, dexc_dummy)
         !
      end if
      !
      if (.not.is_empty_) then
         !
         etmp  = sum( vxc_(1:nnrx,ispin) * orb_rhor(1:nnrx) )
         !
         w2cst = -etxc_ + etmp
         w2cst = w2cst * fact
         !
         call mp_sum(w2cst, intra_bgrp_comm)
         !
         pink = -f*(etxc_ + ehele)
         !
      else
         !
         etmp  = sum( vxc_(1:nnrx,ispin) * orb_rhor(1:nnrx) )
         !
         w2cst = -etxc_ + etmp
         w2cst =  w2cst * fact
         !
         call mp_sum(w2cst, intra_bgrp_comm)
         !
         pink = -(etxc_ + ehele)
         !
      endif 
      !
      pink = pink*fact
      !
      call mp_sum(pink, intra_bgrp_comm)
      !
      vsic(1:nnrx) =  vsic(1:nnrx) - vxc_(1:nnrx,ispin) + w2cst
      !
      pink = pink * nkscalfact
      vsic = vsic * nkscalfact
      !
      deallocate( grhoraux )
      deallocate( rhogaux )
      deallocate( vxc_ )
      !
      call stop_clock( 'nkipz_corr' )
      !
      return
      !
endsubroutine nksic_correction_nkipz
!
!
!
subroutine nksic_eforce( i, nbsp, nx, vsic, deeq_sic, bec, ngw, c1, c2, vsicpsi)
      !
      ! Compute vsic potential for orbitals i and i+1 (c1 and c2)
      !
      use kinds,                    only : dp
      use fft_interfaces,           only : fwfft, invfft
      use fft_base,                 only : dffts, dfftp
      use gvecs
      use gvect,                    only : ngm, nl, nlm                    
      use uspp,                     only : nkb, vkb
      use uspp_param,               only : nhm, nh
      use ions_base,                only : nsp, na, nat
      !
      implicit none

      !
      ! input/output vars
      !
      integer,       intent(in)  :: i, nbsp, nx, ngw
      real(dp),      intent(in)  :: vsic(dffts%nnr,nx)
      real(dp),      intent(in)  :: deeq_sic(nhm,nhm,nat,nx)
      real(dp),      intent(in)  :: bec(nkb,nbsp) 
      complex(dp),   intent(in)  :: c1(ngw), c2(ngw)
      complex(dp),   intent(out) :: vsicpsi(ngw, 2)
      !
      ! local vars
      !
      character(12) :: subname='nksic_eforce'
      integer       :: ir, ig, ierr, j
      integer       :: is, iv, jv, isa, ism, nnrx, nnrsx 
      integer       :: ivoff, jvoff, ia, inl, jnl
      real(dp)      :: wfc(2), dd
      complex(dp)   :: wfc_c(2)
      complex(dp)   :: fm, fp
      complex(dp),  allocatable :: psi1(:), psi2(:)
      real(dp),     allocatable :: aa(:,:)
      complex(dp),     allocatable :: aa_c(:,:)
      complex(dp), parameter :: c_one= CMPLX(1.d0,0.d0)
      !
      !====================
      ! main body
      !====================
      !
      call start_clock( 'nk_eforce' )
      !
      nnrx = dffts%nnr
      nnrsx = dfftp%nnr 
      !
      allocate( psi1(nnrx), stat=ierr )
      if ( ierr/=0 ) call errore(subname,'allocating psi1',abs(ierr))
      !
      vsicpsi(:,:) = (0.0_dp, 0.0_dp)
      !
      ! take advantage of the smooth and the dense grids
      ! being equal (NCPP case)
      !
      if ( nnrsx == nnrx ) then !waring:giovanni we are not using ultrasoft
         !
         ! no need to take care of the double grid.
         ! typically, NCPP case
         !
         CALL c2psi( psi1, nnrx, c1, c2, ngw, 2 ) !warning:giovanni need to change this
         !
         CALL invfft('Dense', psi1, dfftp )
         !
         ! computing the orbital wfcs
         ! and the potentials in real space on the full grid
         !
         do ir = 1, nnrx
            !
            wfc(1)    =  DBLE( psi1(ir) )
            wfc(2)    = AIMAG( psi1(ir) )
            !
            psi1( ir ) = CMPLX( wfc(1) * vsic(ir,i), &
                                wfc(2) * vsic(ir,i+1))
            !
         enddo
         !
         CALL fwfft('Dense', psi1, dfftp )
         !
         vsicpsi(:,:) = (0.0_dp, 0.0_dp)
         !
         do ig=1,ngw
            !
            fp = psi1(nls(ig))+psi1(nlsm(ig))
            fm = psi1(nls(ig))-psi1(nlsm(ig))
            !
            vsicpsi(ig,1)=0.5d0*CMPLX(DBLE(fp), AIMAG(fm), kind=dp)
            vsicpsi(ig,2)=0.5d0*CMPLX(AIMAG(fp),-DBLE(fm), kind=dp)
            !
         enddo
         !
      else
         !
         call errore(subname,'does not work with uspp', 0)
         !
         ! here we take properly into account the
         ! smooth and the dense grids
         ! typically, USPP case
         !
         !CALL nksic_eforce_std(lgam) !warning:giovanni this makes fourier transforms
         !
      endif
      !
      deallocate( psi1 )
      !
      ! add USPP non-local contribution
      ! to the potantial
      ! (this comes from the orbital-dependent piece of
      ! the potential)
      !
      if( nkb > 0 ) then
         call errore(subname,'does not work with uspp', 0)
      endif
      !
      call stop_clock( 'nk_eforce' )
      !
      return
      !
end subroutine nksic_eforce
