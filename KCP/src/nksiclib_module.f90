module nksic
  !
  use kinds
  !
  implicit none
  !
  save
  !
  real(dp) :: fref
  real(dp) :: rhobarfact
  real(dp) :: nkscalfact
  real(dp) :: kfact
  real(dp) :: epsi_cutoff_renorm=1.d-7
  real(dp) :: epsi2_cutoff_renorm=1.d-6
  real(dp) :: vanishing_rho_w
  real(dp) :: f_cutoff
  real(dp) :: eodd
  real(dp) :: etxc_sic
  !
  real(dp),allocatable :: fsic(:)
  real(dp),allocatable :: vsic(:,:)
  real(dp),allocatable :: fion_sic(:,:)
  real(dp),allocatable :: deeq_sic(:,:,:,:)
  real(dp),allocatable :: pink(:)
  real(dp),allocatable :: pink_emp(:)
  real(dp),allocatable :: odd_alpha(:)
  real(dp),allocatable :: wtot(:,:)
  complex(dp),allocatable :: vsicpsi(:,:)
  !
  !complex(dp) :: complexification_index
  !
  real(dp),    allocatable :: alpha0_ref(:)
  real(dp),    allocatable :: alpha0_ref_emp(:)
  complex(dp), allocatable :: swfc_fixed(:,:)
  complex(dp), allocatable :: becwfc_fixed(:,:)
  complex(dp), allocatable :: valpsi(:,:)
  !
  integer :: call_index = 0
  integer :: call_index_emp = 0
  !
  logical :: do_orbdep
  logical :: icompute_spread
  logical :: do_nk
  logical :: do_pz
  logical :: do_pzha
  logical :: do_pz_renorm
  logical :: do_bare_eigs
  logical :: do_nkpz
  logical :: do_nkipz
  logical :: do_nki
  logical :: do_spinsym
  logical :: do_wxd
  logical :: do_wref
  logical :: odd_nkscalfact
  logical :: hartree_only_sic
  !
  logical :: do_innerloop
  integer :: innerloop_cg_nsd
  integer :: innerloop_cg_nreset
  integer :: innerloop_nmax
  integer :: innerloop_init_n
  integer :: innerloop_atleast
  integer :: innerloop_until
  real(dp):: innerloop_cg_ratio
  real(dp):: esic_conv_thr 
  !
  integer :: sizwtot
  !
  ! For empty state calculations
  !
contains
  !
  subroutine allocate_nksic( nnrx, ngw, nspin, nx, nat)
     !
     use funct,            only : dft_is_gradient
     use uspp_param,       only : nhm
     !
     implicit none
     integer, intent(in):: nx,nspin
     integer, intent(in):: nat
     integer, intent(in):: ngw
     integer, intent(in):: nnrx
     !
     integer :: ispin
     logical :: lgam
     !
     allocate( fsic(nx) )
     allocate( vsic(nnrx,nx) )
     allocate( fion_sic(3,nat) )
     !
     allocate( pink(nx) )
     allocate( vsicpsi(ngw,2) )
     !
     if ( nhm > 0 ) then
        !
        allocate( deeq_sic(nhm, nhm, nat, nx) )
        !
     else
        !
        allocate( deeq_sic(1, 1, nat, nx) )
        ! 
     endif
     !
     if ( do_nki .or. do_nkipz) then
        ! 
        allocate(wtot(nnrx,2))
        !
        sizwtot=nnrx
        !
     else
        !
        allocate(wtot(1,2))
        ! 
        sizwtot=1
        !
     endif
     ! 
     wtot=0.0_dp
     !
     allocate( odd_alpha(nx) )
     !
     odd_alpha(:)= 0.d0
     !
     fsic     = 0.0d0
     pink     = 0.0d0
     vsic     = 0.0d0
     deeq_sic = 0.0d0
     vsicpsi  = 0.0d0
     !
  end subroutine allocate_nksic
  !
  real(dp) function nksic_memusage( )
     !
     ! output in MB (according to 4B integers and 8B reals)  
     real(dp) :: cost
     !
     cost = 0.0_dp
     if ( allocated(fsic) )       cost = cost + real( size(fsic) )       *  8.0_dp 
     if ( allocated(vsic) )       cost = cost + real( size(vsic) )       *  8.0_dp 
     if ( allocated(fion_sic) )   cost = cost + real( size(fion_sic) )   *  8.0_dp 
     if ( allocated(deeq_sic) )   cost = cost + real( size(deeq_sic) )   *  8.0_dp 
     if ( allocated(pink) )       cost = cost + real( size(pink) )       *  8.0_dp 
     if ( allocated(vsicpsi) )    cost = cost + real( size(vsicpsi) )    * 16.0_dp 
     if ( allocated(wtot) )       cost = cost + real( size(wtot) )       *  8.0_dp
     !
     nksic_memusage = cost / 1000000.0_dp
     !   
  end function nksic_memusage
  !
  subroutine deallocate_nksic
      !
      !use input_parameters, only: odd_nkscalfact
      if(allocated(fsic))        deallocate(fsic)
      if(allocated(vsic))        deallocate(vsic)
      if(allocated(fion_sic))    deallocate(fion_sic)
      if(allocated(deeq_sic))    deallocate(deeq_sic)
      if(allocated(pink))        deallocate(pink)
      if(allocated(odd_alpha))   deallocate(odd_alpha)
      if(allocated(vsicpsi))     deallocate(vsicpsi)
      if(allocated(wtot))        deallocate(wtot)
         !
      !if (odd_nkscalfact) then
         !
         !if (allocated(valpsi)) deallocate(valpsi)
         !if (allocated(alpha0_ref)) deallocate(alpha0_ref)
         !if (allocated(alpha0_ref_emp)) deallocate(alpha0_ref_emp)
         !if (allocated(swfc_fixed)) deallocate(swfc_fixed)
         !  
         !call deallocate_twin(becwfc_fixed)
         !
      !endif
      !
  end subroutine deallocate_nksic
  !
end module nksic
!
subroutine pc2nc_nksic(a, b, n, ispin)
  !      
  ! this function applies the modified Pc operator which is
  ! equivalent to Lowdin orthonormalization of the revised wavefunctions.
  ! currently implemented only for norm-conserving pseudopotentials. 
  !
  ! this subroutine applies the modified Pc operator
  ! a input :unperturbed wavefunctions
  ! b input :first order wavefunctions
  ! b output:b_i =b_i - |a_j>(<a_j|b_i>+<b_j|a_i>)/2
  !
  use kinds
  !use io_global, only: stdout
  use mp, only: mp_sum
  use mp_global, only: intra_bgrp_comm
  use gvecw, only: ngw
  use gvect, only: gstart
  !
  implicit none
  !
  integer, intent(in) :: n, ispin(n)
  complex(dp), intent(inout) :: a(ngw,n), b(ngw,n)
  ! 
  ! local variables
  !
  complex(dp) :: bold(ngw,n)
  integer i,j,ig
  complex(dp) :: sca_c
  real(dp), allocatable:: scar(:)
  complex(dp), allocatable:: scar_c(:)
  !
  allocate(scar_c(n))
  !
  bold(:,:)=b(:,:)
  !
  do j=1,n
     !
     do i=1,n
        !
        sca_c=CMPLX(0.0d0,0.d0)
        !
        if (ispin(i) == ispin(j)) then
           !
           if (gstart==2) bold(1,i) = CMPLX(DBLE(bold(1,i)),0.0d0)
           !
           do ig=1,ngw
              !
              sca_c=sca_c+CONJG(a(ig,j))*bold(ig,i) !uncomment this for lowdin ortho
              sca_c=sca_c+(a(ig,i))*CONJG(bold(ig,j)) !remove the 2.d0 for lowdin ortho
              !
           enddo
           !sca = sca*2.0d0  !2. for real weavefunctions
           !$$ not necessary: sca = sca*2.0d0  !2. for real weavefunctions
           if (gstart==2) then
              !
              sca_c = CMPLX(DBLE(sca_c),0.d0)-CMPLX(0.5d0*DBLE(CONJG(a(1,j)) & 
                      *(bold(1,i))+(a(1,i))*CONJG(bold(1,j))),0.d0) !usethis one for lowdin ortho
              !
           else
              !
              sca_c = CMPLX(DBLE(sca_c), 0.d0)
              !
           endif
           !
           scar_c(i) = sca_c
           !
        endif
        !  
     enddo
     !
     call mp_sum( scar_c, intra_bgrp_comm )
     !
     do i=1,n
        !
        if (ispin(i) == ispin(j)) then
           !
           sca_c = scar_c(i)
           !
           do ig=1,ngw
              !
              b(ig,i)=b(ig,i)-sca_c*a(ig,j)
              !
           enddo
           !
           if (gstart == 2) b(1,i) = CMPLX(DBLE(b(1,i)),0.0d0)
           !
        endif
        !
     enddo
     !
  enddo
  !
  deallocate(scar_c)
  !
  return
  !
end subroutine pc2nc_nksic
!
!
!
subroutine zdiag(nx,n,amat,dval,cvec,iflag)
  !
  use kinds , only : dp
  use zhpev_module, only: zhpev_drv
  !    
  implicit none
  !
  integer  ::  nx,n,iflag
  real(dp) ::  dval(n)
  complex(dp) :: amat(nx,n), cvec(nx,n)
  !
  ! internal vars
  !
  integer :: ndim,k,i,j
  complex(dp), allocatable :: ap(:)
  ! 
  ndim=(n*(n+1))/2
  allocate(ap(ndim))
  !
  ap(:)=CMPLX(0.d0,0.d0)
  k=0
  do j=1,n
     do i=1,j
        k=k+1
        ap(k)=amat(i,j)
     enddo
  enddo
  !
  call zhpev_drv( 'V', 'U', n, ap, dval, cvec, nx )
  !
  deallocate(ap)
  !
  return
  !
endsubroutine zdiag



