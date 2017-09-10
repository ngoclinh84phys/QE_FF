!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Electrostatic embedding methods 
! Developed and implemented by I. Dabo (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Parallelized by Andrea Ferretti (MIT)
!
!-----------------------------------------------------------------------
subroutine ee_green_0d_init(box)
!-----------------------------------------------------------------------
!
! ... initialize Green's functions for periodic-image correction
! ... in 0D setttings (e.g., isolated molecule, cluster)
!
      use kinds,              only : dp
      use cell_base,          only : at, omega, tpiba2, s_to_r, &
                                     boxdimensions, cell_alat
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use gvect,              only : gstart, gg
      use gvect,              only : ngm, nl, nlm
      use fft_base,           only : dfftp
      use fft_interfaces,     only : fwfft, invfft
      use eecp_mod,           only : gcorr,gcorr_fft
      use mp_global,          only : me_bgrp
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      real(dp),      parameter :: sigma=2.0_dp
      real(dp),      parameter :: vanishing_dist=1.0e-3_dp
      !
      complex(dp), allocatable :: vtemp(:)
      real(dp),    allocatable :: vtempr(:)
      real(dp) :: aux(dfftp%nr1,dfftp%nr2,dfftp%nr3)
      !
      integer             :: ig, ir, i, j, k
      integer             :: index0, ir_end, index
      integer             :: nnrx, nr1,  nr2,  nr3, nr1x, nr2x, &
                             nr3x
      real(dp)            :: sv(3), lv(3) ,dist
      real(dp)            :: a(3,3), alat_
      integer             :: npt(3)
      logical             :: tperiodic(3)
      real(dp),  external :: qe_erf
      !
      interface 
        !
        function afc(a,npt,tperiodic,spreadopt)
          !
          real(8), intent(in), optional :: spreadopt
          real(8), intent(in), dimension(3,3) :: a
          integer, intent(in), dimension(3) :: npt
          logical, intent(in), dimension(3) :: tperiodic
          real(8) :: afc(npt(1),npt(2),npt(3))
          !
       end function
       !
      end interface
      !
      ! main body
      !
      nr1=dfftp%nr1
      nr2=dfftp%nr2
      nr3=dfftp%nr3
      nr1x=dfftp%nr1x
      nr2x=dfftp%nr2x
      nr3x=dfftp%nr3x
      nnrx=dfftp%nnr
      !
      allocate(vtemp(nnrx))
      allocate(vtempr(nnrx))
      !
      vtemp=0.0_dp
      vtempr=0.0_dp
      gcorr=0.0_dp     
      ! 
      npt(1)=dfftp%nr1
      npt(2)=dfftp%nr2
      npt(3)=dfftp%nr3
      tperiodic=.false.
      !
      alat_=cell_alat()
      a(1:3,1)=at(1:3,1)*alat_
      a(1:3,2)=at(1:3,2)*alat_
      a(1:3,3)=at(1:3,3)*alat_
      !
      aux=afc(a,npt,tperiodic,sigma)
      !
      index0=0
      !careful, below me_bgrp starts from zero; the number of Z planes for me is dfftp%npp(me_bgrp+1)

#if defined (__MPI)
       DO i = 1, me_bgrp
         index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
       END DO
#endif
     !
#if defined (__MPI)
       ir_end = MIN(nnrx,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
       ir_end = nnrx
#endif
       !
       DO ir = 1, ir_end 
          !
          ! ... three dimensional indexes
          !
          index = index0 + ir - 1
          k     = index / (dfftp%nr1x*dfftp%nr2x)
          index = index - (dfftp%nr1x*dfftp%nr2x)*k
          j     = index / dfftp%nr1x
          index = index - dfftp%nr1x*j
          i     = index
          !
          gcorr(ir)=aux(i,j,k)
          !
       ENDDO
       !
!       call writetofile(gcorr,nnrx,'afc0d.dat',dfftp,'az')
       vtemp(:)=gcorr(:)
       call fwfft('Dense',vtemp,dfftp)
       do ig=1,ngm
          gcorr_fft(ig)=vtemp(nl(ig))
       enddo
      !
      deallocate(vtempr)
      deallocate(vtemp)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine ee_green_0d_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine ee_green_1d_init(box)
!-----------------------------------------------------------------------
!
! ... initialize Green's functions for periodic-image correction
! ... for 1D setttings (e.g., nanotube, polymer chain)
!
      use kinds,              only : dp
      use cell_base,          only : at, omega, tpiba2, s_to_r, & 
                                     boxdimensions, cell_alat
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use gvect,              only : gstart, gg, g
      use gvect,              only : ngm, nl, nlm
      use fft_base,           only : dfftp
      use fft_interfaces,     only : fwfft, invfft
      use eecp_mod,           only : gcorr1d,gcorr1d_fft
      use mp_global,          only : me_bgrp
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      real(dp),      parameter :: sigma=2.0_dp
      real(dp),      parameter :: vanishing_dist=1.0e-3_dp
      real(dp),      parameter :: vanishing_g=1.0e-3_dp
      real(dp),      parameter :: euler_gamma=0.57721566490153286061d0
      !
      complex(dp), allocatable :: vtemp(:)
      real(dp),    allocatable :: vtempr(:)
      real(dp) :: aux(dfftp%nr1,dfftp%nr2,dfftp%nr3)
      !
      integer             :: ig, ir, i, j, k
      integer             :: index0, ir_end, index
      real(dp)                 :: sv(3), lv(3), dist
      real(dp)            :: a(3,3), alat_
      integer             :: npt(3)
      logical             :: tperiodic(3)
      integer             :: nnrx, nr1, nr2, nr3, nr1x, nr2x, nr3x, &
                                           nr1l, nr2l, nr3l
      !
      real(dp),       external :: qe_erf
      real(dp),       external :: eimlmg
      !
      interface 
        !
        function afc(a,npt,tperiodic,spreadopt)
          !
          real(8), intent(in), optional :: spreadopt
          real(8), intent(in), dimension(3,3) :: a
          integer, intent(in), dimension(3) :: npt
          logical, intent(in), dimension(3) :: tperiodic
          real(8) :: afc(npt(1),npt(2),npt(3))
          !
        end function
        !
      end interface
      !
      ! main body
      !
      nr1=dfftp%nr1
      nr2=dfftp%nr2
      nr3=dfftp%nr3
      nr1x=dfftp%nr1x
      nr2x=dfftp%nr2x
      nr3x=dfftp%nr3x
      nnrx=dfftp%nnr

      allocate(vtemp(nnrx))
      allocate(vtempr(nnrx))
      ! 
      vtemp=0.0_dp
      vtempr=0.0_dp
      gcorr1d=0.0_dp
      !
      nr1=dfftp%nr1
      nr2=dfftp%nr2
      nr3=dfftp%nr3
      nr1x=dfftp%nr1x
      nr2x=dfftp%nr2x
      nr3x=dfftp%nr3x
      nnrx=dfftp%nnr
      !      
      npt(1)=dfftp%nr1
      npt(2)=dfftp%nr2
      npt(3)=dfftp%nr3
      tperiodic(1)=.false.
      tperiodic(2)=.false.
      tperiodic(3)=.true.
      !
      alat_=cell_alat()
      a(1:3,1)=at(1:3,1)*alat_
      a(1:3,2)=at(1:3,2)*alat_
      a(1:3,3)=at(1:3,3)*alat_
      !
      aux=afc(a,npt,tperiodic,sigma)
      !
      index0=0
      !careful, below me_bgrp starts from zero; the number of Z planes for me is dfftp%npp(me_bgrp+1)
#if defined (__MPI)
      DO i = 1, me_bgrp
         index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
#endif
      !
#if defined (__MPI)
      ir_end = MIN(nnrx,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
      ir_end = nnrx
#endif
!
      DO ir = 1, ir_end 
        !
        ! ... three dimensional indexes
        !
        index = index0 + ir - 1
        k     = index / (dfftp%nr1x*dfftp%nr2x)
        index = index - (dfftp%nr1x*dfftp%nr2x)*k
        j     = index / dfftp%nr1x
        index = index - dfftp%nr1x*j
        i     = index
        !
        gcorr1d(ir)=aux(i,j,k)
        !
     ENDDO

!     call writetofile(gcorr1d,nnrx,'afc1d.dat',dfftp, 'ax')
     !
     vtemp=gcorr1d
     call fwfft('Dense',vtemp,dfftp)
     do ig=1,ngm
       gcorr1d_fft(ig)=vtemp(nl(ig))
     enddo
      !
      deallocate(vtempr)
      deallocate(vtemp)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine ee_green_1d_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine ee_green_2d_init(box)
!-----------------------------------------------------------------------
!
! ... initialize Green's functions for periodic-image correction
! ... for 2D setttings (e.g., surface, thin film)
!
      use kinds,              only : dp
      use cell_base,          only : at, omega, tpiba2, s_to_r, &
                                  boxdimensions, cell_alat
      use constants,          only : fpi, pi
      use io_global,          only : stdout
      use gvect,              only : gstart, gg, g
      use gvect,              only : ngm, nl, nlm
      use fft_base,           only : dfftp
      use fft_interfaces,     only : fwfft, invfft
      use eecp_mod,           only : gcorr2d,gcorr2d_fft
      use mp_global,          only : me_bgrp
      !
      implicit none
      !
      type(boxdimensions), intent(in) :: box
      !
      real(dp),      parameter :: sigma=2.0_dp
      real(dp),      parameter :: vanishing_dist=1.0e-3_dp
      real(dp),      parameter :: vanishing_g=1.0e-3_dp
      real(dp),      parameter :: euler_gamma=0.57721566490153286061d0
      !
      complex(dp), allocatable :: vtemp(:)
      real(dp),    allocatable :: vtempr(:)
      real(dp) :: aux(dfftp%nr1,dfftp%nr2,dfftp%nr3)
      !
      integer             :: ig, ir, i, j, k
      integer             :: index0, ir_end, index
      integer             :: nnrx, nr1,  nr2,  nr3, nr1x, nr2x, &
                           nr3x
      real(dp)                 :: sv(3), lv(3), dist
      real(dp)            :: a(3,3), alat_
      integer             :: npt(3)
      logical             :: tperiodic(3)
      !
      real(dp),       external :: qe_erf
      real(dp),       external :: eimlmg
      !
      interface 
      !
      function afc(a,npt,tperiodic,spreadopt)
        !
        real(8), intent(in), optional :: spreadopt
        real(8), intent(in), dimension(3,3) :: a
        integer, intent(in), dimension(3) :: npt
        logical, intent(in), dimension(3) :: tperiodic
        real(8) :: afc(npt(1),npt(2),npt(3))
        !
      end function
      !
      end interface
      !
      ! main body
      !
      allocate(vtemp(nnrx))
      allocate(vtempr(nnrx))
      ! 
      vtemp=0.0_dp
      vtempr=0.0_dp
      gcorr2d=0.0_dp
      !
      nr1=dfftp%nr1
      nr2=dfftp%nr2
      nr3=dfftp%nr3
      nr1x=dfftp%nr1x
      nr2x=dfftp%nr2x
      nr3x=dfftp%nr3x
      nnrx=dfftp%nnr
      !   
      npt(1)=dfftp%nr1
      npt(2)=dfftp%nr2
      npt(3)=dfftp%nr3
      tperiodic(1)=.true.
      tperiodic(2)=.true.
      tperiodic(3)=.false.
      !
      alat_=cell_alat()
      a(1:3,1)=at(1:3,1)*alat_
      a(1:3,2)=at(1:3,2)*alat_
      a(1:3,3)=at(1:3,3)*alat_
      !
      aux=afc(a,npt,tperiodic,sigma)

      index0=0
      !careful, below me_bgrp starts from zero; the number of Z planes for me is dfftp%npp(me_bgrp+1)

#if defined (__MPI)
      DO i = 1, me_bgrp
      index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
      END DO
#endif
      !
#if defined (__MPI)
      ir_end = MIN(nnrx,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
      ir_end = nnrx
#endif
      !
      DO ir = 1, ir_end 
       !
       ! ... three dimensional indexes
       !
       index = index0 + ir - 1
       k     = index / (dfftp%nr1x*dfftp%nr2x)
       index = index - (dfftp%nr1x*dfftp%nr2x)*k
       j     = index / dfftp%nr1x
       index = index - dfftp%nr1x*j
       i     = index
       !
       gcorr2d(ir)=aux(i,j,k)
       !
      ENDDO
      !
!      call writetofile(gcorr2d,nnrx,'afc2d.dat',dfftp, 'ax')
      !
      vtemp=gcorr2d
      call fwfft('Dense',vtemp,dfftp)
      do ig=1,ngm
       gcorr2d_fft(ig)=vtemp(nl(ig))
      enddo
      !
      deallocate(vtempr)
      deallocate(vtemp)
      !
      return
      !
!-----------------------------------------------------------------------
   end subroutine ee_green_2d_init
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_compensation_potential(vcorr_fft,rho_fft,odd_flag)
!-----------------------------------------------------------------------
!
! ... driver for calculating the truncated countercharge (TCC) 
! ... for 0,1,2d periodicity
!
      use kinds,              only: dp
      use gvect,              only: ngm
      use eecp_mod,           only: gcorr_fft, which_compensation, tcc_odd
      use cell_base,          only: omega

      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      logical :: odd_flag !true if compensation is being computed
                          !for self-interaction correction

      select case(trim(which_compensation))
          !
        case('tcc')
          !
          call calc_tcc_potential(vcorr_fft,rho_fft)
          !
        case('tcc1d')
          !
          IF((.not.odd_flag).or.(.not.tcc_odd)) THEN
            call calc_tcc1d_potential(vcorr_fft,rho_fft)
          ELSE IF(odd_flag.and.tcc_odd) THEN
            call calc_tcc_potential(vcorr_fft,rho_fft)
          ENDIF
          !
        case('tcc2d')
          !
          IF((.not.odd_flag).or.(.not.tcc_odd)) THEN
            call calc_tcc2d_potential(vcorr_fft,rho_fft)
          ELSE IF(odd_flag.and.tcc_odd) THEN
            call calc_tcc_potential(vcorr_fft,rho_fft)
          ENDIF
          !
        case('none')
          !
          IF((.not.odd_flag).or.(.not.tcc_odd)) THEN
            continue
          ELSE IF(odd_flag.and.tcc_odd) THEN
            call calc_tcc_potential(vcorr_fft,rho_fft)
          ENDIF
          !
        case default
          !
          call errore('vofrho','Invalid correction: '//TRIM(which_compensation), 10)
          !
        end select
!-----------------------------------------------------------------------
      end subroutine calc_compensation_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc_potential(vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image correction potential in
! ... reciprocal space for 0D settings
!
      use kinds,              only: dp
      use gvect,              only: ngm
      use eecp_mod,           only: gcorr_fft
      use cell_base,          only: omega
      !
      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      integer :: ig
      !
      do ig=1,ngm
        vcorr_fft(ig)=omega*gcorr_fft(ig)*rho_fft(ig)
      end do
      !
      return
!
!-----------------------------------------------------------------------
      end subroutine calc_tcc_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc1d_potential(vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image correction potential in
! ... reciprocal space for 1D settings
!
      use kinds,              only: dp
      use gvect,              only: ngm
      use eecp_mod,           only: gcorr1d_fft
      use cell_base,          only: omega
      !
      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      integer :: ig
      !
      do ig=1,ngm
        vcorr_fft(ig)=omega*gcorr1d_fft(ig)*rho_fft(ig)
      end do
      !
      return
!
!-----------------------------------------------------------------------
      end subroutine calc_tcc1d_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc2d_potential(vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image correction potential in
! ... reciprocal space for 2D settings
!
      use kinds,              only: dp
      use gvect,              only: ngm
      use eecp_mod,           only: gcorr2d_fft
      use cell_base,          only: omega
      !
      implicit none
      complex(dp) :: rho_fft(ngm)
      complex(dp) :: vcorr_fft(ngm)
      integer :: ig
      !
      do ig=1,ngm
        vcorr_fft(ig)=omega*gcorr2d_fft(ig)*rho_fft(ig)
      end do
      !
      return
!
!-----------------------------------------------------------------------
      end subroutine calc_tcc2d_potential
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine calc_tcc_energy(ecomp,vcorr_fft,rho_fft)
!-----------------------------------------------------------------------
!
! ... calculate the truncated countercharge (TCC) 
! ... periodic-image corrective energy in for 0D settings
!
      use kinds,              only : dp
      use gvect,              only : ngm
      use eecp_mod,           only : gcorr_fft
      use cell_base,          only : omega
      use gvect,              only : gstart
      use mp,                 only : mp_sum
      use mp_global,          only : intra_bgrp_comm
      !
      implicit none
      !
      real(dp),    intent(out) :: ecomp
      complex(dp), intent(in)  :: rho_fft(ngm)
      complex(dp), intent(in)  :: vcorr_fft(ngm)
      !
      complex(dp), allocatable :: aux(:)
      integer      :: ig
      complex(dp)  :: zh
      real(dp), parameter :: wz=2.0_dp
      !
      allocate(aux(ngm))
      !
      aux=0.0_dp
      !
      if(gstart.ne.1) then
        aux(1)=0.5d0*omega*vcorr_fft(1)*conjg(rho_fft(1))
      end if
      !
      do ig=gstart,ngm
        aux(ig)=0.5d0*wz*omega*vcorr_fft(ig)*conjg(rho_fft(ig))
      end do
      !
      zh=0.0_dp
      do ig=1,ngm
        zh=zh+aux(ig)
      enddo
      ecomp=dble(zh)
      !
      call mp_sum(ecomp,intra_bgrp_comm)
      !
      deallocate(aux)
      !
      return
      !
!-----------------------------------------------------------------------
      end subroutine calc_tcc_energy
!-----------------------------------------------------------------------
