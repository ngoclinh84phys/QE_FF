module eecp_mod

  USE kinds

  implicit none
  real(dp), allocatable :: gcorr(:)
  complex(dp), allocatable :: gcorr_fft(:)
  real(dp), allocatable :: gcorr1d(:)
  complex(dp), allocatable :: gcorr1d_fft(:)
  real(dp), allocatable :: gcorr2d(:)
  complex(dp), allocatable :: gcorr2d_fft(:)
  real(dp), allocatable :: vcorr(:)
  complex(dp), allocatable :: vcorr_fft(:)
  character (len=256) :: which_compensation
  logical  :: do_comp
  real(dp) :: ecomp
  logical  :: tcc_odd 
  !
  contains
  !
  subroutine allocate_ee(nnrx,ngm)
  implicit none
  integer, intent(in):: nnrx
  integer, intent(in):: ngm
  !
  allocate(gcorr(nnrx))
  allocate(gcorr_fft(nnrx))
  allocate(gcorr1d(nnrx))
  allocate(gcorr1d_fft(nnrx))
  allocate(gcorr2d(nnrx))
  allocate(gcorr2d_fft(nnrx))
  allocate(vcorr(nnrx))
  allocate(vcorr_fft(ngm))
  !
  end subroutine
  !
  subroutine deallocate_ee
  !
  if(allocated(gcorr)) deallocate(gcorr)
  if(allocated(gcorr_fft)) deallocate(gcorr_fft)
  if(allocated(gcorr1d)) deallocate(gcorr1d)
  if(allocated(gcorr1d_fft)) deallocate(gcorr1d_fft)
  if(allocated(gcorr2d)) deallocate(gcorr2d)
  if(allocated(gcorr2d_fft)) deallocate(gcorr2d_fft)
  if(allocated(vcorr)) deallocate(vcorr)
  if(allocated(vcorr_fft)) deallocate(vcorr_fft)
  !
  end subroutine
  !
end module eecp_mod
!
subroutine ee_init
    !
    use eecp_mod,   only : do_comp, which_compensation, allocate_ee, tcc_odd 
    !use input_parameters, only : do_ee, &
    !                             which_compensation_ => which_compensation, &
    !                             tcc_odd_=> tcc_odd
    use io_global,         only : stdout
    use fft_base,          only : dfftp 
    use gvect,             only : ngm
    use cp_main_variables, only : ht0
    !
    implicit none
    !
    do_comp = .False.!do_ee
    which_compensation = 'tcc'!which_compensation_
    tcc_odd = .True.!tcc_odd_
    !
    if (do_comp) then 
       !  
       write(stdout, 2010) which_compensation
       !
       call allocate_ee(dfftp%nnr, ngm)
       !     
       write(stdout,*) "EE using tcc for odd: ", tcc_odd
       ! 
       if (trim(which_compensation)=='tcc1d') then
          !
          call ee_green_1d_init( ht0 )
          !  
          if (tcc_odd) then
             ! 
             call ee_green_0d_init( ht0 )
             ! 
          endif
          !  
       elseif (trim(which_compensation)=='tcc2d') then
          !   
          call ee_green_2d_init( ht0 )
          !
          if (tcc_odd) then
             ! 
             call ee_green_0d_init( ht0 )
             ! 
          endif
          !
       else
          ! 
          call ee_green_0d_init( ht0 )
          !
       endif
       !
    endif
    !
2010 format( 3X,'EE with periodic-image correction method = ',a20)
    !
    return
    !
end subroutine ee_init
