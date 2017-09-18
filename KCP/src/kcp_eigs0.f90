!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine kcp_eigs0( ei, nudx, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nlam, desc )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds,             only : DP
      use io_global,         only : stdout
      use constants,         only : autoev
      use dspev_module,      only : dspev_drv, pdspev_drv
      USE sic_module,        only : self_interaction
      USE descriptors,       ONLY : la_descriptor
      USE mp,                only : mp_sum, mp_bcast
      USE mp_global,         only : intra_bgrp_comm, root_bgrp, me_bgrp
      USE cg_module,         ONLY : tcg
      USE nksic,             ONLY : f_cutoff
      !
      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin), nlam
      type(la_descriptor), intent(in) :: desc( 2 )
      real(DP), intent(in) :: lambda( nlam, nlam, nspin ), f( nx )
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: ap(:), wr(:)
      real(DP) zr(1)
      integer :: iss, j, i, ierr, k, n, ndim, nspin_eig, npaired
      INTEGER :: ir, ic, nr, nc, nrl, nrlx, comm, np, me
      logical :: tsic
      CHARACTER(LEN=80) :: msg
!
      tsic = ( ABS( self_interaction) /= 0 )

      IF( tsic ) THEN
         nspin_eig = 1
         npaired   = nupdwn(2)
      ELSE
         nspin_eig = nspin
         npaired   = 0
      END IF

      do iss = 1, nspin_eig

         IF( nudx < nupdwn(iss) ) THEN 
            WRITE( msg, 100 ) nudx, SIZE( ei, 1 ), nupdwn(iss)
100         FORMAT( ' wrong dimension array ei = ', 3I10 )
            CALL errore( ' eigs0 ', msg, 1 )
         END IF

         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF

         allocate( wr( n ) )

         IF( desc( iss )%active_node > 0 ) THEN

            np = desc( iss )%npc * desc( iss )%npr

            IF( np > 1 ) THEN

               !  matrix is distributed

               CALL qe_pdsyevd( .false., n, desc(iss), lambda(1,1,iss), nlam, wr )

            ELSE

               !  matrix is not distributed

               allocate( ap( n * ( n + 1 ) / 2 ) )

               k = 0
               do i = 1, n
                  do j = i, n
                     k = k + 1
                     ap( k ) = lambda( j, i, iss )
                  end do
               end do

               CALL dspev_drv( 'N', 'L', n, ap, wr, zr, 1 )

               deallocate( ap )

            END IF

         END IF

         call mp_bcast( wr, root_bgrp, intra_bgrp_comm )

         if( lf ) then
            do i = 1, n
               !
               if (tcg) then
                  !
                  wr(i)=wr(i)/max(f(iupdwn(iss)-1+i),f_cutoff)
                  !
               else
                  !
                  if ( f(iupdwn(iss)-1+i).gt.1.e-6) then
                     ! 
                     wr(i)=wr(i)/f(iupdwn(iss)-1+i)
                     !
                  else
                     !
                     wr(i)=wr(i)/2.0d0 * nspin  ! fake occupation factor to print empty states
                     !
                  end if
                  !
               endif
               !
            end do
         end if
         !
         !     store eigenvalues
         !
         ei( 1:n, iss ) = wr( 1:n )

         IF( tsic ) THEN
            !
            !  store unpaired state
            !
            ei( 1:n,       1 ) = ei( 1:n, 1 ) / 2.0d0
            ei( nupdwn(1), 1 ) = 0.0d0
            if( desc( iss )%active_node > 0 ) then
               IF( desc( iss )%myc == desc( iss )%myr ) THEN
                  ir = desc( iss )%ir
                  nr = desc( iss )%nr
                  IF( nupdwn(1) >= ir .AND. nupdwn(1) < ir + nr ) then
                     ei( nupdwn(1), 1 ) = lambda( nupdwn(1)-ir+1, nupdwn(1)-ir+1, 1 )
                  end if
               END IF
            endif
            call mp_sum( ei( nupdwn(1), 1 ), intra_bgrp_comm )
            !
         END IF

         ! WRITE( stdout,*)  '---- DEBUG ----' ! debug
         ! WRITE( stdout,14) ( wr( i ) * autoev / 2.0d0, i = 1, nupdwn(iss) ) ! debug

         deallocate( wr )

      end do
      !
      !
      do iss = 1, nspin

         IF( tsic .AND. iss == 2 ) THEN
            ei( 1:npaired, 2 ) = ei( 1:npaired, 1 )
         END IF

         IF( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0.d0, 0.d0, 0.d0
            WRITE( stdout,14) ( ei( i, iss ) * autoev, i = 1, nupdwn(iss) )

         ENDIF

      end do

      IF( tprint ) WRITE( stdout,*)

   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
      !
      return
      ! 
end subroutine kcp_eigs0

!-----------------------------------------------------------------------
SUBROUTINE kcp_eigs_x( nfi, lambdap, lambda, descla )
!-----------------------------------------------------------------------

      USE kinds,             ONLY: DP
      use ensemble_dft,      only: tens
      use electrons_base,    only: nbspx, f, nspin
      use electrons_base,    only: iupdwn, nupdwn, nudx
      use electrons_module,  only: ei
      use io_global,         only: stdout
      USE descriptors,       ONLY: la_descriptor

      IMPLICIT NONE

      INTEGER :: nfi
      REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
      TYPE(la_descriptor), INTENT(IN) :: descla( : )

      if( .not. tens ) then
         call kcp_eigs0( ei, nudx, .false. , nspin, nupdwn, iupdwn, .true. , f, nbspx, lambda, SIZE(lambda,1), descla )
      else
         call kcp_eigs0( ei, nudx, .false. , nspin, nupdwn, iupdwn, .false. , f, nbspx, lambdap, SIZE(lambdap,1), descla )
      endif

      RETURN
END SUBROUTINE kcp_eigs_x
