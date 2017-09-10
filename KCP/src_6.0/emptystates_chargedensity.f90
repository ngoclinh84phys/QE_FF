!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE rhoofr_cp_generalized_x ( nbspx_bgrp, nbsp_bgrp, ispin_bgrp, f_bgrp, & 
       nfi, c_bgrp, irb, eigrb, bec_bgrp, dbec, rhovan, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin, tstress, ndwwf )
!-----------------------------------------------------------------------
!
!  this routine computes:
!  rhor  = normalized electron density in real space
!  ekin  = kinetic energy
!  dekin = kinetic energy term of QM stress
!
!    rhor(r) = (sum over ib) fi(ib) |psi(r,ib)|^2
!
!    Using quantities in scaled space
!    rhor(r) = rhor(s) / Omega
!    rhor(s) = (sum over ib) fi(ib) |psi(s,ib)|^2 
!
!    fi(ib) = occupation numbers
!    psi(r,ib) = psi(s,ib) / SQRT( Omega ) 
!    psi(s,ib) = INV_FFT (  c0(ig,ib)  )
!
!    ib = index of band
!    ig = index of G vector
!  ----------------------------------------------
!     the normalized electron density rhor in real space
!     the kinetic energy ekin
!     subroutine uses complex fft so it computes two ft's
!     simultaneously
!
!     rho_i,ij = sum_n < beta_i,i | psi_n >< psi_n | beta_i,j >
!     < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                   2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
!
!     e_v = sum_i,ij rho_i,ij d^ion_is,ji
!
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: iprint, iverbosity, thdyn, tpre, trhor
      USE ions_base,          ONLY: nat
      USE gvect,              ONLY: ngm,  nl, nlm
      USE gvecs,              ONLY: ngms, nls, nlsm
      USE smallbox_gvec,              ONLY: ngb
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE uspp,               ONLY: nkb
      USE uspp_param,         ONLY: nh, nhm
      USE cell_base,          ONLY: omega
      USE electrons_base,     ONLY: nspin
      USE constants,          ONLY: pi, fpi
      USE mp,                 ONLY: mp_sum
      USE io_global,          ONLY: stdout, ionode
      USE mp_global,          ONLY: intra_bgrp_comm, nbgrp, inter_bgrp_comm, me_bgrp, nproc_bgrp
      USE funct,              ONLY: dft_is_meta
      USE cg_module,          ONLY: tcg
      USE cp_interfaces,      ONLY: stress_kin, enkin, ennl_new
      USE fft_interfaces,     ONLY: fwfft, invfft
      USE fft_base,           ONLY: dffts, dfftp, dfft3d, dtgs
      USE cp_interfaces,      ONLY: checkrho, calrhovan 
      USE cp_main_variables,  ONLY: iprint_stdout, descla
      USE wannier_base,       ONLY: iwf
      USE exx_module,         ONLY: rhopr 
      USE input_parameters,   ONLY: tcpbo ! BS
      !
      IMPLICIT NONE
      !
      INTEGER  nbspx_bgrp, nbsp_bgrp, ispin_bgrp(:)
      REAL(DP) f_bgrp(:) 
      INTEGER  nfi
      REAL(DP) bec_bgrp(:,:)
      REAL(DP) dbec(:,:,:,:)
      REAL(DP) rhovan(:, :, : )
      REAL(DP) rhor(:,:)
      REAL(DP) drhor(:,:,:,:)
      REAL(DP) rhos(:,:)
      REAL(DP) enl, ekin
      REAL(DP) denl(3,3), dekin(6)
      COMPLEX(DP) eigrb( :, : )
      COMPLEX(DP) rhog( :, : )
      COMPLEX(DP) drhog( :, :, :, : )
      COMPLEX(DP) c_bgrp( :, : )
      INTEGER irb( :, : )
      LOGICAL, OPTIONAL, INTENT(IN) :: tstress
      INTEGER, OPTIONAL, INTENT(IN) :: ndwwf

      ! local variables

      INTEGER  :: iss, isup, isdw, iss1, iss2, ios, i, ir, ig, k
      REAL(DP) :: rsumr(2), rsumg(2), sa1, sa2, detmp(6), mtmp(3,3)
      REAL(DP) :: rnegsum, rmin, rmax, rsum
      COMPLEX(DP) :: ci,fp,fm
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: psi, psis, drhovan
#endif
#endif
      COMPLEX(DP), ALLOCATABLE :: psi(:), psis(:)
      REAL(DP), ALLOCATABLE :: drhovan(:,:,:,:,:)

      LOGICAL, SAVE :: first = .TRUE.
      LOGICAL :: ttstress
      !
      CALL start_clock( 'rhoofr_generalized' )

      ttstress = tpre
      IF( PRESENT( tstress ) ) ttstress = tstress

      ci = ( 0.0d0, 1.0d0 )

      rhor = 0.d0
      rhos = 0.d0
      rhog = (0.d0, 0.d0)
      !
      !  calculation of kinetic energy ekin
      !
      ekin = enkin( c_bgrp, f_bgrp, nbsp_bgrp )
      !
      IF( nbgrp > 1 ) &
         CALL mp_sum( ekin, inter_bgrp_comm )
      !
      IF( ttstress ) THEN
         !
         ! ... compute kinetic energy contribution
         !
         CALL stress_kin( dekin, c_bgrp, f_bgrp )
         !
         IF( nbgrp > 1 ) &
            CALL mp_sum( dekin, inter_bgrp_comm )
         !
      END IF

      IF( PRESENT( ndwwf ) ) THEN
         !
         !     called from WF, compute only of rhovan
         !
         CALL calrhovan( rhovan, bec_bgrp, iwf )
         !
      ELSE
         !
         !     calculation of non-local energy
         !
         CALL ennl_new( enl, rhovan, bec_bgrp, nbspx_bgrp, nbsp_bgrp, ispin_bgrp, f_bgrp )
         !
         IF( nbgrp > 1 ) THEN
            CALL mp_sum( enl, inter_bgrp_comm )
            CALL mp_sum( rhovan, inter_bgrp_comm )
         END IF
         !
      END IF
      !
      !IF( ttstress ) THEN
         !
      !   ALLOCATE( drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 ) )
         !
      !   CALL dennl( bec_bgrp, dbec, drhovan, denl, descla ) 
      !   !
      !   IF( nbgrp > 1 ) THEN
      !      CALL mp_sum( denl, inter_bgrp_comm )
      !      CALL mp_sum( drhovan, inter_bgrp_comm )
      !   END IF
         !
      !END IF
      !    
      !    warning! trhor and thdyn are not compatible yet!   
      !
      COMPUTE_CHARGE: IF( trhor .AND. ( .NOT. thdyn ) ) THEN
         !
         !     ==================================================================
         !     non self-consistent charge: charge density is read from unit 47
         !     ==================================================================
         !
         ! Lingzhu Kong
         !
         IF( first ) THEN
            CALL read_rho( nspin, rhor )
            rhopr = rhor
            first = .FALSE.
         ELSE
            rhor = rhopr
         END IF

         ALLOCATE( psi( dfftp%nnr ) )

         IF(nspin.EQ.1)THEN
            iss=1
            DO ir=1,dfftp%nnr
               psi(ir)=CMPLX(rhor(ir,iss),0.d0,kind=DP)
            END DO
            CALL fwfft('Dense', psi, dfftp )
            DO ig=1,ngm
               rhog(ig,iss)=psi(nl(ig))
            END DO
         ELSE
            isup=1
            isdw=2
            DO ir=1,dfftp%nnr
               psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw),kind=DP)
            END DO
            CALL fwfft('Dense', psi, dfftp )
            DO ig=1,ngm
               fp=psi(nl(ig))+psi(nlm(ig))
               fm=psi(nl(ig))-psi(nlm(ig))
               rhog(ig,isup)=0.5d0*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
               rhog(ig,isdw)=0.5d0*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
            END DO
         ENDIF

         DEALLOCATE( psi )

      ELSE
         !
         !     ==================================================================
         !     self-consistent charge
         !     ==================================================================
         !
         !     important: if n is odd then nx must be .ge.n+1 and c(*,n+1)=0.
         ! 

         IF ( MOD( nbsp_bgrp, 2 ) /= 0 ) THEN
            !
            IF( SIZE( c_bgrp, 2 ) < nbsp_bgrp + 1 ) &
               CALL errore( ' rhoofr_generalized ', ' c second dimension too small ', SIZE( c_bgrp, 2 ) )
            !
            c_bgrp( :, nbsp_bgrp + 1 ) = ( 0.d0, 0.d0 )
            !
         ENDIF
         !
         IF( PRESENT( ndwwf ) ) THEN
            !
            ! Wannier function, charge density from state iwf
            !
            ALLOCATE( psis( dffts%nnr ) ) 
            !
            i = iwf
            !
            psis = 0.D0
            DO ig=1,ngw
               psis(nlsm(ig))=CONJG(c_bgrp(ig,i))
               psis(nls(ig))=c_bgrp(ig,i)
            END DO
            !
            CALL invfft('Wave',psis, dffts )
            !
            iss1=1
            sa1=f_bgrp(i)/omega
            DO ir=1,dffts%nnr
               rhos(ir,iss1)=rhos(ir,iss1) + sa1*( DBLE(psis(ir)))**2
            END DO
            !
            DEALLOCATE( psis )
            !
         ELSE 
            !
            CALL loop_over_states()
            !
         END IF
         !
         !     smooth charge in g-space is put into rhog(ig)
         !
         ALLOCATE( psis( dffts%nnr ) ) 
         !
         IF(nspin.EQ.1)THEN
            iss=1
!$omp parallel do
            DO ir=1,dffts%nnr
               psis(ir)=CMPLX(rhos(ir,iss),0.d0,kind=DP)
            END DO
!$omp end parallel do
            CALL fwfft('Smooth', psis, dffts )
!$omp parallel do
            DO ig=1,ngms
               rhog(ig,iss)=psis(nls(ig))
            END DO
!$omp end parallel do
         ELSE
            isup=1
            isdw=2
!$omp parallel do
             DO ir=1,dffts%nnr
               psis(ir)=CMPLX(rhos(ir,isup),rhos(ir,isdw),kind=DP)
            END DO
!$omp end parallel do
            CALL fwfft('Smooth',psis, dffts )
!$omp parallel do private(fp,fm)
            DO ig=1,ngms
               fp= psis(nls(ig)) + psis(nlsm(ig))
               fm= psis(nls(ig)) - psis(nlsm(ig))
               rhog(ig,isup)=0.5d0*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
               rhog(ig,isdw)=0.5d0*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
            END DO
!$omp end parallel do
         ENDIF
         !
         ALLOCATE( psi( dfftp%nnr ) )
         !
         IF( nspin .EQ. 1 ) THEN
            ! 
            !     case nspin=1
            ! 
            iss=1
            psi (:) = (0.d0, 0.d0)
!$omp parallel do
            DO ig=1,ngms
               psi(nlm(ig))=CONJG(rhog(ig,iss))
               psi(nl (ig))=      rhog(ig,iss)
            END DO
!$omp end parallel do
            CALL invfft('Dense',psi, dfftp )
!$omp parallel do
            DO ir=1,dfftp%nnr
               rhor(ir,iss)=DBLE(psi(ir))
            END DO
!$omp end parallel do
            !
         ELSE 
            !
            !     case nspin=2
            !
            isup=1
            isdw=2
            psi (:) = (0.d0, 0.d0)
!$omp parallel do
            DO ig=1,ngms
               psi(nlm(ig))=CONJG(rhog(ig,isup))+ci*CONJG(rhog(ig,isdw))
               psi(nl(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            END DO
!$omp end parallel do
            CALL invfft('Dense',psi, dfftp )
!$omp parallel do
            DO ir=1,dfftp%nnr
               rhor(ir,isup)= DBLE(psi(ir))
               rhor(ir,isdw)=AIMAG(psi(ir))
            END DO
!$omp end parallel do
         ENDIF
         !
         IF ( dft_is_meta() ) CALL kedtauofr_meta( c_bgrp, psi, SIZE( psi ), psis, SIZE( psis ) ) ! METAGGA
         !
         DEALLOCATE( psi ) 
         DEALLOCATE( psis ) 
         !
         !     add vanderbilt contribution to the charge density
         !     drhov called before rhov because input rho must be the smooth part
         !
         IF ( ttstress ) THEN
            CALL drhov( irb, eigrb, rhovan, drhovan, rhog, rhor, drhog, drhor )
            DEALLOCATE( drhovan )
         END IF
         !
         CALL rhov( irb, eigrb, rhovan, rhog, rhor )

      ENDIF COMPUTE_CHARGE
!
      IF( PRESENT( ndwwf ) ) THEN
         !
         CALL old_write_rho( ndwwf, nspin, rhor )
         !
      END IF
!
!     here to check the integral of the charge density
!
! BS: I have turned off computing and printing integrated electronic density at
! every iprint_stdout steps during CP-BO calculations ... 
!     IF( ( iverbosity > 1 ) .OR. ( nfi == 0 ) .OR. &
!         ( MOD(nfi, iprint_stdout) == 0 ) .AND. ( .NOT. tcg ) ) THEN
      IF( ( iverbosity > 1 ) .OR. ( nfi == 0 ) .OR. &
          ( MOD(nfi, iprint_stdout) == 0 ) .AND. ( .NOT. tcg ) .AND. (.NOT.tcpbo )) THEN

         IF( iverbosity > 1 ) THEN
            CALL checkrho( dfftp%nnr, nspin, rhor, rmin, rmax, rsum, rnegsum )
            rnegsum = rnegsum * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
            rsum    = rsum    * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
            WRITE( stdout,'(a,4(1x,f12.6))')                                     &
     &     ' rhoofr_generalized: rmin rmax rnegsum rsum  ',rmin,rmax,rnegsum,rsum
         END IF

         CALL sum_charge( rsumg, rsumr )

         IF ( nspin == 1 ) THEN
           WRITE( stdout, 10) rsumg(1), rsumr(1)
         ELSE
           WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
         ENDIF

      ENDIF

10    FORMAT( /, 3X, 'from rhoofr_generalized: total integrated electronic density', &
            & /, 3X, 'in g-space = ', f13.6, 3x, 'in r-space =', f13.6 )
20    FORMAT( /, 3X, 'from rhoofr_generalized: total integrated electronic density', &
            & /, 3X, 'spin up', &
            & /, 3X, 'in g-space = ', f13.6, 3x, 'in r-space =', f13.6 , &
            & /, 3X, 'spin down', &
            & /, 3X, 'in g-space = ', f13.6, 3x, 'in r-space =', f13.6 )
!
      CALL stop_clock( 'rhoofr_generalized' )

!
      RETURN


   CONTAINS   
      !
      !
      SUBROUTINE sum_charge( rsumg, rsumr )
         !
         REAL(DP), INTENT(OUT) :: rsumg( : )
         REAL(DP), INTENT(OUT) :: rsumr( : )
         INTEGER :: iss
         !
         DO iss=1,nspin
            rsumg(iss)=omega*DBLE(rhog(1,iss))
            rsumr(iss)=SUM(rhor(:,iss),1)*omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         END DO

         IF (gstart.NE.2) THEN
            ! in the parallel case, only one processor has G=0 !
            DO iss=1,nspin
               rsumg(iss)=0.0d0
            END DO
         END IF

         CALL mp_sum( rsumg( 1:nspin ), intra_bgrp_comm )
         CALL mp_sum( rsumr( 1:nspin ), intra_bgrp_comm )

         RETURN
      END SUBROUTINE

      !

      SUBROUTINE loop_over_states
         !
         USE parallel_include
         USE fft_parallel,           ONLY: pack_group_sticks, fw_tg_cft3_z, fw_tg_cft3_scatter, fw_tg_cft3_xy
         USE fft_scalar, ONLY: cfft3ds
         USE scatter_mod, ONLY: maps_sticks_to_3d
         !
         !        MAIN LOOP OVER THE EIGENSTATES
         !           - This loop is also parallelized within the task-groups framework
         !           - Each group works on a number of eigenstates in parallel
         !
         IMPLICIT NONE
         !
         INTEGER :: from, i, eig_index, eig_offset, ii
         !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: tmp_rhos, aux
#endif
#endif
         REAL(DP), ALLOCATABLE :: tmp_rhos(:,:)
         COMPLEX(DP), ALLOCATABLE :: aux(:)

         ALLOCATE( psis( dtgs%tg_nnr * dtgs%nogrp ) ) 
         ALLOCATE( aux( dtgs%tg_nnr * dtgs%nogrp ) ) 
         !
         ALLOCATE( tmp_rhos ( dffts%nr1x * dffts%nr2x * dtgs%tg_npp( me_bgrp + 1 ), nspin ) )
         !
         tmp_rhos = 0_DP


         do i = 1, nbsp_bgrp, 2*dtgs%nogrp
            !
            !  Initialize wave-functions in Fourier space (to be FFTed)
            !  The size of psis is nnr: which is equal to the total number
            !  of local fourier coefficients.
            !

#if defined(__MPI)

            aux = (0.d0, 0.d0)
            !
            !  Loop for all local g-vectors (ngw)
            !  ci_bgrp: stores the Fourier expansion coefficients
            !     the i-th column of c_bgrp corresponds to the i-th state (in
            !     this band group)
            !  nlsm and nls matrices: hold conversion indices form 3D to
            !     1-D vectors. Columns along the z-direction are stored contigiously
            !
            !  The outer loop goes through i : i + 2*NOGRP to cover
            !  2*NOGRP eigenstates at each iteration
            !
            eig_offset = 0

!!$omp  parallel
!!$omp  single

            do eig_index = 1, 2*dtgs%nogrp, 2   
               !
!!$omp task default(none) &
!!$omp          firstprivate( i, eig_offset, nbsp_bgrp, ngw, eig_index  ) &
!!$omp          shared(  aux, c_bgrp, dffts )
               !
               !  here we pack 2*nogrp electronic states in the psis array
               !  note that if nogrp == nproc_bgrp each proc perform a full 3D
               !  fft and the scatter phase is local (without communication)
               !
               IF ( ( i + eig_index - 1 ) <= nbsp_bgrp ) THEN
                  !
                  !  The  eig_index loop is executed only ONCE when NOGRP=1.
                  !
                  CALL c2psi( aux(eig_offset*dtgs%tg_nnr+1), dtgs%tg_nnr, &
                        c_bgrp( 1, i+eig_index-1 ), c_bgrp( 1, i+eig_index ), ngw, 2 )
                  !
               ENDIF
!!$omp end task
               !
               eig_offset = eig_offset + 1
               !
            end do

!!$omp  end single
!!$omp  end parallel
            !
            !  2*NOGRP are trasformed at the same time
            !  psis: holds the fourier coefficients of the current proccesor
            !        for eigenstates i and i+2*NOGRP-1
            !
            !  now redistribute data
            !
            !
            IF( dtgs%nogrp == dtgs%nproc ) THEN
               CALL pack_group_sticks( aux, psis, dtgs )
               CALL maps_sticks_to_3d( dffts, dtgs, psis, SIZE(psis), aux, 2 )
               CALL cfft3ds( aux, dfft3d%nr1, dfft3d%nr2, dfft3d%nr3, &
                             dfft3d%nr1x,dfft3d%nr2x,dfft3d%nr3x, 1, 1, dfft3d%isind, dfft3d%iplw )
               psis = aux
            ELSE
               !
               CALL pack_group_sticks( aux, psis, dtgs )
               CALL fw_tg_cft3_z( psis, dffts, aux, dtgs )
               CALL fw_tg_cft3_scatter( psis, dffts, aux, dtgs )
               CALL fw_tg_cft3_xy( psis, dffts, dtgs )

            END IF
#else

            psis = (0.d0, 0.d0)

            CALL c2psi( psis, dffts%nnr, c_bgrp( 1, i ), c_bgrp( 1, i+1 ), ngw, 2 )

            CALL invfft('Wave', psis, dffts )

#endif
            !
            ! Now the first proc of the group holds the first two bands
            ! of the 2*nogrp bands that we are processing at the same time,
            ! the second proc. holds the third and fourth band
            ! and so on
            !
            ! Compute the proper factor for each band
            !
            DO ii = 1, dtgs%nogrp
               IF( dtgs%nolist( ii ) == me_bgrp ) EXIT
            END DO
            !
            ! Remember two bands are packed in a single array :
            ! proc 0 has bands ibnd   and ibnd+1
            ! proc 1 has bands ibnd+2 and ibnd+3
            ! ....
            !
            ii = 2 * ii - 1

            IF( ii + i - 1 < nbsp_bgrp ) THEN
               iss1=ispin_bgrp( ii + i - 1 )
               sa1 =f_bgrp( ii + i - 1 )/omega
               iss2=ispin_bgrp( ii + i )
               sa2 =f_bgrp( ii + i )/omega
            ELSE IF( ii + i - 1 == nbsp_bgrp ) THEN
               iss1=ispin_bgrp( ii + i - 1 )
               sa1 =f_bgrp( ii + i - 1 )/omega
               iss2=iss1
               sa2=0.0d0
            ELSE
               iss1=ispin_bgrp( nbsp_bgrp )
               sa1 = 0.0d0
               iss2=iss1
               sa2 =0.0d0
            END IF
            !
            !Compute local charge density
            !
            !This is the density within each orbital group...so it
            !coresponds to 1 eignestate for each group and there are
            !NOGRP such groups. Thus, during the loop across all
            !occupied eigenstates, the total charge density must me
            !accumulated across all different orbital groups.
            !

            !This loop goes through all components of charge density that is local
            !to each processor. In the original code this is nnr. In the task-groups
            !code this should be equal to the total number of planes
            !

            ir =  dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 ) 
            IF( ir > SIZE( psis ) ) &
               CALL errore( ' rhoofr_generalized ', ' psis size too small ', ir )

            do ir = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
               tmp_rhos(ir,iss1) = tmp_rhos(ir,iss1) + sa1*( real(psis(ir)))**2
               tmp_rhos(ir,iss2) = tmp_rhos(ir,iss2) + sa2*(aimag(psis(ir)))**2
            end do
            !
         END DO

         IF( nbgrp > 1 ) THEN
            CALL mp_sum( tmp_rhos, inter_bgrp_comm )
         END IF

         !ioff = 0
         !DO ip = 1, nproc_bgrp
         !   CALL MPI_REDUCE( rho(1+ioff*nr1*nr2,1), rhos(1,1), dffts%nnr, MPI_DOUBLE_PRECISION, MPI_SUM, ip-1, intra_bgrp_comm, ierr)
         !   ioff = ioff + dffts%npp( ip )
         !END DO
         IF ( dtgs%nogrp > 1 ) THEN
            CALL mp_sum( tmp_rhos, gid = dtgs%ogrp_comm )
         ENDIF
         !
         !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
         !
         !If the current processor is not the "first" processor in its
         !orbital group then does a local copy (reshuffling) of its data
         !
         from = 1
         DO ii = 1, dtgs%nogrp
            IF ( dtgs%nolist( ii ) == me_bgrp ) EXIT !Exit the loop
            from = from +  dffts%nr1x*dffts%nr2x*dffts%npp( dtgs%nolist( ii ) + 1 )! From where to copy initially
         ENDDO
         !
         DO ir = 1, nspin
            CALL dcopy( dffts%nr1x*dffts%nr2x*dffts%npp(me_bgrp+1), tmp_rhos(from,ir), 1, rhos(1,ir), 1)
         ENDDO

         DEALLOCATE( tmp_rhos )
         DEALLOCATE( aux ) 
         DEALLOCATE( psis ) 

         RETURN
      END SUBROUTINE loop_over_states
      ! 
END SUBROUTINE rhoofr_cp_generalized_x 

SUBROUTINE ennl_new_x( ennl_val, rhovan, bec_bgrp, nbspx_bgrp, nbsp_bgrp, ispin_bgrp, f_bgrp)
      !
      ! calculation of nonlocal potential energy term and array rhovan
      !
      use kinds,          only : DP
      use uspp_param,     only : nh, ish
      use uspp,           only : dvan
      use electrons_base, only : nspin 
      use ions_base,      only : nsp, na
      !
      implicit none
      !
      ! input
      !
      integer,  intent(in)  :: nbspx_bgrp, nbsp_bgrp, ispin_bgrp(:)
      real(DP), intent(in)  :: f_bgrp(:)
      real(DP), intent(out) :: ennl_val
      real(DP), intent(out) :: rhovan( :, :, : )
      real(DP), intent(in)  :: bec_bgrp( :, : )
      !
      ! local
      !
      real(DP) :: sumt, sums(2), ennl_t
      integer  :: is, iv, jv, ijv, inl, jnl, isa, isat, ism, ia, iss, i
      !
      ennl_t = 0.d0
      !
      !  xlf does not like name of function used for OpenMP reduction
      !
!$omp parallel default(shared), &
!$omp private(is,iv,jv,ijv,isa,isat,ism,ia,inl,jnl,sums,i,iss,sumt),
!reduction(+:ennl_t)
      do is = 1, nsp
         do iv = 1, nh(is)
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               isa = 0
               do ism = 1, is - 1
                  isa = isa + na(ism)
               end do
!$omp do
               do ia = 1, na(is)
                  inl = ish(is)+(iv-1)*na(is)+ia
                  jnl = ish(is)+(jv-1)*na(is)+ia
                  isat = isa+ia
                  sums = 0.d0
                  do i = 1, nbsp_bgrp
                     iss = ispin_bgrp(i)
                     sums(iss) = sums(iss) + f_bgrp(i) * bec_bgrp(inl,i) * bec_bgrp(jnl,i)
                  end do
                  sumt = 0.d0
                  do iss = 1, nspin
                     rhovan( ijv, isat, iss ) = sums( iss )
                     sumt = sumt + sums( iss )
                  end do
                  if( iv .ne. jv ) sumt = 2.d0 * sumt
                  ennl_t = ennl_t + sumt * dvan( jv, iv, is)
               end do
!$omp end do
            end do
         end do
      end do
!$omp end parallel
      !
      ennl_val = ennl_t
      !
      return
      !
endsubroutine ennl_new_x
