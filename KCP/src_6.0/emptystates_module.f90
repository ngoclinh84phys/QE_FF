
! THIS FILE CONTENTS THE SUBROUTINE USED IN EMPTY STATES CALCULATIONS 

SUBROUTINE gram_empty ( tortho, eigr, betae, bec_emp, bec_occ, nkbx, c_emp, c_occ, ngwx, n_emp, n_occ )
   !     
   !   gram-schmidt orthogonalization of the empty states ( c_emp ) 
   !   c_emp are orthogonalized among themself
   !   and to the occupied states c_occ
   !
   use io_global,      only: stdout
   USE kinds,          ONLY : DP
   USE uspp,           ONLY : nkb, nkbus
   USE uspp_param,     ONLY : nvb
   USE gvecw,          ONLY : ngw
   USE mp,             ONLY : mp_sum
   USE mp_global,      ONLY : intra_bgrp_comm
   USE ions_base,      ONLY : nat
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: nkbx, ngwx, n_emp, n_occ
   COMPLEX(DP)   :: eigr(ngwx,nat)
   REAL(DP)      :: bec_emp( nkbx, n_emp )
   REAL(DP)      :: bec_occ( nkbx, n_occ )
   COMPLEX(DP)   :: betae( ngwx, nkb )
   COMPLEX(DP)   :: c_emp( ngwx, n_emp )
   COMPLEX(DP)   :: c_occ( ngwx, n_occ )
   LOGICAL, INTENT(IN) :: tortho
   !
   REAL(DP) :: anorm, cscnorm
   REAL(DP), ALLOCATABLE :: csc_emp( : )
   REAL(DP), ALLOCATABLE :: csc_occ( : )
   INTEGER :: i, k, inl
   EXTERNAL cscnorm
   !
   ALLOCATE( csc_emp( n_emp ) )
   ALLOCATE( csc_occ( n_occ ) )
   !
   ! orthogonalize empty states to the occupied one and among each
   ! other
   !
   DO i = 1, n_emp
      !
      csc_emp = 0.0d0
      csc_occ = 0.0d0
      !
      ! compute scalar product csc_occ(k) = <c_emp(i)|c_occ(k)>
      !
      CALL smooth_csv( c_emp(1:,i), c_occ(1:,1:), ngwx, csc_occ, n_occ )
      !
      IF ( .NOT. tortho ) THEN
         !
         ! compute scalar product csc_emp(k) = <c_emp(i)|c_emp(k)>
         !
         CALL smooth_csv( c_emp(1:,i), c_emp(1:,1:), ngwx, csc_emp, i-1 )
         !
         CALL mp_sum( csc_emp, intra_bgrp_comm )
         !
      ENDIF
      !
      CALL mp_sum( csc_occ, intra_bgrp_comm )
      !
      IF ( nvb > 1 ) THEN
         !
         CALL grabec( bec_emp(1:,i), nkbx, betae, c_emp(1:,i), ngwx)
         !
         CALL mp_sum( bec_emp(1:nkbus,i), intra_bgrp_comm )
         !
         CALL bec_csv( bec_emp(1:,i), bec_occ, nkbx, csc_occ, n_occ )
         !
         IF ( .NOT. tortho ) THEN
            CALL bec_csv( bec_emp(1:,i), bec_emp, nkbx, csc_emp, i-1)
         ENDIF
         !
         DO k = 1, n_occ
            DO inl = 1, nkbx
               bec_emp( inl, i ) = bec_emp( inl, i ) - csc_occ(k)*bec_occ( inl, k )
            ENDDO
         ENDDO
         !
         IF ( .NOT. tortho ) THEN
            DO k = 1, i-1
               DO inl = 1, nkbx
                  bec_emp( inl, i ) = bec_emp( inl, i ) - csc_emp(k)*bec_emp( inl, k )
               ENDDO
            ENDDO
         ENDIF
         !
      ENDIF
      !
      ! calculate orthogonalized c_emp(i) : |c_emp(i)> = |c_emp(i)> -
      ! SUM_k    csv(k)|c_occ(k)>
      !                          c_emp(i) : |c_emp(i)> = |c_emp(i)> -
      !                          SUM_k<i  csv(k)|c_emp(k)>
      !
      DO k = 1, n_occ
         CALL DAXPY( 2*ngw, -csc_occ(k), c_occ(1,k), 1, c_emp(1,i), 1)
      ENDDO
      ! 
      IF ( .NOT. tortho ) THEN
         DO k = 1, i - 1
            CALL DAXPY( 2*ngw, -csc_emp(k), c_emp(1,k), 1, c_emp(1,i), 1 )
         END DO
      ENDIF
      !
      IF ( .NOT. tortho ) THEN
         !
         anorm = cscnorm( bec_emp, nkbx, c_emp, ngwx, i, n_emp )
         !
         CALL DSCAL( 2*ngw, 1.0d0/anorm, c_emp(1,i), 1 )
         !
         IF ( nvb > 1 ) THEN
            CALL DSCAL( nkbx, 1.0d0/anorm, bec_emp(1,i), 1 )
         ENDIF
      ENDIF
      !
   ENDDO
   !
   DEALLOCATE( csc_emp )
   DEALLOCATE( csc_occ )
   !
   RETURN
   !
   CONTAINS
     !   
     SUBROUTINE smooth_csv( c, v, ngwx, csv, n )
        ! 
        USE kinds,              ONLY: DP
        USE gvecw,              ONLY: ngw
        USE gvect,              ONLY: gstart
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: ngwx, n
        COMPLEX(DP)    :: c( ngwx )
        COMPLEX(DP)    :: v( ngwx, n )
        REAL(DP)    :: csv( n )
        INTEGER     :: k, ig
        REAL(DP), ALLOCATABLE :: temp(:)
        !  
        ! calculate csv(k)=<c|v(k)>
        !
        ALLOCATE( temp( ngw ) )
        !  
        DO k = 1, n
           !
           DO ig = 1, ngw
              !
              temp(ig) = CONJG(v(ig,k)) * c(ig)
              ! 
           ENDDO
           ! 
           csv(k) = 2.0d0 * SUM(temp)
           ! 
           IF (gstart == 2) csv(k) = csv(k) - temp(1)
           !
        ENDDO
        !
        DEALLOCATE( temp )
        !  
        RETURN
        !
     ENDSUBROUTINE smooth_csv
     !
     SUBROUTINE grabec( becc, nkbx, betae, c, ngwx )
        !
        !  on output: bec(i) is recalculated
        !
        USE kinds,          ONLY: DP
        USE uspp,           ONLY: nkb, nkbus
        USE gvecw,          ONLY: ngw
        USE gvect,          ONLY: gstart
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nkbx, ngwx
        COMPLEX(DP) :: betae( ngwx, nkb )
        REAL(DP)    :: becc( nkbx )
        COMPLEX(DP) :: c( 1, ngwx ) 
        INTEGER     :: ig, inl
        REAL(DP), ALLOCATABLE :: temp(:)
        !
        ALLOCATE( temp( ngw ) )
        !
        !     calculate becc=<c|beta>
        !
        DO inl=1, nkbus
           DO ig=1,ngw
              temp(ig)=DBLE(CONJG(c(1,ig))* betae(ig,inl))
           ENDDO
           ! 
           becc(inl)=2.d0*SUM(temp)
           !
           IF (gstart == 2) becc(inl)= becc(inl)-temp(1)
           !
        ENDDO
        ! 
        DEALLOCATE( temp )
        ! 
        RETURN
        !  
     ENDSUBROUTINE grabec
     !
     SUBROUTINE bec_csv( becc, becv, nkbx, csv, n )
        !     
        ! requires in input the updated becc and becv(k)
        ! on output: csv is updated
        !
        USE kinds,          ONLY: DP
        USE ions_base,      ONLY: na
        USE uspp,           ONLY: qq
        USE uspp_param,     ONLY: nh, nvb, ish
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nkbx, n
        REAL(DP)    :: becc( nkbx )
        REAL(DP)    :: becv( nkbx, n )
        REAL(DP)    :: csv( n )
        INTEGER     :: k, is, iv, jv, ia, inl, jnl
        REAL(DP)    :: rsum
        ! 
        ! calculate csv(k) = csv(k) + <c| SUM_nm |beta(n)><beta(m)|v(k)>,  k<i
        !
        DO k=1,n
           rsum=0.d0
           DO is=1,nvb
              DO iv=1,nh(is)
                 DO jv=1,nh(is)
                    IF (ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                       DO ia=1,na(is)
                          inl=ish(is)+(iv-1)*na(is)+ia
                          jnl=ish(is)+(jv-1)*na(is)+ia
                          rsum = rsum + qq(iv,jv,is)*becc(inl)*becv(jnl,k)
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
           !
           csv(k)=csv(k)+rsum
           !
        ENDDO
        ! 
        RETURN
        ! 
     ENDSUBROUTINE bec_csv
     !
ENDSUBROUTINE gram_empty
!
SUBROUTINE gram_standard_empty( betae, bec, nkbx, cp, ngwx, n )
     ! 
     ! Gram-schmidt orthogonalization of the set of wavefunctions cp
     !
     USE uspp,           ONLY : nkb, nhsavb=> nkbus
     USE gvecw,          ONLY : ngw
     USE kinds,          ONLY : DP
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nkbx, ngwx, n
     REAL(DP)      :: bec( nkbx, n )
     COMPLEX(DP)   :: cp( ngwx, n ), betae( ngwx, nkb )
     ! 
     REAL(DP)      :: anorm, cscnorm
     COMPLEX(DP), ALLOCATABLE :: csc( : )
     INTEGER :: i,k
     EXTERNAL cscnorm
     !
     CALL start_clock( 'gram_empty' )
     !  
     ALLOCATE( csc( n ) )
     !
     DO i = 1, n
        !
        CALL gracsc( bec, nkbx, betae, cp, ngwx, i, csc, n )
        !
        ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i
        ! csc(k)|cp(k)>
        !
        DO k = 1, i - 1
           CALL DAXPY( 2*ngw, -csc(k), cp(1,k), 1, cp(1,i), 1 )
        ENDDO
        !
        anorm = cscnorm( bec, nkbx, cp, ngwx, i, n )
        !
        CALL DSCAL( 2*ngw, 1.0d0/anorm, cp(1,i), 1 )
        !
        ! these are the final bec's
        !
        CALL DSCAL( nkbx, 1.0d0/anorm, bec(1,i), 1 )
        ! 
     ENDDO
     !
     DEALLOCATE( csc )
     ! 
     CALL stop_clock( 'gram_empty' )
     !
     RETURN
     !
     CONTAINS
        !
        SUBROUTINE gracsc( bec, nkbx, betae, cp, ngwx, i, csc, n )
           !
           ! requires in input the updated bec(k) for k<i
           ! on output: bec(i) is recalculated
           !
           USE ions_base,      ONLY: na
           USE uspp,           ONLY: nkb, nhsavb=>nkbus, qq
           USE uspp_param,     ONLY: nh, nvb, ish
           USE electrons_base, ONLY: ispin
           USE gvecw,          ONLY: ngw 
           USE mp,             ONLY: mp_sum
           USE mp_global,      ONLY: intra_bgrp_comm 
           USE kinds,          ONLY: DP
           USE gvect,       ONLY : gstart
           !
           IMPLICIT NONE
           !
           INTEGER, INTENT(IN) :: i, nkbx, ngwx, n
           COMPLEX(DP):: betae( ngwx, nkb )
           REAL(DP):: bec( nkbx, n )
           COMPLEX(DP) :: cp(ngwx, n)
           COMPLEX(DP) :: csc( n ) 
           INTEGER :: k, kmax,ig, is, iv, jv, ia, inl, jnl
           REAL(DP)    :: rsum
           REAL(DP), ALLOCATABLE :: temp(:)
           !
           ! calculate csc(k)=<cp(i)|cp(k)>,  k<i
           ! 
           kmax = i - 1
           !
           ALLOCATE( temp( ngw ) )
           !
           temp = 0.d0
           DO k = 1, kmax
              csc(k) = CMPLX(0.0d0, 0.d0)
              IF ( ispin(i) .EQ. ispin(k) ) THEN
                 !
                 DO ig = 1, ngw
                    temp(ig) = DBLE(CONJG(cp(ig,i))*cp(ig,k)) 
                 ENDDO
                 csc(k) = CMPLX(2.d0* SUM(temp,ngw),0.d0)
                 IF (gstart == 2) csc(k) = csc(k) - CMPLX(temp(1), 0.d0)
              ENDIF
           ENDDO
           !
           DEALLOCATE( temp ) 
           !
           CALL mp_sum( csc( 1:kmax ), intra_bgrp_comm )
           !  
           ALLOCATE( temp( ngw ) )
           !
           ! calculate bec(i)=  <beta|cp(i)> 
           !
           DO inl=1,nhsavb
              DO ig=1,ngw
                 temp(ig)=DBLE(cp(ig,i)* CONJG(betae(ig,inl))) 
              ENDDO
              bec(inl,i)=2.d0*SUM(temp)
              IF (gstart == 2) bec(inl,i)= bec(inl,i)-temp(1)
           ENDDO
           ! 
           CALL mp_sum( bec( 1:nhsavb, i ), intra_bgrp_comm )
           !
           ! calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
           !
           DO k=1,kmax
              IF (ispin(i).EQ.ispin(k)) THEN
                 rsum=0.d0
                 DO is=1,nvb
                    DO iv=1,nh(is)
                       DO jv=1,nh(is)
                          IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                             DO ia=1,na(is)
                                inl=ish(is)+(iv-1)*na(is)+ia
                                jnl=ish(is)+(jv-1)*na(is)+ia
                                rsum = rsum + qq(iv,jv,is)*bec(inl,i)*bec(jnl,k)
                             ENDDO
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
                 csc(k)=csc(k)+ CMPLX(rsum,0.d0)
              ENDIF
           ENDDO
           !
           ! orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
           ! corresponing bec:  bec(i)=<cp(i)|beta>-csc(k)<cp(k)|beta>
           !
           DO k=1,kmax
              DO inl=1,nkbx
                 bec(inl,i)=bec(inl,i)-DBLE(csc(k))*bec(inl,k)
              ENDDO
           ENDDO
           !
           DEALLOCATE( temp )
           !
           RETURN
           !
        ENDSUBROUTINE gracsc
        !
ENDSUBROUTINE gram_standard_empty
!
FUNCTION cscnorm( bec, nkbx, cp, ngwx, i, n)
        !
        ! Compute the norm of the i-th electronic state = (<c|S|c>)^(1/2) 
        ! requires in input the updated bec(i)
        !
        USE kinds,              ONLY: DP
        USE ions_base,          ONLY: na
        USE gvecw,              ONLY: ngw
        USE gvect,              ONLY: gstart
        USE uspp_param,         ONLY: nh, ish, nvb
        USE uspp,               ONLY: qq, nkb
        USE mp,                 ONLY: mp_sum
        USE mp_global,          ONLY: intra_bgrp_comm
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: i, n
        INTEGER, INTENT(IN) :: ngwx, nkbx
        REAL(DP),    INTENT(IN) :: bec( nkb, n )
        COMPLEX(DP), INTENT(IN) :: cp( ngw, n )
        !
        REAL(DP) :: cscnorm
        !
        INTEGER ig, is, iv, jv, ia, inl, jnl
        REAL(DP) rsum
        REAL(DP), ALLOCATABLE:: temp(:)
        !
        ALLOCATE(temp(ngw))
        !
        DO ig=1,ngw
           temp(ig)=DBLE(CONJG(cp(ig,i))*cp(ig,i))
        ENDDO
        ! 
        rsum=2.d0*SUM(temp)
        !
        IF (gstart == 2) rsum=rsum-temp(1)
        !
        CALL mp_sum( rsum, intra_bgrp_comm )
        ! 
        DEALLOCATE(temp)
        !
        DO is=1,nvb
           DO iv=1,nh(is)
              DO jv=1,nh(is)
                 IF (ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                    DO ia=1,na(is)
                       inl=ish(is)+(iv-1)*na(is)+ia
                       jnl=ish(is)+(jv-1)*na(is)+ia
                       rsum=rsum + qq(iv,jv,is)*bec(inl,i)*bec(jnl,i)
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !
        cscnorm=SQRT(rsum)
        !
        RETURN
        ! 
ENDFUNCTION cscnorm

subroutine pc2_generalized(a,beca,b,becb, n, ispin, nupdwn, iupdwn)      
               
! this function applies the operator Pc
            
!    this subroutine applies the Pc operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i-a_j><a_j|S|b_i>
    
      use kinds, only: dp 
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum
      use electrons_base, only: nspin
      use uspp_param, only: nh, nvb, ish
      use uspp, only :nhsa=>nkb
      use uspp, only :qq
      use parallel_toolkit, only : rep_matmul_drv
      
                           
      implicit none        
                          
      integer :: n, ispin(n), nupdwn(nspin), iupdwn(nspin) 
      complex(kind=DP) a(ngw,n), b(ngw,n)
                     
      real(kind=DP)    beca(nhsa,n),becb(nhsa,n)
! local variables
      integer is, iv, jv, ia, inl, jnl, i, j,ig
      real(kind=DP) sca
      real(DP), allocatable :: bectmp(:,:)
      real(DP), allocatable :: qq_tmp(:,:), qqb_tmp(:,:)
      complex(DP), allocatable :: zbectmp(:,:)
      integer :: nl_max
      integer :: nss,iss, istart

      logical :: mat_par=.true.!if true uses parallel routines

      CALL start_clock( 'pc2_generalized' )

      do iss= 1, nspin
         nss= nupdwn( iss )
         istart= iupdwn( iss )

         allocate(bectmp(nss,nss))
         allocate(zbectmp(nss,nss))
         bectmp(:,:)=0.d0

         call zgemm('C','N',nss,nss,ngw,(1.d0,0.d0),a(:,istart),ngw,b(:,istart),ngw,(0.d0,0.d0),zbectmp,nss)

         do j=1,nss
            do i=1,nss
               bectmp(i,j)=2.d0*dble(zbectmp(i,j))
               if(gstart==2) bectmp(i,j)=bectmp(i,j)-DBLE(a(1,j))*DBLE(b(1,i))
               
            enddo
         enddo
         deallocate(zbectmp)
         call mp_sum( bectmp(:,:), intra_bgrp_comm)
         if(nvb >= 0) then

            nl_max=0
            do is=1,nvb
               nl_max=nl_max+nh(is)*na(is)
            enddo
            allocate (qq_tmp(nl_max,nl_max))
            allocate (qqb_tmp(nl_max,nss))
            qq_tmp(:,:)=0.d0
            do is=1,nvb
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     do ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        qq_tmp(inl,jnl)=qq(iv,jv,is)
                     enddo
                  enddo
               enddo
            enddo
            if(.not. mat_par)  then
               call dgemm('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb(:,istart),nhsa,0.d0,qqb_tmp,nl_max)
               call dgemm('T','N',nss,nss,nl_max,1.d0,beca(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss)
            else
               call para_dgemm & 
('N','N',nl_max,nss,nl_max,1.d0,qq_tmp,nl_max,becb(:,istart),nhsa,0.d0,qqb_tmp,nl_max, intra_bgrp_comm)
               call para_dgemm &
('T','N',nss,nss,nl_max,1.d0,beca(:,istart),nhsa,qqb_tmp,nl_max,1.d0,bectmp,nss, intra_bgrp_comm)
            endif
            deallocate(qq_tmp,qqb_tmp)
         endif
         allocate(zbectmp(nss,nss))
         do i=1,nss
            do j=1,nss
               zbectmp(i,j)=CMPLX(bectmp(i,j),0.d0,kind=dp)
            enddo
         enddo
         call zgemm('N','N',ngw,nss,nss,(-1.d0,0.d0),a(:,istart),ngw,zbectmp,nss,(1.d0,0.d0),b(:,istart),ngw)
         deallocate(zbectmp)
         call dgemm('N','N',nhsa,nss,nss,1.0d0,beca(:,istart),nhsa,bectmp,nss,1.0d0,becb(:,istart),nhsa)
         deallocate(bectmp)
      enddo!on spin
      CALL stop_clock( 'pc2_generalized' )
      return
    end subroutine pc2_generalized


    subroutine pcdaga2_generalized (a,as ,b, n, ispin )

! this function applies the operator Pc

!    this subroutine applies the Pc^dagerr operator
!    a input :unperturbed wavefunctions
!    b input :first order wavefunctions
!    b output:b_i =b_i - S|a_j><a_j|b_i>

      use kinds
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use gvect, only: gstart
      use mp, only: mp_sum
      use uspp_param, only: nh, ish, nvb
      use uspp, only :nhsa=>nkb


      implicit none

      integer, intent(in) :: n, ispin(n)
      complex(dp) a(ngw,n), b(ngw,n), as(ngw,n)
      ! local variables
      integer is, iv, jv, ia, inl, jnl, i, j,ig
      real(dp) sca
      real(DP), allocatable:: scar(:)
      !
      call start_clock('pcdaga2_generalized')
      allocate(scar(n))
      do j=1,n
         do i=1,n
            sca=0.0d0
            if(ispin(i) == ispin(j)) then
               if (gstart==2) b(1,i) = CMPLX(dble(b(1,i)),0.0d0,kind=dp)
               do  ig=1,ngw           !loop on g vectors
                  sca=sca+DBLE(CONJG(a(ig,j))*b(ig,i))
               enddo
               sca = sca*2.0d0  !2. for real weavefunctions
               if (gstart==2) sca = sca - dble(a(1,j))*dble(b(1,i))
            endif
            scar(i) = sca
         enddo
         
          
         call mp_sum( scar, intra_bgrp_comm )


         do i=1,n
            if(ispin(i) == ispin(j)) then
               sca = scar(i)
               do ig=1,ngw
                  b(ig,i)=b(ig,i)-sca*as(ig,j)
               enddo
               ! this to prevent numerical errors
               if (gstart==2) b(1,i) = CMPLX(dble(b(1,i)),0.0d0,kind=dp)
            endif
         enddo
      enddo
      deallocate(scar)
      call stop_clock('pcdaga2_generalized')
      return
      end subroutine pcdaga2_generalized

subroutine xminus1_generalized(n, c0,betae,ema0bg,beck,m_minus1,do_k)
! if (do_k) then
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = K^{-1}|c0>
! else
!-----------------------------------------------------------------------
!     input: c0 , bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s^{-1}|c0> 
! endif
      use kinds, only: dp
      use ions_base, only: na, nsp
      use io_global, only: stdout
      use mp_global, only: intra_bgrp_comm
      use uspp_param, only: nh, nvb, ish
      use uspp, only :nhsa=>nkb, nhsavb=>nkbus, qq
      use gvecw, only: ngw
      use constants, only: pi, fpi
      use mp, only: mp_sum
      use gvect, only: gstart
!
      implicit none
      integer :: n
      complex(dp) c0(ngw,n), betae(ngw,nhsa)
      real(dp)    beck(nhsa,n), ema0bg(ngw)
      real(DP)    :: m_minus1(nhsavb,nhsavb)
      logical :: do_k
! local variables
      complex(dp), allocatable :: phi(:,:)
      real(dp) , allocatable   :: qtemp(:,:)
      integer is, iv, jv, ia, inl, jnl, i, j, js, ja,ig
      real(dp) becktmp

      
      logical :: mat_par=.true.!if true uses parallel routines      

      call start_clock('xminus1_generalized')
      if (nvb.gt.0) then
!calculates beck
         if (do_k) then
            beck(:,:) = 0.d0

            do is=1,nvb
               do iv=1,nh(is)
                  do ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        becktmp = 0.0d0
                        do ig=1,ngw
                           becktmp=becktmp+ema0bg(ig)*DBLE(CONJG(betae(ig,inl))*c0(ig,i))
                        enddo
                        becktmp = becktmp*2.0d0
                        if (gstart==2) becktmp = becktmp-ema0bg(1)*DBLE(CONJG(betae(1,inl))*c0(1,i)) 
                        beck(inl,i) = beck(inl,i) + becktmp
                     enddo
                  enddo
               enddo
            enddo
            call mp_sum( beck, intra_bgrp_comm )
         endif
!
!
      allocate(phi(ngw,n))
      allocate(qtemp(nhsavb,n))
      phi(1:ngw,1:n) = 0.0d0
      qtemp(:,:) = 0.0d0
      if(.not.mat_par) then
         call dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb ,    &
                    beck, nhsa, 0.0d0, qtemp,nhsavb )
      else
         call para_dgemm( 'N', 'N', nhsavb, n, nhsavb, 1.0d0, m_minus1,nhsavb , &
                    beck, nhsa, 0.0d0, qtemp,nhsavb,intra_bgrp_comm )
      endif

!NB  nhsavb is the total number of US projectors
!    it works because the first pseudos are the vanderbilt's ones

         CALL dgemm( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngw,    &
                    qtemp, nhsavb, 0.0d0, phi, 2*ngw )
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=(phi(ig,j)+c0(ig,j))*ema0bg(ig)
               end do
            end do
         else
            do j=1,n
               do i=1,ngw
                  c0(i,j)=(phi(i,j)+c0(i,j))
               end do
            end do
         endif
      deallocate(qtemp,phi)

      else
         if (do_k) then
            do j=1,n
               do ig=1,ngw
                  c0(ig,j)=c0(ig,j)*ema0bg(ig)
               end do
            end do
         endif
      endif
      call stop_clock('xminus1_generalized')
      return
     end subroutine xminus1_generalized
