subroutine wan2odd ( seedname, ispin, iknum, ikstart, extension)
  !
  use io_global,            only : stdout
  use kinds,                only : DP
  use io_files,             only : iunwfc, nwordwfc
  use gvect,                only : g, ngm
  use gvecw,                only : gcutw
  use wavefunctions_module, only : evc
  use wvfct,                only : nbnd, npw, g2kin
  use klist,                only : xk, igk_k, ngk
  use control_flags,        only : gamma_only
  USE mp,                   only : mp_bcast, mp_sum
  use mp_global,            only : intra_pool_comm
  use mp_images,            only : intra_image_comm
  use fft_base,             only : dfftp, dffts, dtgs
  use fft_interfaces,       only : fwfft, invfft
  use gvecs,                ONLY : nls, nlsm
  use klist,                ONLY : igk_k
  use mp,                   ONLY : mp_sum
  use read_wannier_files
  !
  IMPLICIT NONE 
  !
  character(len=50), intent (in) :: seedname
  character(len=*), intent (in)  :: extension
  integer, intent (in) :: iknum, ikstart, ispin
  !
  ! local variables
  ! 
  integer     :: ir, j, nn, ik, ibnd, iw, ikevc, counter, num_inc
  integer     :: loop_b, loop_w, jbnd
  real(DP)    :: summ1, summ2, summ3, summ4, max_value_ibnd
  real(DP)    :: summ_ib_r, summ_ib_i, summ_jb_r, summ_jb_i
  real (DP),    allocatable:: overlap_orbs(:,:)
  complex (DP), allocatable:: evc_disen(:,:)
  complex (DP), allocatable:: evc_tmp(:,:)
  complex (DP), allocatable:: evc_disen_tmp(:,:)
  complex (DP), allocatable:: c0_tmp(:,:)
  complex (DP), allocatable:: psic1(:), psic2(:), c0_in_r(:,:)
  logical,      allocatable:: inc_band(:),  exclude_band(:)
  ! 
  ! main program
  !
  call conv_read_chkpt(seedname)
  !
  write(stdout, *) " **** Transform KS orbitals to Wannier orbs **** " 
  !
  allocate (c0_tmp (npw, num_wann))
  c0_tmp(:,:) = (0.0, 0.0)
  !
  if (have_disentangled) then 
     allocate (inc_band ( maxval(ndimwin)) )
  else
     allocate (inc_band (num_bands))
  endif
  !
  do ik=1, iknum
     !
     inc_band(:) = .false.
     num_inc=num_wann
     if (have_disentangled) then
        !
        inc_band(:) = lwindow(:,ik)
        num_inc = ndimwin(ik)
        !
     endif
     !
     ikevc = ik + ikstart - 1
     call davcio (evc, nwordwfc, iunwfc, ikevc, -1)
     call gk_sort (xk(1,ik), ngm, g, gcutw, ngk(ik), igk_k(1,ik), g2kin)
     !
     allocate (exclude_band(nbnd))
     exclude_band(:) = .false.
     do ibnd = 1, num_exclude_bands
        exclude_band(exclude_bands(ibnd)) = .true.  
     enddo 
     !
     npw = ngk(ik) 
     !
     if (have_disentangled) allocate (evc_disen (npw, maxval(ndimwin)) )
     !
     allocate (evc_tmp (npw, num_bands))
     evc_tmp(:,:) = (0.0, 0.0)
     !
     counter = 0
     do ibnd = 1, nbnd
        !
        if (.not.exclude_band(ibnd)) then
           counter = counter + 1
           evc_tmp(:,counter) = evc(:, ibnd)
        endif
        ! 
     enddo 
     !
     if ( counter.ne.num_bands ) &
        call errore('wan2odd',' counter.ne.num_bands',1)
     !
     if (have_disentangled) then
        !
        counter=1
        ! 
        do loop_b=1, num_bands
           ! 
           if(counter > num_inc) exit
           !
           if (inc_band(loop_b)) then 
              !
              evc_disen (:, counter) = evc_tmp (:, loop_b)  
              !
              counter=counter + 1
              ! 
           endif 
           !
        enddo
        !
        if ( (counter-1).ne. num_inc ) &
           call errore('wan2odd',' counter.ne.num_inc',1)
        !
        allocate (evc_disen_tmp(npw, num_wann))
        evc_disen_tmp(:,:) = (0.0, 0.0)
        !  
        do loop_w=1, num_wann
           !    
           do loop_b=1, num_inc
              !
              evc_disen_tmp (:,loop_w) = evc_disen_tmp (:,loop_w) + &
                      u_matrix_opt(loop_b,loop_w, ik) * evc_disen (:,loop_b)   
              !
           enddo
           !   
        enddo
        ! 
     endif
     !
     ! Here is the core of this routine
     !   
     do ibnd = 1, num_wann
        !
        do jbnd = 1, num_wann
           !
           if (have_disentangled) then
              c0_tmp(:,jbnd) = c0_tmp(:,jbnd) + u_matrix(ibnd,jbnd, ik) * evc_disen_tmp(:,ibnd)
           else
              c0_tmp(:,jbnd) = c0_tmp(:,jbnd) + u_matrix(ibnd,jbnd, ik) * evc_tmp(:,ibnd)
           endif
           !
        enddo
        !
     enddo
     !
     if (have_disentangled) deallocate (evc_disen_tmp)
     if (have_disentangled) deallocate (evc_disen)
     deallocate (evc_tmp)
     deallocate (exclude_band)
     ! 
  enddo ! k-points
  !
  Write(stdout,*) "***** write wannier orbs and rotation matrix to file *****"
  Write(stdout,*) "*****"
  !
  CALL write_u_matrix (num_wann, ispin, 'umatrix'//TRIM(extension), u_matrix(:,:,1))
  call write_orbs (c0_tmp, num_wann, ispin, extension)
  !
  write(stdout, *) " **** Computing overlap between the Wannier orbs **** " 
  WRITE(stdout,'(5x,"[MEM] **Memory** needs per cpu")')
  !
  WRITE( stdout, '(5x,"[MEM] Wannier orbs in R space ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
  DBLE(16*nbnd*dfftp%nnr)/DBLE(1024*1024), dfftp%nnr, nbnd
  WRITE( stdout, '(5x,"[MEM] Wannier orbs in G space ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
  DBLE(16*nbnd*npw)/DBLE(1024*1024), npw, nbnd
  WRITE( stdout, '(5x,"[MEM] Overlap matrix ",f10.4," Mb", 5x,"(",i7,",",i5,")")') &
  DBLE(2*8*nbnd*nbnd)/DBLE(1024*1024), nbnd, nbnd
  WRITE( stdout, '(5x,"[MEM] Allocated aux vars ",f10.2," Mb", 5x,"(",i7,",",i5,")")') &
  DBLE(16*2*dfftp%nnr)/DBLE(1024*1024), dfftp%nnr, 2
  !
  ALLOCATE (psic1(dfftp%nnr), psic2(dfftp%nnr), c0_in_r(dfftp%nnr,num_wann))
  ALLOCATE (overlap_orbs(num_wann, num_wann))
  !
  npw = ngk(1)  
  ! 
  overlap_orbs = 0.0_DP
  IF (gamma_only) THEN
     !
     c0_in_r(:,:) = ( 0.D0, 0.D0 )
     !
     DO ibnd = 1, num_wann, 2
        !
        psic1(:) = ( 0.D0, 0.D0 )
        !
        IF ((ibnd < num_wann).OR.(num_wann < nbnd)) THEN
           ! 
           psic1(nls(1:npw))  = c0_tmp(1:npw,ibnd) + &
                            ( 0.D0, 1.D0 ) * c0_tmp(1:npw,ibnd+1)
           psic1(nlsm(1:npw)) = CONJG( c0_tmp(1:npw,ibnd) - &
                            ( 0.D0, 1.D0 ) * c0_tmp(1:npw,ibnd+1))
        ELSE
           ! 
           psic1(nls (1:npw)) = c0_tmp(1:npw,ibnd)
           psic1(nlsm(1:npw)) = CONJG( c0_tmp(1:npw,ibnd) )
           !     
        ENDIF
        !
        CALL invfft ('Wave', psic1, dfftp)
        !
        c0_in_r(:,ibnd) = psic1(:)
        !   
     ENDDO
     !
     DO ibnd = 1, num_wann, 2
        !
        psic1(:) = c0_in_r(:,ibnd) 
        !
        summ_ib_r = 0.0_DP
        summ_ib_i = 0.0_DP
        !
        DO ir = 1, dfftp%nnr
           !
           summ_ib_r = summ_ib_r +  (DBLE(psic1(ir)))**4 
           summ_ib_i = summ_ib_i +  (AIMAG(psic1(ir)))**4 
           !
        ENDDO
        !
        summ_ib_r = summ_ib_r / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
        summ_ib_i = summ_ib_i / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
        ! 
        CALL mp_sum(summ_ib_r, intra_image_comm)
        CALL mp_sum(summ_ib_i, intra_image_comm)
        !
        DO jbnd = 1, num_wann, 2 
           !
           psic2(:) = c0_in_r(:,jbnd) 
           !
           summ_jb_r = 0.0_DP
           summ_jb_i = 0.0_DP
           !
           DO ir = 1, dfftp%nnr
              !
              summ_jb_r = summ_jb_r +  (DBLE(psic2(ir)))**4
              summ_jb_i = summ_jb_i +  (AIMAG(psic2(ir)))**4
              !
           ENDDO
           !
           summ_jb_r = summ_jb_r / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           summ_jb_i = summ_jb_i / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           ! 
           CALL mp_sum(summ_jb_r, intra_image_comm)
           CALL mp_sum(summ_jb_i, intra_image_comm)
           !
           summ1 = 0.0_DP
           summ2 = 0.0_DP
           summ3 = 0.0_DP
           summ4 = 0.0_DP
           !
           DO ir = 1, dfftp%nnr
              !
              summ1  = summ1  + ((DBLE(psic1(ir)))**2 * (DBLE(psic2(ir)) )**2)
              summ2  =  summ2 + ((DBLE(psic1(ir)))**2 * (AIMAG(psic2(ir)) )**2)
              summ3 = summ3 + ((AIMAG(psic1(ir)))**2 * (DBLE(psic2(ir)) )**2)
              summ4 = summ4 + ((AIMAG(psic1(ir)))**2 * (AIMAG(psic2(ir)) )**2)
              !
           ENDDO
           !
           summ1 = summ1 / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           summ2 = summ2 / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           summ3 = summ3 / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           summ4 = summ4 / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           ! 
           CALL mp_sum(summ1, intra_image_comm)
           CALL mp_sum(summ2, intra_image_comm)
           CALL mp_sum(summ3, intra_image_comm)
           CALL mp_sum(summ4, intra_image_comm)
           !
           !overlap_orbs(ibnd,jbnd) = summ1
           overlap_orbs(ibnd,jbnd) = summ1/sqrt(summ_ib_r*summ_jb_r)
           !
           IF (jbnd < num_wann) THEN 
              !
              !overlap_orbs(ibnd,jbnd+1) = summ2
              overlap_orbs(ibnd,jbnd+1) = summ2/sqrt(summ_ib_r*summ_jb_i)
              !
           ENDIF
           !
           IF (ibnd < num_wann) THEN
              !
              !overlap_orbs(ibnd+1,jbnd) = summ3
              overlap_orbs(ibnd+1,jbnd) = summ3/sqrt(summ_ib_i*summ_jb_r)
              !
           ENDIF
           !
           IF ((ibnd < num_wann).AND.(jbnd < num_wann)) THEN
              !
              !overlap_orbs(ibnd+1,jbnd+1) = summ4
              overlap_orbs(ibnd+1,jbnd+1) = summ4/sqrt(summ_ib_i*summ_jb_i)
              !
           ENDIF
           !
        ENDDO 
        !
     ENDDO 
     !
  ELSE
     !
     DO ibnd = 1, num_wann
        !
        psic1(:) = ( 0.D0, 0.D0 )
        !
        psic1(nls(1:npw)) = c0_tmp(1:npw,ibnd)
        !
        CALL invfft ('Wave', psic1, dffts)
        !
        DO jbnd = 1, num_wann
           !
           psic2(:) = ( 0.D0, 0.D0 )
           !
           psic2(nls(1:npw))  = c0_tmp(1:npw,jbnd)
           !
           CALL invfft ('Wave', psic2, dffts)
           !
           summ1 = 0.0_DP
           DO ir = 1, dfftp%nnr
              summ1 = summ1 + (((DBLE(psic1(ir)))**2 + (AIMAG(psic1(ir)) )**2)) &
                            * (((DBLE(psic2(ir)))**2 + (AIMAG(psic2(ir)) )**2))  
           ENDDO
           !
           summ1 = summ1 / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
           ! 
           CALL mp_sum(summ1, intra_image_comm)
           !   
           overlap_orbs(ibnd,jbnd) = summ1
           !
        ENDDO
        !
     ENDDO 
     !
  ENDIF  
  !
  WRITE( stdout, '(5x," MAX and MAX Overlap values orb_i and the rest")')
  !
  DO ibnd = 1, num_wann
     !
     max_value_ibnd = maxval(overlap_orbs(ibnd,:))
     !
     WRITE( stdout,*), ibnd, maxval(overlap_orbs(ibnd,:)), minval(overlap_orbs(ibnd,:)) 
     !
  ENDDO
  !
  WRITE( stdout, '(5x," Relative Overlap values orb_1 and the rest")')
  !
  DO ibnd =1, num_wann
     write(stdout,*) 1, ibnd, overlap_orbs(1,ibnd)
  ENDDO
  !
  CALL write_ovl_matrix (num_wann, ispin,'overlap_matrix'//TRIM(extension), overlap_orbs)
  !call write_overlapmatrix (num_wann, ispin, extension, overlap_orbs)
  !
  WRITE( stdout,*) "--------------------------------------"
  WRITE( stdout,*) "JOB DONE"
  WRITE( stdout,*) "--------------------------------------"
  !
  deallocate ( psic1, psic2, c0_in_r )
  deallocate ( overlap_orbs )
  deallocate ( c0_tmp )
  !
  if (allocated(inc_band)) deallocate (inc_band)
  if (allocated(exclude_bands)) deallocate (exclude_bands)
  if (allocated(exclude_band)) deallocate (exclude_band)
  if (allocated(lwindow))  deallocate (lwindow)
  if (allocated(ndimwin)) deallocate (ndimwin)
  if (allocated(u_matrix_opt)) deallocate(u_matrix_opt)
  if (allocated(u_matrix)) deallocate(u_matrix)
  if (allocated(wannier_centres)) deallocate(wannier_centres)
  if (allocated(wannier_spreads)) deallocate(wannier_spreads)
  !
end subroutine wan2odd
!
!
subroutine write_orbs (c_emp, nempty, ispin, extension)
  !
  ! ...   This subroutine writes orbs to unit emptyunitc0
  ! 
  USE kinds,              ONLY: DP
  USE mp,                 ONLY: mp_sum
  USE io_global,          ONLY: ionode, ionode_id, stdout
  USE gvect,              ONLY: ig_l2g
  USE klist,              ONLY: ngk, igk_k
  USE wvfct,              ONLY: npw, npwx
  USE mp_images,          ONLY: intra_image_comm 
  USE mp_wave,            ONLY: mergewf
  USE mp,                 ONLY: mp_get, mp_size, mp_rank, mp_sum
  USE io_files,           ONLY: prefix, tmp_dir
  !
  IMPLICIT NONE
  ! 
  complex (DP),       intent(in) :: c_emp(npw, nempty)
  integer,            intent(in) :: nempty, ispin
  character(len=*),   intent(in) :: extension
  !
  LOGICAL :: exst, ierr
  INTEGER :: ig, i, ngw_g, iss, ngw_l, emptyunit
  INTEGER :: io_in_parent, nproc_in_parent, me_in_parent
  COMPLEX(DP), ALLOCATABLE :: ctmp(:)
  INTEGER, ALLOCATABLE     :: igk_l2g(:)
  character(len=256) :: fileempty, dirname
  character(len=4)   :: my_spin
  !
  ! ... Subroutine Body
  !
  ALLOCATE ( igk_l2g( npw ) )
  !
  ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
  !
  igk_l2g = 0
  DO ig = 1, ngk(1)
     igk_l2g(ig) = ig_l2g(igk_k(ig,1))
  ENDDO
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//'/wannier_orbs'//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing wannier wfcs to file :', fileempty
  !
  me_in_parent    = mp_rank( intra_image_comm )
  nproc_in_parent = mp_size( intra_image_comm )
  !    
  io_in_parent = 0
  IF( ionode ) io_in_parent = me_in_parent
  CALL mp_sum( io_in_parent, intra_image_comm )
  !
  emptyunit = 100
  !
  ngw_g    = npw
  ngw_l    = npw
  !
  CALL mp_sum( ngw_g, intra_image_comm )
  !
  ALLOCATE( ctmp( ngw_g ) )
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
     !
     REWIND( emptyunit )
     !
     WRITE ( emptyunit ) ngw_g, nempty
     !
  ENDIF
  !
  DO i = 1, nempty
     !
     ctmp = (0.0d0, 0.0d0)
     !
     CALL mergewf ( c_emp(:,i), ctmp(:), ngw_l, igk_l2g, & 
                    me_in_parent, nproc_in_parent, io_in_parent, intra_image_comm )
     !
     IF ( ionode ) THEN
        ! 
        WRITE ( emptyunit ) ( ctmp(ig), ig=1, ngw_g )
        !
     ENDIF
     !
  ENDDO
  ! 
  IF ( ionode ) CLOSE ( emptyunit )
  !
  DEALLOCATE(ctmp, igk_l2g)
  !
  RETURN
  ! 
ENDSUBROUTINE write_orbs 
!
!
subroutine write_umatrix (num_wan, ispin, extension, u_matrix)
  !
  ! ...   This subroutine writes orbs to unit emptyunitc0
  ! 
  USE kinds,              ONLY: DP
  USE io_global,          ONLY: ionode, stdout
  USE io_files,           ONLY: prefix, tmp_dir
  !
  IMPLICIT NONE
  ! 
  integer,           intent(in) :: num_wan,ispin
  real(DP),          intent(in) :: u_matrix(num_wan,num_wan)
  character(len=*),  intent(in) :: extension
  !
  INTEGER :: iw, jw, emptyunit
  character(len=256) :: fileempty, dirname
  character(len=4)   :: my_spin
  !
  ! ... Subroutine Body
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//'/wannier_umatrix'//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  !
  WRITE(stdout,*) 'Writing wannier_umatrix to file :', fileempty
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
     !
     REWIND( emptyunit )
     !
     WRITE ( emptyunit ) (( u_matrix(iw,jw), iw=1, num_wan), jw=1, num_wan ) 
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  RETURN
  ! 
ENDSUBROUTINE write_umatrix

subroutine write_overlapmatrix (num_wan,ispin,extension,u_matrix)
  !
  ! ...   This subroutine writes orbs to unit emptyunitc0
  ! 
  USE kinds,              ONLY: DP
  USE io_global,          ONLY: ionode, stdout
  USE io_files,           ONLY: prefix, tmp_dir
  !
  IMPLICIT NONE
  ! 
  integer,           intent(in) :: num_wan,ispin
  real(DP),          intent(in) :: u_matrix(num_wan,num_wan)
  character(len=*),  intent(in) :: extension
  !
  INTEGER :: iw, jw, emptyunit
  character(len=256) :: fileempty, dirname
  character(len=4)   :: my_spin
  !
  ! ... Subroutine Body
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//'/overlap_matrix'//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing overlap_matrix to file :', fileempty
  !
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
     !
     REWIND( emptyunit )
     !
     WRITE ( emptyunit ) (( u_matrix(iw,jw), iw=1, num_wan), jw=1, num_wan )
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//'/overlap_matrix_formated'//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing formatted overlap_matrix to file :', fileempty
  !
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'FORMATTED' )
     !
     REWIND( emptyunit )
     !
     DO iw = 1, num_wan
        DO jw = 1, num_wan 
           WRITE ( emptyunit, *) iw, jw, u_matrix(iw,jw)
        ENDDO
     ENDDO
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  RETURN
  ! 
ENDSUBROUTINE

subroutine wan2odd_combine ( nbnd_computed_occ, nbnd_computed_emp, ispin, extension)
  !
  use io_global,            only : stdout
  use kinds,                only : DP
  use io_files,             only : iunwfc, nwordwfc
  use wavefunctions_module, only : evc
  use wvfct,                only : nbnd, npw
  !
  IMPLICIT NONE 
  !
  integer, intent (in) ::  nbnd_computed_occ
  integer, intent (in) ::  nbnd_computed_emp
  integer, intent (in) ::  ispin
  character(len=*),  intent (in) :: extension
  !
  ! local variables
  ! 
  integer :: nbnd_to_save
  logical :: save_ks_wfc,  save_wn_wfc
  complex (DP), allocatable:: c0_tmp(:,:)
  ! 
  ! main program
  !
  save_ks_wfc = .false.
  save_wn_wfc = .false.
  SELECT CASE( trim( extension ) )
  CASE( '.ks.occ.emp' )
      save_ks_wfc = .true.
  CASE( '.wan.occ.emp' )
      save_wn_wfc = .true.
  CASE DEFAULT
        !
        call errore('wan2odd_combine',' extension should be ks.occ.emp or wan.occ.emp', 1)
        !
  END SELECT
  !
  nbnd_to_save = nbnd_computed_occ + nbnd_computed_emp
  !
  allocate (c0_tmp(npw, nbnd_to_save))
  c0_tmp(:,:) = (0.0_DP, 0.0_DP)
  !
  if (save_wn_wfc) then
     !
     CALL read_orbs_spin ( nbnd_computed_occ, npw, c0_tmp(:,1:nbnd_computed_occ), ispin, '.wan.occ' ) 
     !
     CALL read_orbs_spin ( nbnd_computed_emp, npw, c0_tmp(:,nbnd_computed_occ+1:nbnd_to_save), ispin, '.wan.emp' ) 
     ! 
  endif
  !
  if (save_ks_wfc) then
     !
     if (nbnd_to_save > nbnd) then
        !
        call errore('wan2odd_combine',' nbnd_computed_empty + nbnd_computed_occ > nbnd',1)
        !
     endif     
     !
     call davcio (evc, nwordwfc, iunwfc, 1, -1)
     !
     c0_tmp(:,1:nbnd_to_save) = evc(:,1:nbnd_to_save)
     ! 
  endif
  !
  call write_orbs (c0_tmp, nbnd_to_save, ispin, extension)
  !
  deallocate ( c0_tmp )
  !
  return
  !
end subroutine wan2odd_combine
!
!
!
subroutine write_u_matrix (num_wan,ispin,extension,u_matrix)
  !
  ! ...   This subroutine writes orbs to unit emptyunitc0
  ! 
  USE kinds,              ONLY: DP
  USE io_global,          ONLY: ionode, stdout
  USE io_files,           ONLY: prefix, tmp_dir 
  !
  IMPLICIT NONE
  ! 
  integer,           intent(in) :: num_wan,ispin
  complex(DP),       intent(in) :: u_matrix(num_wan,num_wan)
  character(len=*),  intent(in) :: extension
  !
  INTEGER :: iw, jw, emptyunit
  character(len=256) :: fileempty, dirname
  character(len=4)   :: my_spin
  !
  ! ... Subroutine Body
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing matrix to file :', fileempty
  !
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
     !
     REWIND( emptyunit )
     !
     WRITE ( emptyunit ) (( u_matrix(iw,jw), iw=1, num_wan), jw=1, num_wan )
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//TRIM(extension)//'_formated.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing formatted matrix to file :', fileempty
  !
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'FORMATTED' )
     !
     REWIND( emptyunit )
     !
     DO iw = 1, num_wan
        DO jw = 1, num_wan 
           WRITE ( emptyunit, *) iw, jw, u_matrix(iw,jw)
        ENDDO
     ENDDO
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  RETURN
  ! 
ENDSUBROUTINE

subroutine write_ovl_matrix (num_wan,ispin,extension,u_matrix)
  !
  ! ...   This subroutine writes orbs to unit emptyunitc0
  ! 
  USE kinds,              ONLY: DP
  USE io_global,          ONLY: ionode, stdout
  USE io_files,           ONLY: prefix, tmp_dir 
  !
  IMPLICIT NONE
  ! 
  integer,           intent(in) :: num_wan,ispin
  real(DP),          intent(in) :: u_matrix(num_wan,num_wan)
  character(len=*),  intent(in) :: extension
  !
  INTEGER :: iw, jw, emptyunit
  character(len=256) :: fileempty, dirname
  character(len=4)   :: my_spin
  !
  ! ... Subroutine Body
  !
  dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/'
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//TRIM(extension)//'.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing matrix to file :', fileempty
  !
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'UNFORMATTED' )
     !
     REWIND( emptyunit )
     !
     WRITE ( emptyunit ) (( u_matrix(iw,jw), iw=1, num_wan), jw=1, num_wan )
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  WRITE(my_spin,'(i1)') ispin
  fileempty = TRIM(dirname)//TRIM(extension)//'_formated.'//TRIM(my_spin)//'.dat'
  WRITE(stdout,*) 'Writing formatted matrix to file :', fileempty
  !
  emptyunit = 100
  !
  IF ( ionode ) THEN
     ! 
     ! 
     OPEN( UNIT = emptyunit, FILE = TRIM(fileempty), status = 'unknown', FORM = 'FORMATTED' )
     !
     REWIND( emptyunit )
     !
     DO iw = 1, num_wan
        DO jw = 1, num_wan 
           WRITE ( emptyunit, *) iw, jw, u_matrix(iw,jw)
        ENDDO
     ENDDO
     !
  ENDIF
  !
  IF ( ionode ) CLOSE ( emptyunit )
  !
  RETURN
  ! 
ENDSUBROUTINE

