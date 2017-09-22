SUBROUTINE kcp_read_init_orbs ( nbndx, c0 )
  ! 
  !  this routine read the occupied  wannier wave functions from pw2wannier PP code
  ! 
      USE kinds,              ONLY: DP
      USE mp_global,          ONLY: me_image, nproc_image, intra_image_comm
      USE io_global,          ONLY: stdout, ionode, ionode_id
      USE mp,                 ONLY: mp_bcast, mp_sum
      USE mp_wave,            ONLY: splitwf
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: ig_l2g
      !
      IMPLICIT none
      ! 
      INTEGER,     INTENT(IN)  :: nbndx
      COMPLEX(DP), INTENT(OUT) :: c0(ngw,nbndx)
      ! 
      LOGICAL :: exst
      INTEGER :: ierr, ig, ibnd, iss
      INTEGER :: ngw_rd, nbnd_rd, ngw_l, ngw_g
      INTEGER :: emptyunitc0
      ! 
      CHARACTER(LEN=256) :: fileocc
      !
      COMPLEX(DP), ALLOCATABLE :: ctmp(:)
      !
      ! ... Subroutine Body
      !
      emptyunitc0 = 200
      !
      ngw_g    = ngw
      ngw_l    = ngw
      !
      CALL mp_sum( ngw_g, intra_image_comm )
      !
      ALLOCATE( ctmp(ngw_g) )
      !
      fileocc = 'evc0_occupied.dat'
      !
      IF ( ionode ) THEN
         !
         INQUIRE( FILE = TRIM(fileocc), EXIST = exst )
         !
         IF ( exst ) THEN
            !
            OPEN( UNIT=emptyunitc0, FILE=TRIM(fileocc), STATUS='OLD', FORM='UNFORMATTED' )
            !
            READ(emptyunitc0) ngw_rd, nbnd_rd
            !
            IF ( ngw_g .ne. ngw_rd ) THEN
               !
               exst = .false.
               WRITE( stdout,10)  TRIM(fileocc) 
               WRITE( stdout,20)  ngw_g, ngw_rd
               !
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
 10   FORMAT('*** OCCUPIED STATES : wavefunctions dimensions changed  ', A )
 20   FORMAT('*** NGW_G = ', I8, ' NE_READ = ', I4)
      !  
      CALL mp_bcast(exst,   ionode_id, intra_image_comm)
      CALL mp_bcast(nbnd_rd,  ionode_id, intra_image_comm)
      CALL mp_bcast(ngw_rd, ionode_id, intra_image_comm)
      !
      c0(:,:) = (0.0_DP, 0.0_DP)
      !
      IF ( exst ) THEN
         ! 
         DO ibnd = 1, MIN( nbndx, nbnd_rd )
            !
            IF ( ionode ) THEN
               ! 
               READ(emptyunitc0) ( ctmp(ig), ig = 1, MIN( SIZE(ctmp), ngw_rd ) )
               !
            ENDIF
            ! 
            IF ( ibnd <= nbnd_rd ) THEN
               !
               CALL splitwf(c0(:,ibnd), ctmp, ngw_l, ig_l2g, me_image, &
                            nproc_image, ionode_id, intra_image_comm)
               !
            ENDIF
            !
         ENDDO
         !
      ELSE
         !
         IF (.not. exst ) CALL errore( 'wave_init_wannier_pwscf', 'Something wrong with reading evc_occupied file', 1 )
         !
      ENDIF
      ! 
      IF ( ionode .AND. exst ) THEN
         !
         CLOSE(emptyunitc0) 
         !
      ENDIF
      ! 
      DEALLOCATE(ctmp)
      !
      RETURN
      !
END SUBROUTINE
!
!
!
SUBROUTINE read_odd_nkscalfact(nbndx, odd_alpha)
  !
  USE kinds,         ONLY : dp
  USE io_global,     ONLY : stdout, ionode, ionode_id 
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE mp_world,      ONLY : world_comm
  !
  IMPLICIT NONE
  !
  integer,  intent(in)  :: nbndx 
  real(dp), intent(out) :: odd_alpha(nbndx)
  ! 
  integer :: ibnd, nbnd_input
  !
  ! read eigenvalues from file
  !
  call mp_barrier(world_comm)
  !
  if (ionode) then 
     !
     open(unit = 99,file ='odd_nkscalfac.dat',form = 'formatted',status = 'old')
     !
     read(99,*) nbnd_input
     ! 
  endif 
  ! 
  call mp_bcast (nbnd_input,ionode_id,world_comm )
  !
  if (nbnd_input > nbndx) then
     ! 
     call errore ('read_odd_nkscalfact', 'number nbnd_input > nbndx = ', nbndx) 
     !
  endif
  ! 
  if (ionode) then
     !  
     do ibnd = 1, nbnd_input
        ! 
        read (99, * ) odd_alpha(ibnd)
        !
     enddo
     !
     close (99)
     !     
  endif
  !
  call mp_bcast(odd_alpha,ionode_id,world_comm)
  !
  return
  !
ENDSUBROUTINE
