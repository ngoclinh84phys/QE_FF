!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE kcp_electrons_module
   !
   USE kinds,              ONLY : DP
   USE electrons_base,     ONLY : nspin, nupdwn, nudx
   !
   IMPLICIT NONE    
   ! 
   SAVE
   !
   PRIVATE 
   !
   PUBLIC :: kcp_electrons_setup   ! NLN
   PUBLIC :: kcp_deallocate_electrons 
   PUBLIC :: kcp_print_eigenvalues ! NLN
   PUBLIC :: print_centers_spreads 
   !
   ! ... Spread of orbitals 
   !
   REAL(DP), ALLOCATABLE :: wfc_centers(:,:,:) !added:giovanni wfc_centers
   REAL(DP), ALLOCATABLE :: wfc_spreads(:,:,:) !added:giovanni wfc_spreads
   LOGICAL :: icompute_spread = .TRUE.
   PUBLIC  :: wfc_centers, wfc_spreads, icompute_spread
   !
   CONTAINS
   ! 
   SUBROUTINE kcp_electrons_setup ( )
     !
     IMPLICIT NONE
     !
     INTEGER :: ierr, i
     ! 
     ! KOOPMANS COMPLIANT FUNCTIONALS
     !  
     IF( ALLOCATED( wfc_centers ) ) DEALLOCATE( wfc_centers )
     IF(nudx > 0) THEN
        ALLOCATE( wfc_centers( 4, nudx, nspin ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_centers ',ierr)
        wfc_centers = 0.0_DP
     ENDIF
     !
     IF( ALLOCATED( wfc_spreads ) ) DEALLOCATE( wfc_spreads )
     IF(nudx > 0) THEN
        ALLOCATE( wfc_spreads( nudx, nspin, 2 ), STAT=ierr)
        IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_spreads ',ierr)
        wfc_spreads = 0.0_DP
     ENDIF
     !
     ! KOOPMANS COMPLIANT FUNCTIONALS
     !  
     RETURN
     !
   ENDSUBROUTINE kcp_electrons_setup

   !
   SUBROUTINE kcp_print_eigenvalues( ei_unit, tfile, tstdout, nfi, tps )
      !
      USE constants,        ONLY : autoev 
      USE io_global,        ONLY : stdout, ionode
      USE electrons_module, ONLY : ei 
      !
      INTEGER,  INTENT(IN) :: ei_unit
      LOGICAL,  INTENT(IN) :: tfile, tstdout
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps
      !
      INTEGER :: i, j, ik
      !
      IF ( tfile ) THEN
          WRITE(ei_unit,30) nfi, tps
      END IF
      !
      ik = 1
      !
      ! KOOPMANS COMPLIANT FUNCTIONALS, MESS: TO PRINT HOMO EIGENVALUE  
      IF (nspin==1) THEN
         WRITE( stdout,1101)
         WRITE( stdout, 1444) MAXVAL(ei(1:nupdwn(1),1)*autoev, nupdwn(1))
      ELSE
         WRITE( stdout,1101)
         WRITE( stdout, 1444) MAX(MAXVAL(ei(1:nupdwn(1),1)*autoev, nupdwn(1)), MAXVAL(ei(1:nupdwn(2),2)*autoev, nupdwn(2)))
      ENDIF
      !
      ! KOOPMANS COMPLIANT FUNCTIONALS, MESS: TO PRINT HOMO EIGENVALUE  
      ! 
      DO j = 1, nspin
         !
         IF( tstdout ) THEN
            WRITE( stdout,1002) ik, j
            WRITE( stdout,1004) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         END IF
         !
         IF( tfile ) THEN
            WRITE(ei_unit,1010) ik, j
            WRITE(ei_unit,1020) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         END IF
         !
      END DO
      !
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
 1002 FORMAT(/,3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1004 FORMAT(10F8.2)
 1010 FORMAT(3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1020 FORMAT(10F8.2)
 1101 FORMAT(/,3X,'HOMO Eigenvalue (eV)',/) !added_giovanni
 1444 FORMAT(1F10.4)
      !
      RETURN
      !
   ENDSUBROUTINE kcp_print_eigenvalues
   ! 
      ! KOOPMANS COMPLIANT FUNCTIONAL
   SUBROUTINE print_centers_spreads( ei_unit, tfile, tstdout, nfi, tps )
      !
      use constants,  only : autoev, hartree_si, electronvolt_si
      USE io_global,  ONLY : stdout, ionode
      use nksic,      ONLY : pink, do_orbdep, odd_alpha 
      !
      INTEGER,  INTENT(IN) :: ei_unit
      LOGICAL,  INTENT(IN) :: tfile, tstdout
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps
      !
      INTEGER :: i, j, ik
      !
      ik = 1
      !
      DO j = 1, nspin
         !
         WRITE(stdout,1222) ik, j
         !
         IF (.not. do_orbdep) THEN
            !       
            WRITE(stdout,1442) ( i, wfc_centers(1:4, i, j), &
                                    wfc_spreads( i, j, 1 ), & 
                                    wfc_spreads( i, j, 2 ), &
                                    i = 1, nupdwn(j) )
         ELSE
            !  
            WRITE(stdout,1444) ( i, wfc_centers(1:4, i, j), &
                                    wfc_spreads( i, j, 1),  &
                                    wfc_spreads( i, j, 2),  &
                                    pink(i)*hartree_si/electronvolt_si, & ! Linh i index is wrong
                                    odd_alpha(i), i = 1, nupdwn(j) )
         ENDIF 
         !
      ENDDO
      !
 1222 FORMAT(/,3X,'Orb -- Charge  ---   Centers xyz (Bohr)  ---  Spreads (Bohr^2) - SH(eV), kp = ',I3, ' , spin = ',I2,/)
 1442 FORMAT('OCC', I5, ' --',F8.2,'   ---',3F8.2,'   ---',2F8.3  )
 1444 FORMAT('OCC', I5, ' --',F8.2,'   ---',3F8.2,'   ---',4F8.3  )
      ! 
      RETURN
      ! 
   END SUBROUTINE print_centers_spreads 
   !
   SUBROUTINE kcp_deallocate_electrons
      !
      INTEGER :: ierr
      !
      ! KOOPMANS COMPLIANT FUNCTIONAL
      !
      IF (ALLOCATED(wfc_centers)) THEN
         !
         DEALLOCATE(wfc_centers, STAT=ierr)
         IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_centers ',ierr )
         !
      ENDIF
      !
      IF (ALLOCATED(wfc_spreads)) THEN 
         ! 
         DEALLOCATE(wfc_spreads, STAT=ierr)
         IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_spreads ',ierr )
         !
      ENDIF
      !
      RETURN
      !
   ENDSUBROUTINE kcp_deallocate_electrons 
   !     
END MODULE kcp_electrons_module
