!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE emptystates_electrons_module
!=----------------------------------------------------------------------------=!
        USE kinds
        USE dspev_module,       ONLY: pdspev_drv, dspev_drv
        USE electrons_base,     ONLY: nbnd, nbndx, nbsp, nbspx, nspin, nel, nelt, &
                                      nupdwn, iupdwn, telectrons_base_initval, f, &
                                      nudx, nupdwn_bgrp, iupdwn_bgrp, nudx_bgrp, &
                                      nbsp_bgrp, nbspx_bgrp, i2gupdwn_bgrp

        USE cp_electronic_mass, ONLY: ecutmass => emass_cutoff, emass, emass_precond


        IMPLICIT NONE
        SAVE

        PRIVATE 
        !
        INTEGER, PARAMETER :: nspinx  = 2
        !
        INTEGER :: n_emp               =  0  ! number of empty states
        INTEGER :: nudx_emp            =  0  ! maximum number of empty states per spin
        INTEGER :: nupdwn_emp(nspinx)  =  0  ! number of empty states
        INTEGER :: iupdwn_emp(nspinx)  =  0  ! number of empty states
        !
        INTEGER :: n_emp_l(nspinx)     =  0  ! local number of emptystates ( for each spin components )
        !
        INTEGER  :: max_emp = 0    !  maximum number of iterations for empty states
        !  
        REAL(DP), ALLOCATABLE :: ei_emp(:,:)
        REAL(DP), ALLOCATABLE :: wfc_centers_emp(:,:,:) 
        REAL(DP), ALLOCATABLE :: wfc_spreads_emp(:,:,:) 
        !
        PUBLIC :: n_emp, ei_emp, n_emp_l, max_emp
        PUBLIC :: nupdwn_emp, iupdwn_emp, nudx_emp
        PUBLIC :: wfc_centers_emp, wfc_spreads_emp
        !
        PUBLIC :: emptystates_electrons_setup
        PUBLIC :: emptystates_print_eigenvalues 
        PUBLIC :: emptystates_print_centers_spreads 
        PUBLIC :: emptystates_deallocate_electrons 
!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!
     !     
     SUBROUTINE emptystates_electrons_setup( )
        !
        use input_parameters, only : number_emptystates_ => number_emptystates, & 
                                     maxstep_minimize_emptystates_ => maxstep_minimize_emptystates
        !
        IMPLICIT NONE
        INTEGER :: ierr, i
        !
        max_emp = maxstep_minimize_emptystates_
        n_emp = number_emptystates_ 
        !
        nupdwn_emp(1) = n_emp
        iupdwn_emp(1) = 1
        nudx_emp = n_emp
        !
        IF ( nspin == 2 ) THEN
           nupdwn_emp(2) = n_emp
           iupdwn_emp(2) = 1 + n_emp
        ENDIF
        ! 
        IF ( ALLOCATED( ei_emp ) ) DEALLOCATE( ei_emp )
        IF ( n_emp > 0 ) THEN
           ALLOCATE( ei_emp( n_emp, nspin ), STAT=ierr)
           IF( ierr/=0 ) CALL errore( ' electrons ',' allocating ei_emp ',ierr)
           ei_emp = 0.0_DP
        ENDIF
        !  
        IF ( ALLOCATED( wfc_centers_emp ) ) DEALLOCATE( wfc_centers_emp )
        IF ( nudx_emp > 0 ) THEN
           ALLOCATE( wfc_centers_emp(4, nudx_emp,nspin ), STAT=ierr)
           IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_centers_emp',ierr)
           wfc_centers_emp = 0.0_DP
        ENDIF
        ! 
        IF ( ALLOCATED( wfc_spreads_emp ) ) DEALLOCATE( wfc_spreads_emp )
        IF ( nudx_emp > 0 ) THEN
           ALLOCATE( wfc_spreads_emp( nudx_emp, nspin, 2 ), STAT=ierr)
           IF( ierr/=0 ) CALL errore( ' electrons ',' allocating wfc_spreads_emp',ierr)
           wfc_spreads_emp = 0.0_DP
        ENDIF
        !
        RETURN
        !
     ENDSUBROUTINE emptystates_electrons_setup
     !
     SUBROUTINE emptystates_print_eigenvalues( )
        !
        use constants,  only : autoev 
        USE io_global,  ONLY : stdout, ionode
        USE electrons_module, ONLY : ei
        !
        INTEGER :: i, j, ik
        !
        ik = 1
        !
        IF (n_emp.gt.0) THEN
           ! 
           IF (nspin==1) THEN
              WRITE( stdout,1201)
              WRITE( stdout,1444) MINVAL(ei_emp(1:n_emp,1)*autoev, n_emp)
           ELSE
              WRITE( stdout,1201)
              WRITE( stdout,1444) MIN(MINVAL(ei_emp(1:n_emp,1)*autoev, n_emp), MINVAL(ei_emp(1:n_emp,2)*autoev, n_emp))
           ENDIF
           ! 
           DO j = 1, nspin
              !
              WRITE( stdout,1004) ik, j
              WRITE( stdout,1005) ( ei_emp( i, j ) * autoev , i = 1, n_emp )
              IF (nupdwn(j)>0) &
                 WRITE( stdout,1006) ( ei_emp( 1, j ) - ei( nupdwn(j), j ) ) * autoev
           ENDDO
           !
        ENDIF   
        !
 1201   FORMAT(/,3X,'LUMO Eigenvalue (eV)',/) !added_giovanni
 1444   FORMAT(1F10.4)
 1004   FORMAT(/,3X,'Empty States Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1005   FORMAT(10F8.2)
 1006   FORMAT(/,3X,'Electronic Gap (eV) = ',F10.4,/)
        !
        RETURN
        !
     ENDSUBROUTINE emptystates_print_eigenvalues
     !
     SUBROUTINE emptystates_print_centers_spreads( )
        !
        use constants,      only : autoev, hartree_si, electronvolt_si
        USE io_global,      only : stdout, ionode
        use nksic,          only : pink_emp, odd_alpha_emp
        !
        INTEGER :: i, j, ik
        !
        ik = 1
        !
        IF ( n_emp > 0 ) THEN
           !
           DO j = 1, nspin
              !
              WRITE( stdout,1333) ik, j
              !
              WRITE( stdout,1446) ( i, wfc_centers_emp(1:4, i, j ), & 
                                       wfc_spreads_emp( i, j, 1),   &
                                       wfc_spreads_emp( i, j, 2),   &
                                       pink_emp(i)*hartree_si/electronvolt_si, & ! Linh check i index, it's wrong
                                       odd_alpha_emp (i), & ! Linh check i index, it's wrong
                                       i = 1, nupdwn_emp(j)) 
              !  
           ENDDO
           !
        ENDIF
        !
 1333   FORMAT(/,3X,'Orb -- Empty Charge  ---   Centers xyz (Bohr)  ---  Spreads (Bohr^2) - SH(eV), kp = ',I3, ' , spin = ',I2,/)
 1446   FORMAT('EMP', I5, ' --',F8.2,'   ---',3F8.2,'   ---',4F8.3)
        !
     ENDSUBROUTINE emptystates_print_centers_spreads
     !
     SUBROUTINE emptystates_deallocate_electrons
        !
        INTEGER :: ierr
        !
        IF (ALLOCATED(ei_emp)) THEN
           DEALLOCATE(ei_emp, STAT=ierr)
           IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating ei_emp ',ierr )
        ENDIF
        !  
        IF (ALLOCATED(wfc_centers_emp)) THEN
           DEALLOCATE(wfc_centers_emp, STAT=ierr)
           IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_centers_emp ',ierr )
        ENDIF
        ! 
        IF (ALLOCATED(wfc_spreads_emp)) THEN
           DEALLOCATE(wfc_spreads_emp, STAT=ierr)
           IF( ierr/=0 ) CALL errore( ' deallocate_electrons ',' deallocating wfc_spreads_emp ',ierr )
        ENDIF
        !
        RETURN
        ! 
     ENDSUBROUTINE emptystates_deallocate_electrons
     !      
ENDMODULE emptystates_electrons_module
