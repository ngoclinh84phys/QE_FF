
subroutine init_nksic()
   !
   ! this routine is called anyway, even if do_nk=F
   !
   use nksic,            only : do_orbdep, do_nk, do_nkipz, do_nkpz, do_pz, &
                                do_nki, do_bare_eigs, do_pz_renorm, kfact, &
                                do_wref, do_wxd, fref, rhobarfact, &
                                vanishing_rho_w, &
                                do_spinsym, f_cutoff, &
                                nkscalfact, nksic_memusage, allocate_nksic, odd_alpha
   !
   use nksic,            only : esic_conv_thr,do_innerloop,do_innerloop_empty, do_innerloop_cg, &
                                innerloop_dd_nstep, innerloop_cg_nsd, innerloop_cg_nreset, innerloop_nmax, &
                                innerloop_cg_ratio, innerloop_init_n, innerloop_until, &
                                innerloop_atleast
   !
   use input_parameters, only : which_orbdep, &
                                do_innerloop_ => do_innerloop_,&
                                nkscalfact_ => nkscalfact,     &
                                innerloop_step, nkscalfact_odd & 
   !
   use io_global,        only : meta_ionode, stdout
   use electrons_base,   only : nspin, nbspx
   use gvecw,            only : ngw
   use fft_base,         only : dffts, dfftp
   use ions_base,        only : nat
   !
   implicit none
   !
   logical       :: found, do_hybrid=.FALSE.
   integer       :: i
   character(10) :: subname='nksic_init'
   character(1), external :: lowercase
   !
   ! overwriten by which_orbdep, if not empty
   !
   f_cutoff      = 0.01
   !
   nkscalfact          = nkscalfact_
   do_innerloop        = do_innerloop_
   innerloop_nmax      = innerloop_step
   odd_alpha           = nkscalfact_odd 
   !
   ! Optimal defaults
   ! 
   esic_conv_thr       = 1.0e-4*nkscalfact_
   innerloop_dd_nstep  = .TRUE.
   innerloop_cg_nsd    = innerloop_cg_nsd_
   innerloop_cg_nreset = innerloop_cg_nreset_
   innerloop_nmax      = innerloop_nmax_
   innerloop_init_n    = innerloop_init_n_
   innerloop_atleast   = innerloop_atleast_
   innerloop_cg_ratio  = innerloop_cg_ratio_
   innerloop_until     = innerloop_until_
   !
  
   !
   ! use the collective var which_orbdep
   !
   do i = 1, LEN_TRIM( which_orbdep_ )
      which_orbdep_(i:i) = lowercase( which_orbdep_(i:i) )
   enddo
   !
   do_orbdep = .true.
   !   
   SELECT CASE ( TRIM(which_orbdep_) )
      !
      CASE ( "", "none" )
         !  
         do_orbdep = .false.
         ! do nothing
      CASE ( "nki" )
         do_nki   = .TRUE.
         do_wxd   = .TRUE.
         fref     = 1.0
      CASE ( "pz", "perdew-zunger" )
         do_pz    = .TRUE.
      CASE ( "pzha")
         do_pzha  = .TRUE.
      CASE ( "nkipz", "pznki" )
         do_nkipz = .TRUE.
         do_wxd   = .TRUE.
         fref     = 1.0
      CASE DEFAULT
         call errore(subname,"invalid which_orbdep = "//TRIM(which_orbdep_),10)
   END SELECT
   !
   if ( do_pz .and. meta_ionode ) then
      !
      write(stdout,2001) do_pz
      !
   else if ( (do_nki .or. do_nkipz) .and. meta_ionode ) then
      !
      write(stdout,2002) do_nki
      write(stdout,2004) rhobarfact, nkscalfact
      ! 
   endif
   !
   ! read referece alpha from file, if any | linh
   ! wherein, the information of n_evc0_fixed, ref_alpha0,
   ! broadening of orbitals will be readed.   
   ! 
   ! call readfile_refalpha()
   ! 
   if ( do_orbdep ) call allocate_nksic( dfftp%nnr, ngw, nspin, nbspx, nat)
   !
   if ( do_orbdep ) odd_alpha(:) = 1.d0
   !
   if ( (do_nk .or. do_nkpz ) .and. meta_ionode ) then
       write(stdout,2010) do_wxd, do_wref, do_nkpz
   endif
   !
   if ( do_orbdep ) then
      !
      write(stdout,2005) vanishing_rho_w
      !  
      write( stdout, "(3x, 'NK memusage = ', f10.3, ' MB', /)" ) &
           nksic_memusage( )
      !    
   endif
   !
   if ( do_orbdep .and. innerloop_until<1) then
      !
      innerloop_until=-1
      !
   endif
   !
2001  format( 3X,'PZ sic = ',l4 ,/)
2002  format( 3X,'NK sic with integral ref = ',l4, / )
2004  format( 3X,'NK background density factor = ',f7.4, /, &
              3X,'NK scaling factor = ',f7.4 )
2005  format( 3X,'rhothr = ',e8.1 )
2006  format( 3X,'f_cutoff = ',f7.4, /, &
              3X,'do_spinsym   = ',l4 )
2010  format( 3X,'NK cross-derivatives = ',l4, /, &
              3X,'NK reference derivatives = ',l4, /, &
              3X,'NK on top of PZ = ',l4 )
2030  format( 3X,'NK applied up to orbital',i7)
   !
   return
   ! 
end subroutine init_nksic
