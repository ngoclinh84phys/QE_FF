
subroutine nksic_init ()
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
   use input_parameters, only : which_orbdep_ => which_orbdep, &
                                fref_ => fref, &
                                rhobarfact_ => rhobarfact, &
                                vanishing_rho_w_ => vanishing_rho_w, &
                                nkscalfact_ => nkscalfact
   !
   use input_parameters, only : esic_conv_thr_ => esic_conv_thr, &
                                do_innerloop_ => do_innerloop, &
                                do_innerloop_empty_ => do_innerloop_empty, &
                                do_innerloop_cg_ => do_innerloop_cg, &
                                innerloop_dd_nstep_ => innerloop_dd_nstep, &
                                innerloop_cg_nsd_ => innerloop_cg_nsd, &
                                innerloop_cg_nreset_ => innerloop_cg_nreset, &
                                innerloop_nmax_ => innerloop_nmax, &
                                innerloop_init_n_ => innerloop_init_n, &
                                innerloop_atleast_ => innerloop_atleast, &
                                innerloop_cg_ratio_ => innerloop_cg_ratio, &
                                innerloop_until_ => innerloop_until
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
   esic_conv_thr       = esic_conv_thr_
   do_innerloop        = do_innerloop_
   do_innerloop_empty  = do_innerloop_empty_
   do_innerloop_cg     = do_innerloop_cg_
   innerloop_dd_nstep  = innerloop_dd_nstep_
   innerloop_cg_nsd    = innerloop_cg_nsd_
   innerloop_cg_nreset = innerloop_cg_nreset_
   innerloop_nmax      = innerloop_nmax_
   innerloop_init_n    = innerloop_init_n_
   innerloop_atleast   = innerloop_atleast_
   innerloop_cg_ratio  = innerloop_cg_ratio_
   innerloop_until     = innerloop_until_
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
      CASE ( "nk", "non-koopmans" )
         do_nk   = .TRUE.
         do_wref = .TRUE.
         do_wxd  = .TRUE.
      CASE ( "nk0" )
         do_nk   = .TRUE.
         do_wref = .FALSE.
         do_wxd  = .FALSE.
      CASE ( "nki" )
         do_nki  = .TRUE.
         do_wxd  = .TRUE.
         fref    = 1.0
      CASE ( "pz", "perdew-zunger" )
         do_pz   = .TRUE.
      CASE ( "nkpz", "pznk" )
         do_nkpz = .TRUE.
      CASE ( "nkipz", "pznki" )
         do_nkipz = .TRUE.
         do_wxd  = .TRUE.
         fref    = 1.0
      CASE DEFAULT
         call errore(subname,"invalid which_orbdep = "//TRIM(which_orbdep_),10)
   END SELECT
   !
   found = .false.
   !
   if (do_orbdep) then
      !
      if ( do_nk.or.do_pz.or.do_nki.or.do_nkpz.or.do_nkipz.or.do_hybrid ) found=.true.
      !   
      if (.not. found ) CALL errore(subname,'no compatible orbital-dependent scheme specified',1)
      !
   endif
   !
   ! check only one orbital dependent scheme is used
   !
   found = .FALSE.
   !
   if (do_nk .and. (do_pz.or.do_nki.or.do_nkpz.or.do_nkipz)) found=.TRUE.
   if (do_nki.and. (do_pz.or.do_nk.or.do_nkpz.or.do_nkipz)) found=.TRUE.
   if (do_pz.and. ( do_nk .or. do_nki .or. do_nkpz .or. do_nkipz ) ) found=.TRUE.
   if (do_nkpz.and. ( do_nk .or. do_nki .or. do_pz   .or. do_nkipz ) ) found=.TRUE.
   if (do_nkipz.and. ( do_nk .or. do_nki .or. do_pz   .or. do_nkpz  ) ) found=.TRUE.
   !
   if ( found ) CALL errore(subname,'more than one orb-dependent schme used',1)
   !
   vanishing_rho_w = vanishing_rho_w_
   rhobarfact    = rhobarfact_
   nkscalfact    = nkscalfact_
   !
   if ( do_nki .and. fref /= 1.0 ) CALL errore(subname,'nki and fref /= 1.0 ',1)
   !
   if ( (do_nk .or. do_nkpz) .and. meta_ionode ) then
      write(stdout,2000) fref
      write(stdout,2004) rhobarfact, nkscalfact
   else if ( do_pz .and. meta_ionode ) then
      write(stdout,2001) do_pz
   else if ( (do_nki .or. do_nkipz) .and. meta_ionode ) then
      write(stdout,2002) do_nki
      write(stdout,2004) rhobarfact, nkscalfact
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
   if (do_orbdep) odd_alpha(:) = 1.d0
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
2000  format( 3X,'NK sic with reference occupation = ',f7.4, /)
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
end subroutine nksic_init
