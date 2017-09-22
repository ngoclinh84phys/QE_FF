
subroutine init_nksic()
   !
   ! this routine is called anyway, even if do_nk=F
   !
   use kinds,               only : dp
   use io_global,           only : stdout
   use nksic,               only : do_orbdep, do_nkipz, do_pz, do_pzha, do_nki, &
                                   do_wxd, fref, f_cutoff, &
                                   nkscalfact, odd_nkscalfact, &
                                   nksic_memusage, allocate_nksic, odd_alpha
   !
   use nksic,               only : esic_conv_thr, do_innerloop, &
                                   innerloop_cg_nsd, innerloop_cg_nreset, innerloop_nmax, &
                                   innerloop_cg_ratio, innerloop_init_n, innerloop_until, &
                                   innerloop_atleast
   !
   use input_parameters,    only : which_orbdep, &
                                   l_read_orbs_from_file, &
                                   do_innerloop_ => do_innerloop,&
                                   nkscalfact_ => nkscalfact,     &
                                   innerloop_step, nkscalfact_odd 
   !
   use electrons_base,      only : nspin, nbspx
   use gvecw,               only : ngw
   use fft_base,            only : dffts, dfftp
   use ions_base,           only : nat
   use wavefunctions_module,only : c0_bgrp
   !
   implicit none
   !
   logical       :: found, do_hybrid=.FALSE.
   integer       :: i
   character(10) :: subname='nksic_init'
   character(1), external :: lowercase
   !
   ! pass input kcp_vars to global kcp_vars 
   !
   nkscalfact          = nkscalfact_
   do_innerloop        = do_innerloop_
   innerloop_nmax      = innerloop_step
   odd_nkscalfact      = nkscalfact_odd 
   !
   ! read input orbs from file, if any
   !
   if ( l_read_orbs_from_file ) then
      !
      c0_bgrp(:,:) = (0.0_dp, 0.0_dp)
      !
      call kcp_read_init_orbs(nbspx, c0_bgrp)
      !
   endif
   !
   ! Optimal defaults, modify for testing only 
   ! 
   f_cutoff            = 0.1
   esic_conv_thr       = 1.0E-4*nkscalfact_
   innerloop_cg_nsd    = 20 
   innerloop_cg_nreset = 10 
   innerloop_init_n    = 10000
   innerloop_atleast   = 0
   innerloop_cg_ratio  = 1.d-3 
   innerloop_until     = -1
   !
   ! use the collective var which_orbdep
   !
   do i = 1, LEN_TRIM( which_orbdep )
      !
      which_orbdep(i:i) = lowercase( which_orbdep(i:i) )
      !
   enddo
   !
   do_orbdep = .true.
   !   
   SELECT CASE ( TRIM(which_orbdep) )
      !
      CASE ( "", "none" )
         !  
         do_orbdep = .false.
         ! 
      CASE ( "pz", "perdew-zunger" )
         !
         do_pz    = .TRUE.
         !
      CASE ( "pzha")
         !
         do_pzha  = .TRUE.
         ! 
      CASE ( "nki" )
         !
         do_nki   = .TRUE.
         do_wxd   = .TRUE.
         fref     = 1.0
         !
      CASE ( "nkipz", "pznki" )
         ! 
         do_nkipz = .TRUE.
         do_wxd   = .TRUE.
         fref     = 1.0
         !
      CASE DEFAULT
         !
         call errore(subname,"invalid which_orbdep = "//TRIM(which_orbdep),10)
         !
   END SELECT
   !
   if ( do_pz ) then
      !
      write(stdout,2001) do_pz
      !
   else if ( do_nki .or. do_nkipz ) then
      !
      write(stdout,2002) do_nki
      write(stdout,2004) nkscalfact
      ! 
   endif
   !
   if ( do_orbdep ) call allocate_nksic( dfftp%nnr, ngw, nspin, nbspx, nat)
   ! 
   ! read referece alpha from file, if any | linh
   ! wherein, the information of n_evc0_fixed, ref_alpha0,
   ! broadening of orbitals will be readed.   
   !
   if (do_orbdep .and. odd_nkscalfact) then
      ! 
      odd_alpha(:) = 0.d0
      !
      call read_odd_nkscalfact(nbspx, odd_alpha)
      !
   endif
   !
   if ( do_nki .or. do_nkipz ) then
      ! 
      write(stdout,2010) do_wxd
      ! 
   endif
   !
   if ( do_orbdep ) then
      !
      write( stdout, "(3x, 'Koopmans memory = ', f10.3, ' MB', /)" ) &
           nksic_memusage( )
      !    
   endif
   !
2001  format(/, 3X,'PZ sic = ',l4 ,/)
2002  format(/, 3X,'NK sic with integral ref = ',l4, / )
2004  format(/, 3X,'NK scaling factor = ',f7.4, / )
2005  format(/, 3X,'Koopmans memory = ',f10.3, ' MB',/)
2006  format(/, 3X,'f_cutoff     = ',f7.4, /)
2010  format(/, 3X,'NK cross-derivatives = ',l4 )
2030  format(/, 3X,'NK applied up to orbital',i7)
   !
   return
   ! 
end subroutine init_nksic
