# Makefile for FF_kc_bse

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules $(MOD_FLAG).
IFLAGS=

PWFFOBJS = \
io_bsepot_xml.o \
finite_field_pert.o \
input_ff.o \
run_pwscf_ff.o 

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs pw-lib lrmods

all : tldeps libpwff.a pw_ff.x 

libpwff.a : $(PWFFOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

pw_ff.x : libpwff.a pwscf_ff.o $(PWOBJS) $(LRMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ pwscf_ff.o libpwff.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../FF_kc_bse/src/$@ .)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/pw_ff.x 

#include make.depend
# DO NOT DELETE
