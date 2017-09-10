# Makefile for CP/FPMD

include	../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../Modules $(MOD_FLAG)../../CPV/src 
IFLAGS=

FOBJS = \
kcp_main_routine.o 

QEMODS=../../Modules/libqemod.a
QEFFT=../../FFTXlib/libqefft.a
QELA=../../LAXlib/libqela.a
QECP=../../CPV/src/libcp.a 

TLDEPS= bindir mods libs libiotk

all : tldeps libkcp.a kcp.x

libkcp.a : $(FOBJS) $(LOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

kcp.x : libkcp.a kcp_main.o $(QECP) $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o kcp.x kcp_main.o libkcp.a \
		$(QECP) $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../KCP/src/kcp.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/pw_ff.x 

include make.depend
