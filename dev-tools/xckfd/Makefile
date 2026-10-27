# Test tool: finite-difference check of the noncollinear DFPT xc kernel.
# Build after the main tree (make pw ph): make -C dev-tools/xckfd
include ../../make.inc

MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

xckfd.x : xckfd.o dgradcorr_cons.o $(LRMODS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ xckfd.o dgradcorr_cons.o \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)

clean :
	- /bin/rm -f *.x *.o *.mod

include ../../install/tldeps
