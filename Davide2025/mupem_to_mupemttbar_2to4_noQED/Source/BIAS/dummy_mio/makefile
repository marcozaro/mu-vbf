include ../../make_opts

all: dummy_mio

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module dummy_mio
#
dummy_mio.o: dummy_mio.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o dummy_mio.o dummy_mio.f

dummy_mio: dummy_mio.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:	
	@echo "VALID"

