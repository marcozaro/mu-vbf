include ./Source/make_opts
SHELL = /bin/bash
LIBDIR = ./lib/
PROG   = driver
LINKLIBS_ME =  -Llib/ -ldhelas -lmodel 
#LINKLIBS_ALL =  -L../lib/ -lmatrix -ldhelas -lmodel
#LIBS = $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext)

emela = /home/mg_user/private_MG5/eMELA_muon_install/bin/eMELA-config
lhapdf = lhapdf-config
#LIBEMELA = $(shell $(emela) --ldflags) $(shell $(lhapdf) --ldflags) -lstdc++
 



MEFILES = $(patsubst %.f,%.o,$(wildcard */matrix_mod.f))
FILES = ranmar.o boostwdir2.o pdf.o vegas.o cuts.o me_wrapper.o qprime.o HwU.o analysis_HwU_template.o kinematics.o integrands.o $(MEFILES) 
DRIVER = driver.o $(FILES)
CHECK1 = check_1coll.o $(FILES)



driver: libs $(DRIVER)
	$(FC) -o driver $(DRIVER) $(LINKLIBS_ME) $(LIBEMELA) 

check_1coll: libs $(CHECK1)
	$(FC) -o check_1coll $(CHECK1) $(LINKLIBS_ME)  $(LIBEMELA) 

cuts.o: input.inc
driver.o: input.inc
pdf.o: input.inc
qprime.o: input.inc
analysis_HwU_template.o : input.inc
me_wrapper.o : input.inc
integrands.o : input.inc

libs:
	cd Source; make
	
$(LIBDIR)/$(LIBRARY): $(patsubst %.f,%.o,$(wildcard */matrix.f))  all_matrix.o
	$(call CREATELIB, $@, $^)  


clean:
	rm *.o */*.o
	cd Source/MODEL ; make clean
	cd Source/DHELAS ; make clean
