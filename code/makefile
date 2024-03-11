include ./Source/make_opts
SHELL = /bin/bash
LIBDIR = ./lib/
PROG   = driver
LINKLIBS_ME =  -Llib/ -ldhelas -lmodel
#LINKLIBS_ALL =  -L../lib/ -lmatrix -ldhelas -lmodel
#LIBS = $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext)



MEFILES = $(patsubst %.f,%.o,$(wildcard */matrix_mod.f))
FILES = ranmar.o driver.o me_wrapper.o $(MEFILES) 

MELIBS: $(MEFILES)


driver: libs $(FILES)
	$(FC) -o driver $(FILES) $(LINKLIBS_ME) 

libs:
	cd Source; make
	
$(LIBDIR)/$(LIBRARY): $(patsubst %.f,%.o,$(wildcard */matrix.f))  all_matrix.o
	$(call CREATELIB, $@, $^)  


clean:
	rm *.o */*.o
	cd Source/MODEL ; make clean
	cd ../DHELAS ; make clean
