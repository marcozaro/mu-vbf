include ./Source/make_opts
SHELL = /bin/bash
LIBDIR = ./lib/
PROG   = driver
LINKLIBS_ME =  -L../lib/ -ldhelas -lmodel
#LINKLIBS_ALL =  -L../lib/ -lmatrix -ldhelas -lmodel
#LIBS = $(LIBDIR)libdhelas.$(libext) $(LIBDIR)libmodel.$(libext)
	

driver: driver.o
	gfortran -o driver driver.o
	
$(LIBDIR)/$(LIBRARY): $(patsubst %.f,%.o,$(wildcard */matrix.f))  all_matrix.o
	$(call CREATELIB, $@, $^)  

