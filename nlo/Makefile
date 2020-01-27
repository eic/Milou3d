#############################################################
#                                                           #
# Makefile         dvcsnlo                                  #
#                                                           # 
#############################################################

OBJECTS= dvcsnlo.o functions.o amplitudes.o crosssection.o tintcross.o funct.o crossh1.o funct1.o bhold.o funct2.o phi1.o phi2.o

$(HOME)/h1/dvcs/nlo/dvcsnlo.exe: Makefile $(OBJECTS) 
	$(LOADER) -o  $(HOME)/h1/dvcs/nlo/dvcsnlo.exe \
	$(OBJECTS)   \
	-L /cern/pro/lib -lkernlib -lmathlib -lpacklib  $(LD_FLAGS) \
	-L$(H1_LIBS) \
	-L/usr/local/lib -lnag_u \


.f.o:
	$(F_COMP) $*.f $(F_FLAGS)








