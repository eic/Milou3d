# Modified 15th April 2010 by Thomas Burton

BITS = 32

F_COMP = gfortran

F_FLAGS = -m$(BITS) -g -fno-automatic -fdollar-ok -fno-backslash -finit-local-zero -fno-second-underscore # for use with gfortran43

SRC = \
	benno_code.f bsinit.f calc_kinem2.f calc_kinem2_fixed.f calcphirb.f smear.f decdel.f deceta.f \
	decnst.f decpi0.f decrho.f define_reso.f dmy.f dvcs_end.f dvcs_ini.f \
	f2allm.f fdvcs.f filllujet.f forpaw.f fragpj.f gener_evt.f kc_free.f lorenb8.f \
	readsingle.f \
	mcdvcsv1.f pdgcode_reso.f pxmass.f qed_isr.f rambo.f ranbw.f \
	readsteer.f rotations.f splitp.f h1stuff.f ludata.f print_asc.f

CERN_LIBS = /cern/pro/lib

BASES = libbases.a

LOADER = $(F_COMP)

SUBDIRECTORIES = bases51
#   -@for package in $(DIRECTORY); do \
#	     ( cd $$package ; $(MAKE) ) || exit ; done
#BASESOBJ = ( ls bases51_32/*.o )		  
BASESOBJ = bases51/*.o

dvcs:	gdvcs.o $(SRC:.f=.o) libbases.a
	$(LOADER) \
	gdvcs.o $(SRC:.f=.o) \
	-m$(BITS) \
	$(GCCLIBS) \
	-o milou \
	$(BASES) \
	-L $(CERN_LIBS) -lpythia6205 -ljetset74 -lpacklib -lmathlib -lkernlib

libbases.a:
	( cd bases51 ; make ) || exit ;
	ar rv $@ bases51/*.o
#	@echo $(BASESOBJ)


objects: $(SRC:.f=.o)

%.o : %.f
	$(F_COMP) -c $(F_FLAGS) $< -o $@

.PHONY: clean

clean:
	rm *.o milou libbases.a bases51/*.o

