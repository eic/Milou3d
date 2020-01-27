
F_COMP = g77

F_FLAGS =-c -g -fno-automatic -fdollar-ok -fno-backslash -finit-local-zero -fno-second-underscore -fugly-logint -ftypeless-boz

CERN = /products/cern/2002/lib

SRC = 	\
	bases.f benno_code.f bsinit.f calc_kinem2.f calc_kinem2_fixed.f calcphirb.f smear.f decdel.f deceta.f \
	decnst.f decpi0.f decrho.f define_reso.f dmy.f dvcs_end.f dvcs_ini.f \
	fdvcs.f filllujet.f forpaw.f fragpj.f gener_evt.f kc_free.f lorenb8.f \
	mcdvcsv1.f pdgcode_reso.f pxmass.f qed_isr.f rambo.f ranbw.f \
	readsteer.f rotations.f splitp.f h1stuff.f ludata.f

dvcs:	gdvcs.o $(SRC:.f=.o)
	$(LOADER) -o dvcs \
	gdvcs.o $(SRC:.f=.o) \
	-L $(H1_LIBS) -lbases \
	-L $(H1_LIBS) -lpythia57 \
	-L $(H1_LIBS) -ljetset74 \
	-L$(CERN_LIBS) -lpacklib -lmathlib -lkernlib \
	$(LD_FLAGS)

.f.o:	
	$(F_COMP) $(F_FLAGS) $*.f

.c.o:	
	gcc -c $*.c



