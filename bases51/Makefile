# Modified 15th April 2010 by Thomas Burton

BITS = 32
#BITS = 64

#F_COMP = g77
#F_COMP = gfortran41
F_COMP = gfortran

#F_FLAGS = -m$(BITS) -g -fno-automatic -fdollar-ok -fno-backslash -finit-local-zero -fno-second-underscore -fugly-logint -ftypeless-boz # for use with g77
#F_FLAGS = -m$(BITS) -g -fno-automatic -fdollar-ok -fno-backslash -fno-second-underscore # For use with gfortran41
F_FLAGS = -m$(BITS) -g -fno-automatic -fdollar-ok -fno-backslash -finit-local-zero -fno-second-underscore # for use with gfortran43

SRC = \
bases.f   bhsave.f  bsetgu.f  bshbok.f	bslist.f  bsputw.f  bsutim.f  dhplot.f	shfill.f  shupdt.f  spinfo.f  xhfill.f	xhrnge.f \
bhinit.f  bschck.f  bsetgv.f  bsinfo.f	bsordr.f  bsread.f  bswrit.f  drn.f	shinit.f  spchck.f  sprgen.f  xhinit.f	xhscle.f \
bhplot.f  bsdate.f  bsgetw.f  bsinit.f	bsparm.f  bstcnv.f  dhfill.f  drnset.f	shplot.f  sphbok.f  spring.f  xhordr.f \
bhrset.f  bsdims.f  bsgrid.f  bsintg.f	bsprnt.f  bstime.f  dhinit.f  shcler.f	shrset.f  sphist.f  xhchck.f  xhplot.f

objects: $(SRC:.f=.o)

%.o : %.f
	$(F_COMP) -c $(F_FLAGS) $< -o $@

.PHONY: clean
clean:
	rm *.o
