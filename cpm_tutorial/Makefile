#
# Simple C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#

OBJPATH = obj/
BINPATH = bin/

# Targets
all: 
	@cd $(OBJPATH); \
	make -f Makefile all

.PHONY: clean realclean distclean
clean:
	@cd $(OBJPATH); make clean

distclean: 
	@cd $(OBJPATH); make distclean; \
	cd ../$(BINPATH); rm -f *

realclean:
	@cd $(OBJPATH); make realclean
