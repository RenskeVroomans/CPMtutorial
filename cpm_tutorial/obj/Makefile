#
# General C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#

#CXX = g++4
CXX = g++


#CPPFLAGS =  -Wall -O3 -march=pentium4 -mmmx -msse -msse2 -fprefetch-loop-arrays
CPPFLAGS =  -Wall -O3 -march=x86-64 -mmmx -msse -msse2 -fprefetch-loop-arrays -DDSFMT_MEXP=521 -DHAVE_SSE2
#CPPFLAGS =  -g -Wall -march=pentium4
#CPPFLAGS = -g -Wall -march=x86-64 -DDSFMT_MEXP=521 -DHAVE_SSE2


#LNFLAGS = -pg -Wl,-rpath,/usr/lib/gcc-4.0.2/lib:/home/renske/local/lib
LNFLAGS = -Wl,-rpath,/usr/lib/gcc-4.0.2/lib:/home/renske/local/lib

MYPATH = /home/renske/local
INCDIR = -I../include -I$(MYPATH)/include -I/usr/include 
LIBDIR = -L$(MYPATH)/lib 

LIBS = -lpng -lgsl -lgslcblas -lboost_program_options 
# Source/object paths
vpath %.cc ../src ../test
vpath %.d .
vpath %.hh ../include
vpath %.o .

BINPATH = ../bin
LIBPATH = ../lib

# Project vars
PROJECT = CPM
LIBRARY = mor
ALL = Agent.o Basic.o Cell.o Neighbours.o dSFMT.o Graphics.o CellAdhesion.o CellLattice.o CellDivision.o IO.o Misc.o SmallWorld.o


OBJECTS = $(ALL)

# Targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS) 
	$(CXX) $(LNFLAGS) $(LIBDIR) $^ -o $(BINPATH)/$@ $(LIBS)

$(LIBRARY): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@.so \
	$(LIBDIR) $(LIBS) $^ -o $(LIBPATH)/$@.so

$(OBJECTS): %.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(INCDIR) $< -o $@

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) $(INCDIR) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(OBJECTS:.o=.d)
endif
endif

.PHONY: clean realclean 
clean:
	@rm -f *.d.* *.o *.d 

realclean: clean 
	@cd ../src; rm -f *~ .*.swp *.bak; \
	cd ../include; rm -f *~ .*.swp *.bak
