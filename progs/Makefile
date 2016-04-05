#note that this makefile is not optimal. Every time you make a change to one of the sample source files
#each source file is recomplied. One solution would be to have a unique target for each source file.

INCDIR=../inc
LIBDIR=../lib
OPENCVLIB=/usr/local/lib
BINDIR=./bin

LIBNAME=laserRacoon

SOURCES=	animatePhys.cpp \
			rootSearchApp.cpp \
			aniPhys.cpp \
			comet.cpp \
			numerov.cpp \
			dataFit.cpp \
			em_main.cpp \
			fresnel.cpp \
			interpolation.cpp \
			laguerre.cpp \
			matrixFactor.cpp \
			dartBoard.cpp \
			nuclearDecay.cpp \
			odeSolvers.cpp \
			quadrat.cpp \
			vanDerPolFFT.cpp \
			BvpODEchk.cpp

CXX=g++
CXXFLAGS=-Wall -O2 -I${INCDIR}/ -std=c++11
LIBFLAGS=-L${LIBDIR}/ -l${LIBNAME} -L${OPENCVLIB}/ -lopencv_core -lopencv_highgui
OPTIONS=

LIB=$(addprefix $(LIBDIR)/lib, $(addsuffix .a, $(LIBNAME)))
BIN=$(addprefix $(BINDIR)/, $(SOURCES:.cpp=))

default: $(BIN)

${BINDIR}/%: %.cpp ${LIB}
	$(CXX) $(CXXFLAGS) $(OPTIONS) $< $(LIBFLAGS) -o $@ 

.PHONY: clean
clean:
	@rm -f $(BIN) ./*~

