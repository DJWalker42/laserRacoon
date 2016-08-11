INCDIR=../inc
LIBDIR=../lib
OPENCVLIB=/usr/local/lib
BINDIR=./bin

LIBNAME=laserRacoon

SOURCES=	shm.cpp \
			rootSearchApp.cpp \
			duffing.cpp \
			comet.cpp \
			numerov.cpp \
			dataFit.cpp \
			earth2moon.cpp \
			fresnel.cpp \
			interpolation.cpp \
			laguerre.cpp \
			matrixFactor.cpp \
			dartBoard.cpp \
			nuclearDecay.cpp \
			odeSolvers.cpp \
			quadrat.cpp \
			vanDerPolFFT.cpp \
			BvpODEchk.cpp \
			resonance.cpp \
			double_pendulum.cpp

CXX=g++
CXXFLAGS=-Wall -O2 -I${INCDIR}/ -std=c++11
LIBFLAGS=-L${LIBDIR}/ -l${LIBNAME} -L${OPENCVLIB}/ -lopencv_core -lopencv_highgui
OPTIONS=

LIB=$(addprefix $(LIBDIR)/lib, $(addsuffix .a, $(LIBNAME)))
BIN=$(addprefix $(BINDIR)/, $(SOURCES:.cpp=))

default: $(BIN)

${BINDIR}/%: %.cpp ${LIB}
	$(CXX) $(CXXFLAGS) $(OPTIONS) $< $(LIBFLAGS) -o $@ 
	
$(BIN): | $(BINDIR)

$(BINDIR):
	mkdir -p $(BINDIR) && mkdir -p $(BINDIR)/data

.PHONY: clean
clean:
	@rm -f $(BIN) ./*~

