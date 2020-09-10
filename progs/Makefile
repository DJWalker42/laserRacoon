# When making individual programs you have to use the command
#	make bin/<prog_name>
#

INCDIR=../inc
LIBDIR=../lib
OCVINC=/opt/local/include
OCVLIB=/opt/local/lib
BINDIR=./bin

LIBNAME=laserRacoon
OCVLIBNAMES=-lopencv_core -lopencv_highgui -lopencv_imgproc

SOURCES=$(wildcard *.cpp)

INCFLAGS=-I${INCDIR}

CXX=g++
CXXFLAGS=-Wall -O2 -std=c++11

LIBFLAGS=-L${LIBDIR} -l${LIBNAME}
OCVINCFLAGS=-I${OCVINC}
LIBFLAGS_OCV=${LIBFLAGS} -L${OCVLIB} ${OCVLIBNAMES}
OPTIONS=

PREFIX=$(CXX) $(CXXFLAGS) ${INCFLAGS} ${OPTIONS}
PREFIX_OCV=$(PREFIX) $(OCVINCFLAGS)

LIB=$(addprefix $(LIBDIR)/lib, $(addsuffix .a, $(LIBNAME)))
BIN=$(addprefix $(BINDIR)/, $(SOURCES:.cpp=))

default: $(BIN)

${BINDIR}/shm: shm.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/rootSearch: rootSearch.cpp ${LIB}
	$(PREFIX) $< $(LIBFLAGS) -o $@

${BINDIR}/duffing: duffing.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/comet: comet.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/numerov: numerov.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -lopencv_imgcodecs -o $@

${BINDIR}/dataFit: dataFit.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/earth2moon: earth2moon.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/fresnel: fresnel.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/laguerre: laguerre.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/matrixFactor: matrixFactor.cpp ${LIB}
	$(PREFIX) $< $(LIBFLAGS) -o $@

${BINDIR}/dartBoard: dartBoard.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/nuclearDecay: nuclearDecay.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/odeSolvers: odeSolvers.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/quadrat: quadrat.cpp ${LIB}
	$(PREFIX) $< $(LIBFLAGS) -o $@

${BINDIR}/vanDerPolFFT: vanDerPolFFT.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/BvpODEchk: BvpODEchk.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/resonance: resonance.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@

${BINDIR}/doublePendulum: doublePendulum.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@
	
${BINDIR}/spring: spring.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@
	
${BINDIR}/car: car.cpp ${LIB}
	$(PREFIX_OCV) $< $(LIBFLAGS_OCV) -o $@
	
${BINDIR}/linearInterp: linearInterp.cpp ${LIB}
	$(PREFIX) $< $(LIBFLAGS) -o $@
	
${BINDIR}/leastLinearSquares: leastLinearSqaures.cpp ${LIB}
	$(PREFIX) $< $(LIBFLAGS) -o $@

$(BIN): | $(BINDIR)

$(BINDIR):
	mkdir -p $(BINDIR) && mkdir -p $(BINDIR)/data

.PHONY: clean
clean:
	@rm -f $(BIN) ./*~
