# When making individual programs you have to use the command
#	make <src_file_name>
# where the .cpp extension is omitted
#
# Any of the variables in this Makefile can be overridden on the command line
#

# add to this list if creating more programs that do NOT use OpenCV
# can we automate this in some way? (This is a genuine question, I don't know the answer)
non_ocv_programs = rootSearch \
				matrixFactor \
				quadrat \
				linearInterp \
				leastLinearSqaures
				
sources = $(wildcard *.cpp)
all_programs = $(sources:%.cpp=%)
ocv_programs = $(filter-out ${non_ocv_programs}, ${all_programs})

CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11

INCDIR=../inc
LIBDIR=../lib
OCVINC=/opt/local/include
OCVLIB=/opt/local/lib
BINDIR=./bin

LIBNAME = -llaserRacoon
OCVLIBNAMES = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs

#Note: opencv_imgcodecs is only used by the numerov program so including it with
# all opencv programs is redundant

INCFLAGS = -I${INCDIR}
LIBFLAGS = -L${LIBDIR} ${LIBNAME}

OCVINCFLAGS = -I${OCVINC}
OCVLIBFLAGS = -L${OCVLIB} ${OCVLIBNAMES}

# to add any additional options to the compiler command
OPTIONS =

PREFIX = ${CXX} ${CXXFLAGS} ${INCFLAGS} ${OPTIONS}

.PHONY: all
all: ${ocv_programs} ${non_ocv_programs} 

# as we can make individual programs we have to check for the
# existence of the 'bin' directory before building the executable

${ocv_programs}: %: %.cpp | ${BINDIR}
	@${PREFIX} ${OCVINCFLAGS} $< ${LIBFLAGS} ${OCVLIBFLAGS} -o $@ 
	@mv $@ ${BINDIR}
	@echo '*** $@ built and installed ***'
	
${non_ocv_programs}: %: %.cpp | ${BINDIR}
	@${PREFIX} $< ${LIBFLAGS} -o $@
	@mv $@ ${BINDIR}
	@echo '*** $@ built and installed ***'
	
.PHONY: clean
clean:
	@rm -rf ${BINDIR}

$(BINDIR):
	@mkdir -p $(BINDIR) && mkdir -p $(BINDIR)/data
	@echo '*** ${BINDIR} created ***'
