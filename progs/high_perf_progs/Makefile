# When making individual programs you have to use the command
#	make <src_file_name>
# where the .cpp extension is omitted
#
# If writing a program that uses the Viewer class (i.e. OpenCV) you have to
# append '_ocv' to the source filename in-order for it to compile correctly
#
# Any of the variables in this Makefile can be overridden on the command line
#
		
all_sources := $(wildcard *.cpp)
all_programs := $(all_sources:%.cpp=%)		
				
ocv_sources := $(wildcard *_ocv.cpp)
ocv_programs_temp := $(ocv_sources:%.cpp=%)			
non_ocv_programs := $(filter-out ${ocv_programs_temp}, ${all_programs})

#remove the _ocv from the binary name
ocv_programs := $(ocv_programs_temp:%_ocv=%)



#-Xpreprocessor is required on Mac OSX to get clang to recognise -fopenmp
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11 -Xpreprocessor -fopenmp

INCDIR=../../inc
LIBDIR=../../lib

# These are locations of OpenCV headers and libraries for Mac OSX
# For Unix/Linux typically these are /usr/include, /usr/lib, or /usr/local/include, /usr/local/lib
OCVINC=/opt/local/include
OCVLIB=/opt/local/lib

OMPINC=${OCVINC}
OMPLIB=/opt/local/lib/libomp


BINDIR=./bin

LIBNAME = -llaserRacoon
OCVLIBNAMES = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs
OMPLIBNAME = -lomp

#Note: opencv_imgcodecs is only used by the numerov program so including it with
# all opencv programs is redundant

INCFLAGS = -I${INCDIR} -I${OMPINC}
LIBFLAGS = -L${LIBDIR} ${LIBNAME} -L${OMPLIB} ${OMPLIBNAME}

OCVINCFLAGS = -I${OCVINC}
OCVLIBFLAGS = -L${OCVLIB} ${OCVLIBNAMES}


# to add any additional options to the compiler command
OPTIONS =

PREFIX = ${CXX} ${CXXFLAGS} ${INCFLAGS} ${OPTIONS}

.PHONY: all
all: ${ocv_programs} ${non_ocv_programs} 

# as we can make individual programs we have to check for the
# existence of the 'bin' directory before building the executable

${ocv_programs}: %: %_ocv.cpp | ${BINDIR}
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
	@mkdir -p $(BINDIR)
	@echo '*** ${BINDIR} created ***'
