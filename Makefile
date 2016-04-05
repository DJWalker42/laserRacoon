#set up some variables for the various directories
INCDIR=./inc
SRCDIR=./src
OBJDIR=./obj
LIBDIR=./lib

#list the header files
HEADERS= 	StaticVector.h \
			StaticVector.inl\
			StaticMatrix.h \
			StaticMatrix.inl \
			DynVector.h \
			DynVector.inl \
			DynMatrix.h\
			DynMatrix.inl \
			BoundaryConditions.h \
			GaussKnotsWeights.h \
			SecondOrderODE.h \
			BvpODE.h \
			Helpers.h \
			Spline.h \
			Complex.h \
			Interpolation.h \
			State.h \
			DataFit.h \
			LinearSolvers.h \
			Differentials.h \
			Maths.h \
			Node.h \
			ODESolvers.h \
			PascalTri.h \
			Storage.h \
			Storage.inl \
			PhysicalUnits.h \
			FiniteDifferenceGrid.h \
			Quadrature.h \
			Timer.h \
			FormatOutput.h \
			Fourier.h \
			RootSearch.h \
			Visualise.h
			
#list the source files
SOURCES=	BvpODE.cpp \
			FormatOutput.cpp \
			Maths.cpp \
			State.cpp \
			Complex.cpp \
			Fourier.cpp \
			ODESolvers.cpp \
			StaticMatrix.cpp \
			DataFit.cpp \
			Quadrature.cpp \
			Storage.cpp \
			Interpolation.cpp \
			RootSearch.cpp \
			Timer.cpp \
			FiniteDifferenceGrid.cpp \
			LinearSolvers.cpp \
			Spline.cpp \
			Visualise.cpp

#give the library a name
LIBNAME=laserRacoon

#note that any variable declared in this Makefile can be overwritten in the call to make
#For example we could use the command: make CXX=gcc to use a different compiler.
CXX=g++
CXXFLAGS=-Wall -I$(INCDIR)/ -std=c++11
OPTIONS=

#addprefix and add addsuffix are built in functions and do exactly waht they say on the tin
#the syntax in the OBJ assignment means for every source file with the .cpp extenstion replace 
#it with a .o extension
#in unix systems library filenames have to be prefixed with lib and end in a .a extension
INC=$(addprefix $(INCDIR)/, $(HEADERS))
SRC=$(addprefix $(SRCDIR)/, $(SOURCES))
OBJ=$(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))
LIB=$(addprefix $(LIBDIR)/lib, $(addsuffix .a, $(LIBNAME)))

#the following is a suggestion as to how to begin to introduce a debug build version
#it needs finishing
OBJ_DEBUG=$(addprefix $(OBJDIR)/, $(SOURCES:.cpp=_d.o))
LIB_DEBUG=$(addprefix $(LIBDIR)/lib, $(addsuffix _d.a, $(LIBNAME)))

#this line tells make to look in the source directory for the .cpp files to compile
#we need it because this Makefile does not reside in the source directory
#what changes need to be made if you wanted to have this Makefile in the source directory
#and thus can remove this line?
vpath %.cpp $(SRCDIR)

#all is the default target for this Makefile depending on the existance of the library file
all: $(LIB)

#this target actually makes the library using the archive command 'ar'
#it requires that the object files exist and are up-to-date, if not go to that target.
#the options rcs mean (r)eplace object files with the same name, (c)reate the archive.
#the s flag is equivalent to running ranlib on the archive.
#$@ refers to the target. In this context $^ means all object files found in the object directory.
$(LIBDIR)/lib$(LIBNAME).a: $(OBJ)
	ar rcs $@ $^

#here the pipe command '|' means check on the existance of the library directory and
#do this before you try to make the library.
#if it doesn't exist go to the target to make the directory else ignore.	
$(LIB): | $(LIBDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

#this target complies the source files into object files.
#The % symbol is a wildcard meaning do this for each source file found.
#The $< means the first prerequisite i.e. the source file but not the headers
#The headers are a prerequisite as we have templated functions which may change
#the functionality in the library code.
$(OBJDIR)/%.o: %.cpp $(INC)
	$(CXX) $(CXXFLAGS) $(OPTIONS) -O2 -c $< -o $@ 

#note this is the same syntax as for the code that checks for and makes the library directory
#but in this case is for the object directory.	
$(OBJ): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

#clean doesn't actually make anything, it just removes files, so must be "declared" phony.
#we do this as if there is a file called "clean" within this directory then make may get confused.
#the @ symbol suppressed the normal 'echo' response of the command it prefixes.
.PHONY: clean
clean:
	@rm -f $(OBJ) $(LIB) $(SRCDIR)/*~ $(INCDIR)/*~ ./*~
	
