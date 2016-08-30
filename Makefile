#note that any variable declared in this Makefile can be overwritten in the call to make
#For example we could use the command: make CXX=gcc to use a different compiler.

#set up some variables for the various directories
INCDIR=./inc
SRCDIR=./src
OBJDIR=./obj
LIBDIR=./lib

#OpenCV include directory - Visualise module uses OpenCV
#can be omitted if the Visulise module (header and source) removed from project
#note this is for Mac OS X El Capitan. For Unix distros try /usr/local/include
#you can locate OpenCV on your machine using the command: locate opencv
OPENCV_INC=/opt/local/include

INCLUDE=-I$(INCDIR) -I$(OPENCV_INC)

#grab a list of the header files; note we have both .h and .inl file extensions
#we want the library build to depend on the headers.
#actually do we? Maybe on the .inl files that contain templated defintions... to think about.
HEADERS=$(wildcard $(INCDIR)/*.h)
INLINES=$(wildcard $(INCDIR)/*.inl)
HDR=$(HEADERS) $(INLINES)

#grab a list the source files - use wildcard function to avoid the expansion pitfall
SOURCES=$(wildcard $(SRCDIR)/*.cpp)

#give the library a name
LIBNAME=laserRacoon

#complier, flags, and options
CXX=g++
CXXFLAGS=-Wall -std=c++11
OPTIONS=-O2

#patsubst - pattern substitution $(patsubst pattern_to_replace, replace_with_this, text)
# i.e. results in OBJ = (./src/*.cpp) -> ./obj/*.o
OBJ=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

#addprefix and add addsuffix are built in functions and do exactly waht they say on the tin
#in unix systems library filenames have to be prefixed with lib and end in a .a extension
LIB=$(addprefix $(LIBDIR)/lib, $(addsuffix .a, $(LIBNAME)))

#the following is a suggestion as to how to begin to introduce a debug build version
#it needs finishing
OBJ_D=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%_d.o, $(SOURCES))
LIB_D=$(addprefix $(LIBDIR)/lib, $(addsuffix _d.a, $(LIBNAME)))

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
#note the -p option in mkdir command is thus belt and braces
$(LIB): | $(LIBDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

#this target complies the source files into object files.
#The % symbol is a wildcard meaning do this for each source file found.
#The $< means the first prerequisite i.e. the source file but not the headers
#The headers are a prerequisite as we have templated functions which may change
#the functionality in the library code.
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $(OPTIONS) $(INCLUDE) -c $< -o $@

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
