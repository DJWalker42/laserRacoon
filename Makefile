INCDIR=./inc
SRCDIR=./src
OBJDIR=./obj
LIBDIR=./lib

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

LIBNAME=RareMagic

CXX=g++
CXXFLAGS=-Wall -I$(INCDIR)/ -std=c++11
OPTIONS=

INC=$(addprefix $(INCDIR)/, $(HEADERS))
SRC=$(addprefix $(SRCDIR)/, $(SOURCES))
OBJ=$(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))
LIB=$(addprefix $(LIBDIR)/lib, $(addsuffix .a, $(LIBNAME)))

OBJ_DEBUG=$(addprefix $(OBJDIR)/, $(SOURCES:.cpp=_d.o))
LIB_DEBUG=$(addprefix $(LIBDIR)/lib, $(addsuffix _d.a, $(LIBNAME)))

vpath %.cpp $(SRCDIR)

all: lib

lib: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

$(OBJDIR)/%.o: %.cpp $(INC)
	$(CXX) $(CXXFLAGS) $(OPTIONS) -O2 -c $< -o $@ 
	
$(OBJ): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean
clean:
	@rm -f $(OBJ) $(LIB) $(SRCDIR)/*~ $(INCDIR)/*~ ./*~
