SRC = $(wildcard *.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = lattice_variables.h
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG) -std=c++0x
LFLAGS = $(DEBUG) -O2 -Wall 

$LATTICE : $(OBJS)
	$(CXX) -o LATTICE $(OBJS) $(LFLAGS)

lattice.o : $(DEPS) lattice.h  lattice.cpp 
	$(CXX) $(CXXFLAGS) lattice.cpp

main.o : $(DEPS) lattice.h main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean:
	\rm *.o *~ LATTICE

