# Makefile alternativo para compilar os códigos do Capítulo 10

CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O2
LDFLAGS = 
LIBS = 

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

all: Static2D Eigen2D

Static2D: Static2D.cpp $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

Eigen2D: Eigen2D.cpp $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o src/*.o Static2D Eigen2D
