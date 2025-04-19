# Makefile alternativo para compilar os códigos do Capítulo 10

CXX = g++
CXXFLAGS = -std=c++17 -I include -O2 -I /usr/include/eigen3
LDFLAGS = 
LIBS = 

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
sudo apt install cmake
all: Static2D Eigen2D

Static2D: Static2D.cpp $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

Eigen2D: Eigen2D.cpp $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o src/*.o Static2D Eigen2D
