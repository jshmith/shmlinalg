# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -g -Iinclude/

# Source files
SRC = src/utils.cpp example.cpp

# Object files
OBJ = utils.obj example.obj

# Target executable
TARGET = example.exe

# Build rules
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

utils.obj: src/utils.cpp include/utils.hpp
	$(CXX) $(CXXFLAGS) -c src/utils.cpp -o utils.obj

example.obj: example.cpp include/matrix.hpp include/operations.hpp include/utils.hpp
	$(CXX) $(CXXFLAGS) -c example.cpp -o example.obj

clean:
	del $(OBJ) $(TARGET)