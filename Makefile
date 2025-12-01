CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall

TARGET = round_robin
SOURCE = round_robin.cpp

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

clean:
	rm -f $(TARGET) $(TARGET).exe

.PHONY: all clean

