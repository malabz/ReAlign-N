# Define compiler and compile options
CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall

all: realign_n

realign_n:
    $(MAKE) -C src/*.cpp

global_realignment:
    $(MAKE) -C src/global_realignment

clean:
    $(MAKE) -C src/*.cpp clean
    $(MAKE) -C src/global_realignment clean

.PHONY: all realign_n global_realignment clean