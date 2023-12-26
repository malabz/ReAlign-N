# Define compiler and compile options
CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall

# Define targets and dependencies
TARGET = realign_n
SRCS = src/main.cpp src/Fasta.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET) global_realignment 

# Link object files
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile the source code
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

global_realignment:
	$(MAKE) -C src/global_realignment

# Clean
clean:
	rm -f $(TARGET) $(OBJS)
	$(MAKE) -C src/global_realignment clean

.PHONY: all $(TARGET) global_realignment clean