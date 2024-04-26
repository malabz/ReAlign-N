# Define compiler and compile options
CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall

# Define targets and dependencies
TARGET = realign_n
SRCS = src/main.cpp src/Fasta.cpp
OBJS = $(SRCS:.cpp=.o)
INSTALL_DIR = $(HOME)/.realign_n/bin

# Detect user's shell
SHELL_NAME := $(shell basename $$SHELL)

all: $(TARGET) global_realignment 

# Link object files
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile the source code
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

global_realignment:
	mkdir -p $(INSTALL_DIR)
	$(MAKE) -C src/global_realignment
	cp src/global_realignment/bin/Release/global_realignment $(INSTALL_DIR)
ifeq ($(SHELL_NAME),bash)
	echo 'export PATH="$(INSTALL_DIR):$$PATH"' >> ~/.bashrc
else ifeq ($(SHELL_NAME),zsh)
	echo 'export PATH="$(INSTALL_DIR):$$PATH"' >> ~/.zshrc
else ifeq ($(SHELL_NAME),fish)
	echo 'set -gx PATH $(INSTALL_DIR) $$PATH' >> ~/.config/fish/config.fish
else ifeq ($(SHELL_NAME),csh)
	echo 'setenv PATH $(INSTALL_DIR):$$PATH' >> ~/.cshrc
else ifeq ($(SHELL_NAME),tcsh)
	echo 'setenv PATH $(INSTALL_DIR):$$PATH' >> ~/.tcshrc
else
	@echo "Unsupported shell: $(SHELL_NAME). Please add the following line manually to your shell configuration file:"
	@echo "export PATH=\"$(INSTALL_DIR):$$PATH\""
endif

# Clean
clean:
	rm -f $(TARGET) $(OBJS)
	$(MAKE) -C src/global_realignment clean

.PHONY: all $(TARGET) global_realignment clean
