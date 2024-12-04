CXX ?= g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS)

TARGET = gfalibs
BUILD = build/bin
SOURCE = src
INCLUDE = include

SOURCES = $(addsuffix .o, input-filters input-gfa input-agp gfa gfa-lines log stream-obj uid-generator struct output memory)

all: $(SOURCES)

%.o: $(SOURCE)/%.cpp $(INCLUDE)/%.h
	$(CXX) $(CXXFLAGS) -c $(SOURCE)/$(basename $@).cpp -o $@
    
clean:
	$(RM) *.o
