CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalibs
BUILD = build/bin
SOURCE = src
INCLUDE = include
LDFLAGS :=

BINS = $(addsuffix .o, gfa-lines log functions bed stream-obj struct)

all: $(BINS)
	@

%.o: $(SOURCE)/%.cpp $(INCLUDE)/%.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(basename $@).cpp -o $@
    
clean:
	$(RM) *.o
