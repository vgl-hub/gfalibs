CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalibs
BUILD = build/bin
SOURCE = src
INCLUDE = include
LDFLAGS :=

BINS = $(addsuffix .o, input-filters input-gfa input-agp gfa gfa-lines log functions bed stream-obj uid-generator struct output)

all: $(BINS)
	@

%.o: $(SOURCE)/%.cpp $(INCLUDE)/*.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(basename $@).cpp -o $@
    
clean:
	$(RM) *.o
