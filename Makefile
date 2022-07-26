CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalibs
BUILD = build/bin
SOURCE = src
INCLUDE = include
LDFLAGS :=

BINS = gfa-lines log

$(BINS): %: $(SOURCE)/%.cpp $(INCLUDE)/%.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$@.cpp -o $@.o
    
clean:

	$(RM) *.o
