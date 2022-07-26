CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = gfalibs
BUILD = build/bin
SOURCE = src
INCLUDE = include
LDFLAGS :=

all:
    
$(BUILD):
    -mkdir -p $@
    
clean:
    $(MAKE) -j -C $(GFALIBS_DIR) clean
    $(RM) -r build
