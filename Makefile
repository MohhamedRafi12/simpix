# ---------------- Configuration ----------------
CXX       := g++
CXXFLAGS  := -O3 -march=native -std=c++17 -Wall -Wextra

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# Image libraries (needed for TASImage)
IMGLIBS := -lASImage 

TARGET := simpix
SRC    := simpix_start.cpp

# ---------------- Rules ----------------
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $^ -o $@ $(ROOTLIBS) $(IMGLIBS)

run: $(TARGET)
	./$(TARGET) imgA.png imgB.png out_AtoB.png

clean:
	rm -f $(TARGET) *.o *.png *.pdf

.PHONY: all run clean
