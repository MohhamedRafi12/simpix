# -------- Compiler --------
CXX := g++
CXXFLAGS := -O3 -std=c++20 -Wall -Wextra

# -------- ROOT (conda) --------
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs --glibs)

# Explicit ASImage libraries (not included in root-config output)
ASIMAGELIBS := -lASImage -lASImageGui

# -------- Target --------
TARGET := simpix
SRC    := simpix_start.cpp

# -------- Default target --------
.DEFAULT_GOAL := run

# -------- Build --------
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) \
		$(ROOTCFLAGS) \
		$(ROOTLIBS) \
		$(ASIMAGELIBS)

# -------- Run --------
run: $(TARGET)
	./$(TARGET) imgA.png imgB.png out

# -------- Clean --------
clean:
	rm -f $(TARGET) *.o

.PHONY: run clean
