# -------- Configuration --------
CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs --glibs)

TARGET := simpix
SRC    := simpix_start.cpp

# -------- Default target --------
.DEFAULT_GOAL := run

# -------- Build rules --------
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@ $(ROOTLIBS)

# -------- Run (default) --------
run: $(TARGET)
	./$(TARGET) imgA.png imgB.png out

# -------- Optional targets --------
debug:
	$(CXX) -g -std=c++17 $(ROOTCFLAGS) $(SRC) -o $(TARGET)_dbg $(ROOTLIBS)

clean:
	rm -f $(TARGET) $(TARGET)_dbg *.o

.PHONY: run debug clean
