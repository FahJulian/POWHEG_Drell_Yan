CXX := g++

COMMONFLAGS := -std=c++20 -Wall -Wextra -MMD -MP
CXXFLAGS := $(COMMONFLAGS) -O2
DEBUGFLAGS := $(COMMONFLAGS) -g -O0 -DPOWHEG_DEBUG

LHAPDF_CXXFLAGS := $(shell lhapdf-config --cflags)
LHAPDF_LDFLAGS  := $(shell lhapdf-config --libs)

CXXFLAGS += $(LHAPDF_CXXFLAGS)
DEBUGFLAGS += $(LHAPDF_CXXFLAGS)
LDFLAGS  += $(LHAPDF_LDFLAGS)

SRC := $(shell find src -name "*.cpp")
OBJ := $(SRC:src/%.cpp=build/release/%.o)
DBG_OBJ := $(SRC:src/%.cpp=build/debug/%.o)

TARGET := build/release/main
DBG_TARGET := build/debug/main

INCLUDES := -Isrc 

release: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

build/release/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

debug: $(DBG_TARGET)

$(DBG_TARGET): $(DBG_OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

build/debug/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(DEBUGFLAGS) $(INCLUDES) -c $< -o $@

DEP := $(OBJ:.o=.d) $(DBG_OBJ:.o=.d)
-include $(DEP)

clean:
	rm -rf build

.PHONY: all debug clean
