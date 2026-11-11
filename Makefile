CXX ?= g++
PKG_CONFIG ?= pkg-config

SRC_DIR := Euler21/src
BUILD_DIR := build
TARGET := $(BUILD_DIR)/euler21
DEBUG_BUILD_DIR := $(BUILD_DIR)/debug
DEBUG_TARGET := $(DEBUG_BUILD_DIR)/euler21
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPFILES := $(OBJECTS:.o=.d)
DEBUG_OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(DEBUG_BUILD_DIR)/%.o,$(SOURCES))
DEBUG_DEPFILES := $(DEBUG_OBJECTS:.o=.d)

PKG_CFLAGS := $(shell $(PKG_CONFIG) --cflags sdl2 SDL2_gfx jsoncpp 2>/dev/null)
PKG_LIBS := $(shell $(PKG_CONFIG) --libs sdl2 SDL2_gfx jsoncpp 2>/dev/null)

CPPFLAGS += -I$(SRC_DIR) $(PKG_CFLAGS)
CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++17 -fopenmp
LDFLAGS += -fopenmp
LDLIBS += $(if $(strip $(PKG_LIBS)),$(PKG_LIBS),-lSDL2_gfx -lSDL2 -ljsoncpp)

.PHONY: all debug clean

all: $(TARGET)

debug: $(DEBUG_TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(DEBUG_TARGET): $(DEBUG_OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(@:.o=.d) -c $< -o $@

$(DEBUG_BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -O0 -g3 -std=c++17 -fopenmp -MMD -MP -MF $(@:.o=.d) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

-include $(DEPFILES) $(DEBUG_DEPFILES)
