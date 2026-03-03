# ─────────────────────────────────────────────────────────────────────────────
# Makefile — genome_mapper
# ─────────────────────────────────────────────────────────────────────────────

CXX      := g++
CXXFLAGS := -O2 -std=c++11 -Wall -Wextra -Wpedantic

SRC_DIR  := src
BIN_DIR  := bin
OBJ_DIR  := build

TARGET   := $(BIN_DIR)/genome_mapper

SRCS     := $(wildcard $(SRC_DIR)/*.cpp)
OBJS     := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR)/*o $(BIN_DIR)/genome_mapper

.PHONY: rebuild
rebuild: clean all

.PHONY: info
info:
	@echo "Sources : $(SRCS)"
	@echo "Objects : $(OBJS)"
	@echo "Target  : $(TARGET)"
