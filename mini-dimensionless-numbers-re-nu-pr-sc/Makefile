# Makefile — mini-dimensionless-numbers-re-nu-pr-sc
# C library for dimensionless number computation and transport correlations.
#
# Reference: Bird-Stewart-Lightfoot (2007), Incropera & DeWitt (2007),
#            White (2016), Buckingham (1914)
#
# Targets:
#   make          — build library and tests
#   make test     — build and run all tests
#   make examples — build all examples
#   make clean    — remove build artifacts

CC       := gcc
CFLAGS   := -std=c11 -Wall -Wextra -O2 -fPIC
LDFLAGS  := -lm

SRC_DIR  := src
INC_DIR  := include
TST_DIR  := tests
EX_DIR   := examples
BLD_DIR  := build

# Source files
SRCS     := $(SRC_DIR)/dimensionless_numbers.c \
            $(SRC_DIR)/reynolds_number.c \
            $(SRC_DIR)/nusselt_number.c \
            $(SRC_DIR)/prandtl_schmidt.c \
            $(SRC_DIR)/buckingham_pi.c \
            $(SRC_DIR)/transport_correlations.c

OBJS     := $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)

# Test files
TESTS    := $(TST_DIR)/test_dimensionless_numbers.c \
            $(TST_DIR)/test_correlations.c

TEST_EXE := $(TESTS:$(TST_DIR)/%.c=$(BLD_DIR)/test_%.exe)

# Example files
EXAMPLES := $(EX_DIR)/example_pipe_flow_regime.c \
            $(EX_DIR)/example_flat_plate_heat_transfer.c \
            $(EX_DIR)/example_natural_convection.c \
            $(EX_DIR)/example_mass_transfer.c

EX_EXE   := $(EXAMPLES:$(EX_DIR)/%.c=$(BLD_DIR)/example_%.exe)

# Default target
.PHONY: all
all: lib tests examples

# Build directory
$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# Library objects
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/%.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

# Special: objects that depend on multiple headers
$(BLD_DIR)/transport_correlations.o: $(SRC_DIR)/transport_correlations.c \
	$(INC_DIR)/transport_correlations.h \
	$(INC_DIR)/dimensionless_numbers.h \
	$(INC_DIR)/nusselt_number.h \
	$(INC_DIR)/reynolds_number.h \
	$(INC_DIR)/prandtl_schmidt.h | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

# Static library
.PHONY: lib
lib: $(OBJS) | $(BLD_DIR)
	ar rcs $(BLD_DIR)/libdimless.a $(OBJS)

# Tests
.PHONY: tests
tests: lib $(TEST_EXE)

$(BLD_DIR)/test_%.exe: $(TST_DIR)/%.c $(OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< $(OBJS) $(LDFLAGS) -o $@

# Examples
.PHONY: examples
examples: lib $(EX_EXE)

$(BLD_DIR)/example_%.exe: $(EX_DIR)/%.c $(OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< $(OBJS) $(LDFLAGS) -o $@

# Run tests
.PHONY: test
test: tests
	@echo "=== Running test_dimensionless_numbers ==="
	@$(BLD_DIR)/test_test_dimensionless_numbers.exe
	@echo ""
	@echo "=== Running test_correlations ==="
	@$(BLD_DIR)/test_test_correlations.exe
	@echo ""
	@echo "=== All tests passed ==="

# Run all examples
.PHONY: run-examples
run-examples: examples
	@echo "=== Example: Pipe Flow Regime ==="
	@$(BLD_DIR)/example_pipe_flow_regime.exe
	@echo ""
	@echo "=== Example: Flat Plate Heat Transfer ==="
	@$(BLD_DIR)/example_flat_plate_heat_transfer.exe
	@echo ""
	@echo "=== Example: Natural Convection ==="
	@$(BLD_DIR)/example_natural_convection.exe
	@echo ""
	@echo "=== Example: Mass Transfer ==="
	@$(BLD_DIR)/example_mass_transfer.exe

# Code quality
.PHONY: lint
lint:
	@echo "Checking for filler patterns..."
	@grep -rn "_fn[0-9]" $(SRC_DIR) $(INC_DIR) || echo "  No _fnN patterns found (OK)"
	@grep -rn "_aux[0-9]" $(SRC_DIR) $(INC_DIR) || echo "  No _auxN patterns found (OK)"
	@grep -rn "_ext[0-9][^a-z]" $(SRC_DIR) $(INC_DIR) || echo "  No _extN patterns found (OK)"
	@grep -rn "algorithm variant" $(SRC_DIR) $(INC_DIR) || echo "  No filler patterns found (OK)"
	@grep -rn "extension point" $(SRC_DIR) $(INC_DIR) || echo "  No filler patterns found (OK)"
	@echo "Checking for TODO/FIXME/stub..."
	@grep -rn "TODO\|FIXME\|stub\|placeholder" $(SRC_DIR) $(INC_DIR) $(TST_DIR) || echo "  None found (OK)"

# Line count
.PHONY: count
count:
	@echo "include/ files:"
	@wc -l $(INC_DIR)/*.h
	@echo ""
	@echo "src/ files:"
	@wc -l $(SRC_DIR)/*.c
	@echo ""
	@echo "Total include/ + src/:"
	@find $(INC_DIR) $(SRC_DIR) \( -name "*.h" -o -name "*.c" \) -exec cat {} + | wc -l

# Clean
.PHONY: clean
clean:
	rm -rf $(BLD_DIR)

# Full clean including object files that might be outside build/
.PHONY: distclean
distclean: clean
	rm -f $(SRC_DIR)/*.o $(TST_DIR)/*.o $(EX_DIR)/*.o
