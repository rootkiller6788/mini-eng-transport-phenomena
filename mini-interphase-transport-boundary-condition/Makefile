CC       = gcc
CFLAGS   = -Wall -std=c11 -O2 -g
LDFLAGS  = -lm

INCLUDE_DIR = include
SRC_DIR     = src
TEST_DIR    = tests
EXAMPLE_DIR = examples
BUILD_DIR   = build

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

TESTS = $(wildcard $(TEST_DIR)/*.c)
TEST_BINS = $(TESTS:$(TEST_DIR)/%.c=$(BUILD_DIR)/test_%)

EXAMPLES = $(wildcard $(EXAMPLE_DIR)/*.c)
EXAMPLE_BINS = $(EXAMPLES:$(EXAMPLE_DIR)/%.c=$(BUILD_DIR)/example_%)

.PHONY: all test examples clean

all: test examples

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

$(BUILD_DIR)/test_%: $(TEST_DIR)/%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< $(OBJS) $(LDFLAGS) -o $@

$(BUILD_DIR)/example_%: $(EXAMPLE_DIR)/%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< $(OBJS) $(LDFLAGS) -o $@

test: $(TEST_BINS)
	@echo "Running tests..."
	@failed=0; for t in $(TEST_BINS); do echo "  $$t"; $$t || failed=1; done; if [ $$failed -eq 0 ]; then echo "All tests passed."; else echo "Some tests FAILED."; exit 1; fi

examples: $(EXAMPLE_BINS)
	@echo "Examples built successfully."

clean:
	rm -rf $(BUILD_DIR)
	@echo "Cleaned."

lines:
	@echo "include/ lines: $$(cat include/*.h | wc -l)"
	@echo "src/ lines:     $$(cat src/*.c | wc -l)"
