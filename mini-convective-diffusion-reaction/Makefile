# Makefile for mini-convective-diffusion-reaction
# Convection-Diffusion-Reaction Transport Phenomena Module

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm
AR       = ar
ARFLAGS  = rcs

SRC = src/cdr_core.c \
      src/cdr_diffusion.c \
      src/cdr_reaction.c \
      src/cdr_reactor.c \
      src/cdr_mass_transfer.c \
      src/cdr_effectiveness.c

HDR = include/cdr_core.h \
      include/cdr_diffusion.h \
      include/cdr_reaction.h \
      include/cdr_reactor.h \
      include/cdr_mass_transfer.h \
      include/cdr_effectiveness.h

OBJ = $(SRC:.c=.o)

TESTS = tests/test_core \
        tests/test_diffusion \
        tests/test_reaction \
        tests/test_reactor \
        tests/test_effectiveness

EXAMPLES = examples/example_cstr_design \
           examples/example_catalyst_pellet \
           examples/example_dispersion \
           examples/example_rtd_analysis

LIB = libcdr.a

.PHONY: all
all: $(LIB) tests examples

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

src/%.o: src/%.c $(HDR)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: tests
tests: $(TESTS)

tests/test_core: tests/test_core.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

tests/test_diffusion: tests/test_diffusion.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

tests/test_reaction: tests/test_reaction.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

tests/test_reactor: tests/test_reactor.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

tests/test_effectiveness: tests/test_effectiveness.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

.PHONY: examples
examples: $(EXAMPLES)

examples/example_cstr_design: examples/example_cstr_design.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

examples/example_catalyst_pellet: examples/example_catalyst_pellet.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

examples/example_dispersion: examples/example_dispersion.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

examples/example_rtd_analysis: examples/example_rtd_analysis.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

.PHONY: bench
bench: benches/bench_cdr

benches/bench_cdr: benches/bench_cdr.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

.PHONY: demo
demo: demos/demo_conversion

demos/demo_conversion: demos/demo_conversion.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) $< -L. -lcdr -o $@ $(LDFLAGS)

.PHONY: test
test: $(TESTS)
	@echo "=== Running CDR Core Tests ==="
	./tests/test_core
	@echo ""
	@echo "=== Running CDR Diffusion Tests ==="
	./tests/test_diffusion
	@echo ""
	@echo "=== Running CDR Reaction Tests ==="
	./tests/test_reaction
	@echo ""
	@echo "=== Running CDR Reactor Tests ==="
	./tests/test_reactor
	@echo ""
	@echo "=== Running CDR Effectiveness Tests ==="
	./tests/test_effectiveness
	@echo ""
	@echo "=== All tests passed ==="

.PHONY: run-examples
run-examples: $(EXAMPLES)
	@echo "=== Example: CSTR Design ==="
	./examples/example_cstr_design
	@echo ""
	@echo "=== Example: Catalyst Pellet ==="
	./examples/example_catalyst_pellet
	@echo ""
	@echo "=== Example: Dispersion ==="
	./examples/example_dispersion
	@echo ""
	@echo "=== Example: RTD Analysis ==="
	./examples/example_rtd_analysis

.PHONY: count
count:
	@echo "include/ + src/ (.c + .h) line count:"
	@cat include/*.h src/*.c | wc -l

.PHONY: clean
clean:
	rm -f $(OBJ) $(LIB)
	rm -f $(TESTS) $(EXAMPLES)
	rm -f benches/bench_cdr demos/demo_conversion

.PHONY: help
help:
	@echo "Targets: all test run-examples count clean"
