CC = gcc
CFLAGS = -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
SRCDIR = src
INCDIR = include
TESTSDIR = tests
EXAMPLESDIR = examples
BUILDDIR = build

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c, $(BUILDDIR)/%.o, $(SRCS))
TEST_SRCS = $(wildcard $(TESTSDIR)/*.c)
TEST_BINS = $(patsubst $(TESTSDIR)/%.c, $(BUILDDIR)/test_%, $(TEST_SRCS))
EXAMPLE_SRCS = $(wildcard $(EXAMPLESDIR)/*.c)
EXAMPLE_BINS = $(patsubst $(EXAMPLESDIR)/%.c, $(BUILDDIR)/example_%, $(EXAMPLE_SRCS))

.PHONY: all clean test examples

all: $(OBJS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/test_%: $(TESTSDIR)/%.c $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

$(BUILDDIR)/example_%: $(EXAMPLESDIR)/%.c $(OBJS)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

test: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo "Running $$t..."; ./$$t && echo "PASS" || echo "FAIL"; done

examples: $(EXAMPLE_BINS)
	@for e in $(EXAMPLE_BINS); do echo "Running $$e..."; ./$$e || true; done

clean:
	rm -rf $(BUILDDIR)
