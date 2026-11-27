CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -march=native -Wall -Wextra -Wshadow
BIN      := build
ALGOS    := greedy1_kmb greedy2_sph greedy3_zelikovsky dp_dreyfus_wagner bnb_steiner
TOOLS    :=
HDRS     := $(wildcard include/*.hpp)

.PHONY: all test bench bench-small report zip clean

all: $(addprefix $(BIN)/,$(ALGOS) $(TOOLS))

$(BIN)/%: src/%.cpp $(HDRS) | $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN)/test_fib_heap: tests/test_fib_heap.cpp $(HDRS) | $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN):
	mkdir -p $(BIN)

# Heap tests + amortized-bound check (writes results/amortized_trace.csv)
test: $(BIN)/test_fib_heap
	mkdir -p results
	./$(BIN)/test_fib_heap results/amortized_trace.csv

# Full benchmark (all instances, all five algorithms); see scripts/run_bench.sh
bench: all
	bash scripts/run_bench.sh

# Quick check on the 18 SteinLib B instances only (about a minute)
bench-small: all
	bash scripts/run_bench.sh --small

report:
	cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex

clean:
	rm -rf $(BIN)
