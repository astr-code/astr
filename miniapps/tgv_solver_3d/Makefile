# Compiler
FC = nvfortran

# Executable files
EXEC1 = tgv_cpu
EXEC2 = tgv_gpu
EXEC3 = time_stats

# Flags
FLAGS = -O1

SRC1 = tgvsolver_cpu.f90
SRC2 = tgvsolver_gpu.cuf
SRC3 = compstats.f90

# all: $(EXEC1) $(EXEC2) $(EXEC3)
# runs cpu & gpu, then compares
all: cpu gpu

# Compile
# nvfortran -O3 tgvsolver.cuf -o tgv
# $^ is the source file
# $@ is the exec file

# $(EXEC1): $(SRC1)
# 	$(FC) $^ -o $@

cpu: $(SRC1)
	$(FC) $(FLAGS) $^ -o $(EXEC1)
	./$(EXEC1)

gpu: $(SRC2)
	$(FC) $(FLAGS) $^ -o $(EXEC2)
	./$(EXEC2)

compare: $(SRC3)
	$(FC) $^ -o $(EXEC3)
	./$(EXEC3)

clean:
	rm -f $(EXEC1) $(EXEC2) $(EXEC2)

# Phony targets
.PHONY: all clean cpu gpu compare
