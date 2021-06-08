#CXX = clang++
CXX = g++

CXXFLAGS += -I/opt/local/include          # Include path to Eigen
#CXXFLAGS += -fopenmp                      # Enable openmp
CXXFLAGS += -O3                           # Enable optimization

.PHONY: all clean
all: step_potential_rk4

clean:
	rm step_potential_rk4
