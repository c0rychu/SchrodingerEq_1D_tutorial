has_icc = $(shell which icc >/dev/null; echo $$?)
isAppleClang = $(shell g++ --version |grep "Apple" >/dev/null; echo $$?)
ARCH = $(shell uname)

ifeq ($(ARCH),Darwin) # on macOS
  CXX ?= g++
  CXXFLAGS += -I/opt/local/include # Include path to Eigen (MacPort installation via $sudo port install eigen3)
  ifneq ($(isAppleClang),0)
    CXXFLAGS += -fopenmp # Enable openmp for gcc
  endif
else # on Linux
  ifeq ($(has_icc),0) # using Intel icc with mkl
    CXX = icc
    CXXFLAGS += -qopenmp # Enable openmp for icc
    CXXFLAGS += -I"${MKLROOT}/include"
    CXXFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl
    CXXFLAGS += -DEIGEN_USE_MKL_ALL 
  else
    CXX = g++
    CXXFLAGS += -fopenmp # Enable openmp for gcc
  endif
endif
 

CXXFLAGS += -O3 # Enable optimization


.PHONY: all clean
all: step_potential_rk4

clean:
	rm step_potential_rk4
