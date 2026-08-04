# forgeodyn Makefile
# Target platform: Linux with Intel oneAPI (mpiifx, Intel MPI and oneMKL).

SHELL := /bin/sh

FC := mpiifx
TARGET := forgeodyn
BUILD_DIR := build/make
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod

HDF5_ROOT ?=
HDF5_INC ?= $(HDF5_ROOT)/include
HDF5_LIB ?= $(HDF5_ROOT)/lib

BUILD ?= release

BASE_FFLAGS := -fpp -qopenmp -module $(MOD_DIR) -I$(MOD_DIR) \
	-I$(HDF5_INC) -I$(MKLROOT)/include \
	-I$(MKLROOT)/include/mkl/intel64/lp64

ifeq ($(BUILD),debug)
  CONFIG_FFLAGS := -O0 -g -traceback -check bounds -check stack
else ifeq ($(BUILD),release)
  CONFIG_FFLAGS := -O2
else
  $(error BUILD must be either release or debug)
endif

FFLAGS := $(BASE_FFLAGS) $(CONFIG_FFLAGS) $(EXTRA_FFLAGS)
LDLIBS := -L$(HDF5_LIB) -lhdf5_fortran -lhdf5 \
	-L$(MKLROOT)/lib/intel64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 \
	-qmkl=parallel -qopenmp $(EXTRA_LIBS)

# Keep this list in module-dependency order. The build is deliberately serial
# because Fortran .mod files must exist before dependent sources are compiled.
SOURCES := \
	src/constants.f90 \
	src/mkl_vsl.f90 \
	src/deterministic_rng.f90 \
	src/fortran/Integral_evals.f90 \
	src/fortran/Legendre_evals.f90 \
	src/linear_interpolation_module.F90 \
	src/utilities.f90 \
	src/inout/priors.f90 \
	src/inout/reads.f90 \
	src/inout/config.f90 \
	src/common.f90 \
	src/generic/computer.f90 \
	src/inout/observations.f90 \
	src/pca.f90 \
	src/corestate.f90 \
	src/generic/generic_algo.f90 \
	src/glassofast.f90 \
	src/augkf/forecaster.f90 \
	src/augkf/analyser.f90 \
	src/augkf/augkf_algo.f90 \
	src/writes.f90 \
	src/run.f90 \
	src/arguments_message.f90 \
	src/main.f90

OBJECTS := $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(SOURCES)))

.PHONY: all clean check-env help
.NOTPARALLEL:

all: check-env $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(OBJECTS) $(LDLIBS) -o $@

$(OBJ_DIR)/%.f90.o: %.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.F90.o: %.F90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

check-env:
	@test -n "$(MKLROOT)" || { echo "Error: MKLROOT is not set. Source oneAPI setvars.sh first."; exit 1; }
	@test -n "$(I_MPI_ROOT)" || { echo "Error: I_MPI_ROOT is not set. Source oneAPI setvars.sh first."; exit 1; }
	@test -n "$(HDF5_ROOT)" || { echo "Error: HDF5_ROOT is not set."; exit 1; }
	@test -d "$(HDF5_INC)" || { echo "Error: HDF5 include directory not found: $(HDF5_INC)"; exit 1; }
	@test -d "$(HDF5_LIB)" || { echo "Error: HDF5 library directory not found: $(HDF5_LIB)"; exit 1; }

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

help:
	@echo "Usage:"
	@echo "  source /path/to/oneapi/setvars.sh"
	@echo "  make HDF5_ROOT=/path/to/hdf5"
	@echo "  make BUILD=debug HDF5_ROOT=/path/to/hdf5"
	@echo "  make clean"
