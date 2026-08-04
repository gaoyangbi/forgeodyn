cmake_minimum_required(VERSION 3.25)

project(forgeodyn
    VERSION 1.0.0
    DESCRIPTION "Fortran implementation of pygeodyn"
    LANGUAGES C Fortran
)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR
        "forgeodyn currently requires the Intel ifx compiler "
        "(CMAKE_Fortran_COMPILER_ID=IntelLLVM)."
    )
endif()

set(MKL_INTERFACE lp64)
set(MKL_LINK dynamic)
set(MKL_THREADING intel_thread)
set(ENABLE_BLAS95 ON)
set(ENABLE_LAPACK95 ON)

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
find_package(MKL CONFIG REQUIRED)
find_package(hdf5 CONFIG REQUIRED)

set(FORGEODYN_SOURCES
    src/arguments_message.f90
    src/constants.f90
    src/deterministic_rng.f90
    src/fortran/Integral_evals.f90
    src/fortran/Legendre_evals.f90
    src/generic/generic_algo.f90
    src/glassofast.f90
    src/inout/priors.f90
    src/inout/reads.f90
    src/linear_interpolation_module.F90
    src/mkl_vsl.f90
    src/utilities.f90
    src/inout/config.f90
    src/common.f90
    src/generic/computer.f90
    src/inout/observations.f90
    src/pca.f90
    src/corestate.f90
    src/augkf/forecaster.f90
    src/augkf/analyser.f90
    src/augkf/augkf_algo.f90
    src/writes.f90
    src/run.f90
    src/main.f90
)

add_executable(forgeodyn ${FORGEODYN_SOURCES})

set_target_properties(forgeodyn PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

target_link_libraries(forgeodyn PRIVATE
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
    MKL::MKL
    hdf5::hdf5_fortran-shared
)

target_compile_options(forgeodyn PRIVATE
    $<$<CONFIG:Debug>:/Od;/traceback;/check:bounds;/check:stack>
    $<$<CONFIG:Release>:/O2>
)

if(WIN32)
    target_link_options(forgeodyn PRIVATE
        "LINKER:/STACK:268435456"
    )
endif()

install(TARGETS forgeodyn RUNTIME DESTINATION bin)
