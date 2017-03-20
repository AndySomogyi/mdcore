project(mdcore)

cmake_minimum_required(VERSION 3.7)

option(MDCORE_DOUBLE "build double precision (libmdcore_d) version (single precision is always built)" NO)

option(MDCORE_SHARED "build shared library (static library is alwasy built)" NO)

option(MDCORE_USE_MPI "build the MPI version" NO)

option(MDCORE_USE_FFTW3 "use FFTW" NO)

option(MDCORE_USE_OPENMP "used OpenMP (only available with GCC)" NO)

# Add our module finders to cmake  |Demos/StaticBoundaryDemo/CMakeLists.txt:13:     ${GLUT_LIBRARIES}
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(MDCORE_USE_FFTW3)
  find_package(FFTW REQUIRED)
  message("FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
  message("FFTW_INCLUDES: ${FFTW_INCLUDES}")
endif()

# all products access the public header directory
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)

add_subdirectory(src)
add_subdirectory(examples)
