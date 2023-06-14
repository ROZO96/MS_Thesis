cmake_minimum_required(VERSION 3.17)
project(rd)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

include_directories(.)
include_directories("/usr/local/Cellar/openblas/0.3.21/include")
include_directories("/usr/local/Cellar/libomp/13.0.0/include")
set(OpenBLAS_DIR "/usr/local/Cellar/openblas/0.3.21/lib/cmake/openblas")


find_package(OpenBLAS REQUIRED)
find_package(OpenMP REQUIRED)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_executable(rd
        base.cpp
        constants.h
        inverse.cpp
        io2D.cpp
        main.cpp
        setup2D.cpp
        source2D.cpp
        timestep.cpp
        triangle2D.h
        vertex2D.h
        all_functions.h
        base.h
        inverse.h)
target_link_libraries(rd ${OpenBLAS_LIBRARIES} OpenMP::OpenMP_CXX)

