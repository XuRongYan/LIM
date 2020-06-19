cmake_minimum_required(VERSION 3.16)
project(LIM)

set(GOOGLETEST_VERSION 1.8.0)
set(CMAKE_CXX_STANDARD 11)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")

set(LIB_SURFACE_MESH /usr/local/lib/libSurfaceMesh.1.0.dylib)
set(BOOST_FILESYSTEM /usr/local/lib/libboost_filesystem.dylib)
set(BOOST_PROGRAM_OPTIONS /usr/local/lib/libboost_program_options.dylib)
set(OPENMP /usr/local/lib/libomp.dylib)

add_subdirectory(external/googletest)

include_directories(/usr/local/include)
include_directories(/external/googletest/include)

find_package(Eigen3 3.3.0)
if (Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else ()
    message(can't find eigen3)
    return()
endif ()

find_package(GTest REQUIRED)

aux_source_directory(. DIR_SRC)
aux_source_directory(./tests TEST_SRC)
aux_source_directory(./utils UTILS_SRC)
aux_source_directory(./energy ENERGY_SRC)
aux_source_directory(./solver SOLVER_SRC)

add_executable(${PROJECT_NAME} ${DIR_SRC} ${UTILS_SRC} ${ENERGY_SRC} ${SOLVER_SRC})
add_executable(test_main ${TEST_SRC} external/googletest ${UTILS_SRC} ${ENERGY_SRC} ${SOLVER_SRC})

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES)
endif (MSVC)

target_link_libraries(${PROJECT_NAME} ${LIB_SURFACE_MESH} ${BOOST_FILESYSTEM} ${BOOST_PROGRAM_OPTIONS} ${OPENMP})
target_link_libraries(test_main ${GTEST_MAIN_LIBRARY} ${GTEST_LIBRARY} ${LIB_SURFACE_MESH} ${OPENMP})