cmake_minimum_required(VERSION 3.16.3)

project(Optimize_L2Discrepancy_For_2D_Tiles)

# add extra include directories
include_directories(/usr/local/include)


# add extra lib directories
#link_directories(/usr/local/opt/libomp/lib /usr/local/lib)

add_executable(${PROJECT_NAME} OptimTiles2D.cpp)

add_subdirectory(Lib/Tile)
add_subdirectory(Lib/Tile/Points_2D)

if(APPLE)
    message("-- Setting OpenMP flags on MacOs. Assuming `brew install libomp`")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    endif()
   if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    endif()
    set(CMAKE_OSX_ARCHITECTURES "x86_64")
#else()
#    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp")
endif()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} Points_2D)
target_link_libraries(${PROJECT_NAME} Tiles)

set(EXECUTABLE_OUTPUT_PATH ./bin)

add_executable(OptimTiles2D OptimTiles2D.cpp )
target_link_libraries(OptimTiles2D PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS OptimTiles2D DESTINATION ~/bin)
