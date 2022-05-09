cmake_minimum_required(VERSION 3.16.3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(OptimTilesND)

add_executable(OptimTiles3D OptimTilesND.cpp)


add_subdirectory(Lib/Tiles/DirectionnalVector)

if(APPLE)
    message("-- Setting OpenMP flags on MacOs. Assuming `brew install libomp`")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
    #IF(EXISTS "/opt/homebrew/")
    #    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp  -I/opt/homebrew/include" CACHE INTERNAL "OpenMP flags for #Xcode toolchain.")
        #set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        #set(OpenMP_omp_LIBRARY "/opt/homebrew/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
      #else()
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
      #endif()
    endif()
   if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
      #IF(EXISTS "/opt/homebrew/")
      #  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp  -I/opt/homebrew/include"  CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
       # set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
       # set(OpenMP_omp_LIBRARY "/opt/homebrew/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
      #else()
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
      #endif()
    endif()
    set(CMAKE_OSX_ARCHITECTURES "x86_64")
endif()

target_compile_options(${PROJECT_NAME} PRIVATE -DDIM=3)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} DirectionnalVector)


set(EXECUTABLE_OUTPUT_PATH ./bin)
