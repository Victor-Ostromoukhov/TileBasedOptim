project(Optimize_MSE_2DTiles)

cmake_minimum_required(VERSION 3.16.3)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# add extra include directories
include_directories(/usr/local/include)


# add extra lib directories
#link_directories(/usr/local/opt/libomp/lib /usr/local/lib)

add_subdirectory(Lib/Math)

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
    #set(CMAKE_OSX_ARCHITECTURES "x86_64")
else()
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp")
endif()

find_package(OpenMP REQUIRED)

#set(INTEGRANDS Lib/Integrands/Heaviside2D_nIntegrands589824_optimSet.cpp Lib/Integrands/SoftEllipses2D_9x256x256_optimSet.cpp )
#add_library(Integrands STATIC ${INTEGRANDS} )
#set(INTEGRANDS)

set(EXECUTABLE_OUTPUT_PATH ./bin)

add_executable(${PROJECT_NAME} OptimME.cpp)

#target_link_libraries(${PROJECT_NAME} PRIVATE MatXDynamic Integrands OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE MatXDynamic OpenMP::OpenMP_CXX)

install(TARGETS ${PROJECT_NAME} DESTINATION ~/bin)
