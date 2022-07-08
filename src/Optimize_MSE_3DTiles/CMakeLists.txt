project(Optimize_MSE_3DTiles)

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
    set(CMAKE_OSX_ARCHITECTURES "x86_64")
else()
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp")
endif()
add_definitions( -DDIM=3 )
find_package(OpenMP REQUIRED)

set(INTEGRANDS Lib/Integrands/Heaviside/Heaviside2D_nIntegrands1048576_optimSet.cpp Lib/Integrands/Ellipses/SoftEllipses/SoftEllipses2D_nIntegrands524288_optimSet.cpp Lib/Integrands/Ellipses/HardEllipses/Ellipses2D_nIntegrands16384_optimSet.cpp Lib/Integrands/Rectangles/HardRectangles/Rectangles2D_nIntegrands16384_optimSet.cpp Lib/Integrands/Rectangles/SoftRectangles/SoftRectangles2D_nIntegrands16384_optimSet.cpp Lib/Integrands/3D/SoftEllipses/Gauss3D.cpp Lib/Integrands/3D/Heaviside/Heaviside3D.cpp)

add_library(Integrands STATIC ${INTEGRANDS} )

set(INTEGRANDS)

set(EXECUTABLE_OUTPUT_PATH ./bin)

add_executable(${PROJECT_NAME} OptimME.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE MatXDynamic Integrands OpenMP::OpenMP_CXX)

install(TARGETS ${PROJECT_NAME} DESTINATION ~/bin)
