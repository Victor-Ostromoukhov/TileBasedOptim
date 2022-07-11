project(OptimMSE2DPipeLine)

cmake_minimum_required (VERSION 3.13)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add extra include directories

# add extra lib directories
#link_directories(/usr/local/opt/libomp/lib /usr/local/lib)


add_subdirectory(Lib/Math)

set(CMAKE_OSX_ARCHITECTURES "x86_64")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")


find_package(OpenMP REQUIRED)

set(INTEGRANDS Lib/Integrands/2D/Heaviside/Heaviside2D_nIntegrands262144_optimSet.cpp 
Lib/Integrands/2D/Ellipses/SoftEllipses/SoftEllipses2D_nIntegrands524288_optimSet.cpp 
Lib/Integrands/2D/Ellipses/HardEllipses/Ellipses2D_nIntegrands16384_optimSet.cpp 
Lib/Integrands/2D/Rectangles/HardRectangles/Rectangles2D_nIntegrands16384_optimSet.cpp 
Lib/Integrands/2D/Rectangles/SoftRectangles/SoftRectangles2D_nIntegrands16384_optimSet.cpp 
Lib/Integrands/3D/Heaviside/Heaviside3D.cpp
Lib/Integrands/3D/SoftEllipses/Gauss3D.cpp
)

add_library(Integrands STATIC ${INTEGRANDS} )

set(INTEGRANDS)
set(EXECUTABLE_OUTPUT_PATH ./bin)



add_executable(${PROJECT_NAME} OptimME.cpp )

target_link_libraries(${PROJECT_NAME} PUBLIC MatXDynamic Integrands OpenMP::OpenMP_CXX)


