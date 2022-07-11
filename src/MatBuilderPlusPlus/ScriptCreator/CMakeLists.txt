#project(ScriptCreator)

#cmake_minimum_required(VERSION 3.16.3)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)


#if(${CMAKE_VERSION} VERSION_LESS 3.14)
    #macro(FetchContent_MakeAvailable NAME)
        #FetchContent_GetProperties(${NAME})
        #if(NOT ${NAME}_POPULATED)
            #FetchContent_Populate(${NAME})
            #add_subdirectory(${${NAME}_SOURCE_DIR} ${${NAME}_BINARY_DIR})
        #endif()
    #endmacro()
#endif()




#set(CMAKE_OSX_ARCHITECTURES "x86_64")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

## add extra include directories
##include_directories(/usr/local/include)


## add extra lib directories
##link_directories(/usr/local/opt/libomp/lib /usr/local/lib)
#include(galois)
#include_directories(${galois_SOURCE_DIR}/include)


#add_subdirectory(include/Configurations/ConfigurationMB)
#add_subdirectory(include/Configurations/ConfigurationConfigFile)

#if(APPLE)
    #message("-- Setting OpenMP flags on MacOs. Assuming `brew install libomp`")
    #if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        #set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        #set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        #set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    #endif()
   #if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        #set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        #set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        #set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    #endif()
    #set(CMAKE_OSX_ARCHITECTURES "x86_64")
#else()
    #SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -fopenmp")
#endif()

#set(EXECUTABLE_OUTPUT_PATH ./bin)
#message(STATUS "CPLEX inc path: ${CPLEX_INC} and ${CPLEX_INC2}")

#add_executable(${PROJECT_NAME} ScriptCreator.cpp MatrixTools.cpp Scrambling.cpp MatrixSamplerClass.cpp)

#target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX_INC} ${CPLEX_INC2})
#target_link_directories(${PROJECT_NAME} PRIVATE ${CPLEX_LIB} ${CPLEX_LIB2} )
#target_link_libraries(${PROJECT_NAME} PUBLIC galois++ concert ilocplex cplex m pthread dl)

##install(TARGETS ${PROJECT_NAME} DESTINATION ~/bin)





project(ScriptCreator)

cmake_minimum_required (VERSION 3.13)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(${CMAKE_VERSION} VERSION_LESS 3.14)
    macro(FetchContent_MakeAvailable NAME)
        FetchContent_GetProperties(${NAME})
        if(NOT ${NAME}_POPULATED)
            FetchContent_Populate(${NAME})
            add_subdirectory(${${NAME}_SOURCE_DIR} ${${NAME}_BINARY_DIR})
        endif()
    endmacro()
endif()

# add extra include directories



# add extra lib directories
#link_directories(/usr/local/opt/libomp/lib /usr/local/lib)


#add_subdirectory(Lib/Math)

set(CMAKE_OSX_ARCHITECTURES "x86_64")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")


include(galois)
include_directories(${galois_SOURCE_DIR}/include)


add_subdirectory(include/Configurations/ConfigurationMB)
add_subdirectory(include/Configurations/ConfigurationConfigFile)


if(APPLE)
  set(CPLEX_INC "/Applications/CPLEX_Studio_Community201/cplex/include")
  set(CPLEX_INC2 "/Applications/CPLEX_Studio_Community201/concert/include")
  set(CPLEX_LIB "/Applications/CPLEX_Studio_Community201/cplex/lib/x86-64_osx/static_pic")
  set(CPLEX_LIB2 "/Applications/CPLEX_Studio_Community201/concert/lib/x86-64_osx/static_pic")
else()
  set(CPLEX_INC "/opt/ibm/ILOG/CPLEX_Studio201/cplex/include")
  set(CPLEX_INC2 "/opt/ibm/ILOG/CPLEX_Studio201/concert/include")
  set(CPLEX_LIB "/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic")
  set(CPLEX_LIB2 "/opt/ibm/ILOG/CPLEX_Studio201/concert/lib/x86-64_linux/static_pic")
endif()

#find_package(OpenMP REQUIRED)

#set(INTEGRANDS Lib/Integrands/2D/Heaviside/Heaviside2D_nIntegrands262144_optimSet.cpp 
#Lib/Integrands/2D/Ellipses/SoftEllipses/SoftEllipses2D_nIntegrands524288_optimSet.cpp 
#Lib/Integrands/2D/Ellipses/HardEllipses/Ellipses2D_nIntegrands16384_optimSet.cpp 
#Lib/Integrands/2D/Rectangles/HardRectangles/Rectangles2D_nIntegrands16384_optimSet.cpp 
#Lib/Integrands/2D/Rectangles/SoftRectangles/SoftRectangles2D_nIntegrands16384_optimSet.cpp 
#Lib/Integrands/3D/Heaviside/Heaviside3D.cpp
#Lib/Integrands/3D/SoftEllipses/Gauss3D.cpp
#)

#add_library(Integrands STATIC ${INTEGRANDS} )

#set(INTEGRANDS)
set(EXECUTABLE_OUTPUT_PATH ./bin)
message(STATUS "CPLEX inc path: ${CPLEX_INC} and ${CPLEX_INC2}")

#add_executable(${PROJECT_NAME} MatBuilder.cpp cplexMatrices.cpp MatrixTools.cpp MatrixSamplerClass.cpp Scrambling.cpp)
add_executable(${PROJECT_NAME} ScriptCreator.cpp MatrixTools.cpp Scrambling.cpp MatrixSamplerClass.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CPLEX_INC} ${CPLEX_INC2})
target_link_directories(${PROJECT_NAME} PRIVATE ${CPLEX_LIB} ${CPLEX_LIB2} )
target_link_libraries(${PROJECT_NAME} PUBLIC galois++ concert ilocplex cplex m pthread dl ConfigurationMB ConfigurationConfigFile)




