cmake_minimum_required(VERSION 3.16.3)

project(OptimTiles2D)

add_executable(${PROJECT_NAME} OptimTiles2D.cpp)

add_subdirectory(Lib/Tile)
add_subdirectory(Lib/Tile/Points_2D)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} Points_2D)
target_link_libraries(${PROJECT_NAME} Tiles)
