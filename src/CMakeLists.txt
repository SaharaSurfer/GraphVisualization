cmake_minimum_required(VERSION 3.8)
project(GraphVisualization LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall")

set(SRC_LIST main.cc graph_visualizator.cc graph_layout.cc graph.cc bmp.cc)
include_directories("header/")

add_executable(${PROJECT_NAME} ${SRC_LIST})