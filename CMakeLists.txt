cmake_minimum_required(VERSION 3.16)
project(span_core)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-O3)
add_executable(span_core main.cpp Graph.cpp Graph.h)
