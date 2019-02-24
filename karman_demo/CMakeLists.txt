cmake_minimum_required(VERSION 3.12)
project(karman_filter)

set(CMAKE_CXX_STANDARD 14)

# 添加Eigen头文件
include_directories( "/usr/include/eigen3" )

add_executable(karman_filter kf_demo.cpp)
