cmake_minimum_required(VERSION 3.6)
project(proj)

add_executable(main main.cpp)

target_include_directories(main PRIVATE third-party)