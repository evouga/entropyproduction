cmake_minimum_required(VERSION 3.1)
project(entropyproduction)

file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
