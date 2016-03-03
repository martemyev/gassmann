cmake_minimum_required(VERSION 2.8)

# Declare the name of the project
project(gassmann CXX)

# Extra compilation keys
add_definitions("-std=c++0x -Wall -Wextra -Wshadow -Weffc++")

# Headers and sources
file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/src/*.hpp")
include_directories("${PROJECT_SOURCE_DIR}/src/")

# Declare the application
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
