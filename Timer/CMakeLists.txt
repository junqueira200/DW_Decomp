cmake_minimum_required(VERSION 3.15)

# 1. Project Declaration
project(
    Timer
    VERSION 1.0.0
    DESCRIPTION "A modern C++ library template"
    LANGUAGES CXX
)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 23)

# 2. Define Library Target
# Options: STATIC, SHARED, or MODULE. Omit to let BUILD_SHARED_LIBS decide.
add_library(Timer STATIC
    src/Alarm.cpp
)

# 3. Target Include Directories
# PUBLIC headers are exposed to targets that link against this library.
# PRIVATE headers are only used to compile this library itself.
target_include_directories(Timer
    PUBLIC
        include/Alarm.h
    PRIVATE
        src/Alarm.cpp
)

# 4. Compiler Warnings and Flags

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fPIC")



