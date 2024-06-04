cmake_minimum_required(VERSION 3.28)
project(cpp_structures)

set(CMAKE_CXX_STANDARD 20)

add_executable(cpp_structures main.cpp
        lib/string.h
        lib/stackallocator.h
        lib/list.h
        lib/unordered_map.h
        lib/deque.h
        lib/geometry.h
        lib/biginteger.h
        lib/matrix.h)
