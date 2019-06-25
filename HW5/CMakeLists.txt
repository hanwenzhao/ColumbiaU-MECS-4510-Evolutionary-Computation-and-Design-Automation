cmake_minimum_required(VERSION 3.12)
# setup the project
project(HW5 CXX)
# checking openGL
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenMP REQUIRED)
# set c++ standard to 17
set(CMAKE_CXX_STANDARD 17)
# adding cuda as library
add_library(HW5 SHARED
        main.cpp
        main.h)
# add target
add_executable(HW5_test main.cpp)
# link library
target_link_libraries(HW5_test
        ${GLUT_LIBRARIES}
        ${GLU_LIBRARIES}
        ${OPENGL_LIBRARIES}
        OpenMP::OpenMP_CXX
        HW5)