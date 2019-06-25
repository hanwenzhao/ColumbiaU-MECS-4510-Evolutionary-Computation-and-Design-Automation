cmake_minimum_required(VERSION 3.12)
# setup the project
#project(HW4 LANGUAGES CXX CUDA)
project(HW4)
# checking openGL
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenMP REQUIRED)
# set c++ standard to 17
set(CMAKE_CXX_STANDARD 17)
# adding cuda as library
add_library(HW4 SHARED
        main.cpp
        #cudaTest.cu
        main.h)
# separate cuda compilation
#set_target_properties( HW4
#        PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# add target
add_executable(HW4_test main.cpp)
#set_property(TARGET HW4_test
#        PROPERTY CUDA_SEPARABLE_COMPILATION ON)
# link library
target_link_libraries(HW4_test
        ${GLUT_LIBRARIES}
        ${GLU_LIBRARIES}
        ${OPENGL_LIBRARIES}
        OpenMP::OpenMP_CXX
        HW4)