# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion-2018.2.6/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2018.2.6/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/hanwen/CLionProjects/HW5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/hanwen/CLionProjects/HW5/cmake-build-default

# Include any dependencies generated for this target.
include CMakeFiles/HW5_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/HW5_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/HW5_test.dir/flags.make

CMakeFiles/HW5_test.dir/main.cpp.o: CMakeFiles/HW5_test.dir/flags.make
CMakeFiles/HW5_test.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/hanwen/CLionProjects/HW5/cmake-build-default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/HW5_test.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HW5_test.dir/main.cpp.o -c /home/hanwen/CLionProjects/HW5/main.cpp

CMakeFiles/HW5_test.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HW5_test.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/hanwen/CLionProjects/HW5/main.cpp > CMakeFiles/HW5_test.dir/main.cpp.i

CMakeFiles/HW5_test.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HW5_test.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/hanwen/CLionProjects/HW5/main.cpp -o CMakeFiles/HW5_test.dir/main.cpp.s

# Object files for target HW5_test
HW5_test_OBJECTS = \
"CMakeFiles/HW5_test.dir/main.cpp.o"

# External object files for target HW5_test
HW5_test_EXTERNAL_OBJECTS =

HW5_test: CMakeFiles/HW5_test.dir/main.cpp.o
HW5_test: CMakeFiles/HW5_test.dir/build.make
HW5_test: /usr/lib/x86_64-linux-gnu/libglut.so
HW5_test: /usr/lib/x86_64-linux-gnu/libXmu.so
HW5_test: /usr/lib/x86_64-linux-gnu/libXi.so
HW5_test: /usr/lib/x86_64-linux-gnu/libGL.so
HW5_test: /usr/lib/x86_64-linux-gnu/libGLU.so
HW5_test: libHW5.so
HW5_test: /usr/lib/gcc/x86_64-linux-gnu/5/libgomp.so
HW5_test: /usr/lib/x86_64-linux-gnu/libpthread.so
HW5_test: CMakeFiles/HW5_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/hanwen/CLionProjects/HW5/cmake-build-default/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable HW5_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/HW5_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/HW5_test.dir/build: HW5_test

.PHONY : CMakeFiles/HW5_test.dir/build

CMakeFiles/HW5_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/HW5_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/HW5_test.dir/clean

CMakeFiles/HW5_test.dir/depend:
	cd /home/hanwen/CLionProjects/HW5/cmake-build-default && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/hanwen/CLionProjects/HW5 /home/hanwen/CLionProjects/HW5 /home/hanwen/CLionProjects/HW5/cmake-build-default /home/hanwen/CLionProjects/HW5/cmake-build-default /home/hanwen/CLionProjects/HW5/cmake-build-default/CMakeFiles/HW5_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/HW5_test.dir/depend

