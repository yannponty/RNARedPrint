# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.4

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.4.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.4.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master

# Include any dependencies generated for this target.
include test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/depend.make

# Include the progress variables for this target.
include test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/progress.make

# Include the compile flags for this target's objects.
include test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/flags.make

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/flags.make
test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o: test/htd/HypertreeDecompositionAlgorithmTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o"
	cd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o -c /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd/HypertreeDecompositionAlgorithmTest.cpp

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.i"
	cd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd/HypertreeDecompositionAlgorithmTest.cpp > CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.i

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.s"
	cd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd/HypertreeDecompositionAlgorithmTest.cpp -o CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.s

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.requires:

.PHONY : test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.requires

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.provides: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.requires
	$(MAKE) -f test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/build.make test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.provides.build
.PHONY : test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.provides

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.provides.build: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o


# Object files for target HypertreeDecompositionAlgorithmTest
HypertreeDecompositionAlgorithmTest_OBJECTS = \
"CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o"

# External object files for target HypertreeDecompositionAlgorithmTest
HypertreeDecompositionAlgorithmTest_EXTERNAL_OBJECTS =

test/htd/HypertreeDecompositionAlgorithmTest: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o
test/htd/HypertreeDecompositionAlgorithmTest: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/build.make
test/htd/HypertreeDecompositionAlgorithmTest: lib/libhtd.0.0.0.dylib
test/htd/HypertreeDecompositionAlgorithmTest: test/googletest/googletest-release-1.7.0/libgtest_main.dylib
test/htd/HypertreeDecompositionAlgorithmTest: test/googletest/googletest-release-1.7.0/libgtest.dylib
test/htd/HypertreeDecompositionAlgorithmTest: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable HypertreeDecompositionAlgorithmTest"
	cd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/build: test/htd/HypertreeDecompositionAlgorithmTest

.PHONY : test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/build

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/requires: test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/HypertreeDecompositionAlgorithmTest.cpp.o.requires

.PHONY : test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/requires

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/clean:
	cd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd && $(CMAKE_COMMAND) -P CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/cmake_clean.cmake
.PHONY : test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/clean

test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/depend:
	cd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd /Users/wei/Desktop/PhD_experiment/RNARedPrint-master/lib/htd-master/test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/htd/CMakeFiles/HypertreeDecompositionAlgorithmTest.dir/depend

