# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.21.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.21.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build

# Include any dependencies generated for this target.
include CMakeFiles/new_integrateND_from_file.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/new_integrateND_from_file.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/new_integrateND_from_file.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/new_integrateND_from_file.dir/flags.make

CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o: CMakeFiles/new_integrateND_from_file.dir/flags.make
CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o: ../new_integrateND_from_file.cpp
CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o: CMakeFiles/new_integrateND_from_file.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o"
	/Applications/Xcode.app/Contents/Developer/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o -MF CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o.d -o CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o -c /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/new_integrateND_from_file.cpp

CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.i"
	/Applications/Xcode.app/Contents/Developer/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/new_integrateND_from_file.cpp > CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.i

CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.s"
	/Applications/Xcode.app/Contents/Developer/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/new_integrateND_from_file.cpp -o CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.s

# Object files for target new_integrateND_from_file
new_integrateND_from_file_OBJECTS = \
"CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o"

# External object files for target new_integrateND_from_file
new_integrateND_from_file_EXTERNAL_OBJECTS =

new_integrateND_from_file: CMakeFiles/new_integrateND_from_file.dir/new_integrateND_from_file.cpp.o
new_integrateND_from_file: CMakeFiles/new_integrateND_from_file.dir/build.make
new_integrateND_from_file: /usr/local/opt/libomp/lib/libomp.dylib
new_integrateND_from_file: CMakeFiles/new_integrateND_from_file.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable new_integrateND_from_file"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/new_integrateND_from_file.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/new_integrateND_from_file.dir/build: new_integrateND_from_file
.PHONY : CMakeFiles/new_integrateND_from_file.dir/build

CMakeFiles/new_integrateND_from_file.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/new_integrateND_from_file.dir/cmake_clean.cmake
.PHONY : CMakeFiles/new_integrateND_from_file.dir/clean

CMakeFiles/new_integrateND_from_file.dir/depend:
	cd /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build /Users/ostrom/TileBasedOptim/src/src_new_integrateND_from_file/build/CMakeFiles/new_integrateND_from_file.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/new_integrateND_from_file.dir/depend

