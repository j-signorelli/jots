# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /opt/apps/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /opt/apps/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work2/08335/jms26/frontera/jots

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work2/08335/jms26/frontera/jots/build

# Include any dependencies generated for this target.
include src/CMakeFiles/jotsLib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/jotsLib.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/jotsLib.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/jotsLib.dir/flags.make

src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.o: src/CMakeFiles/jotsLib.dir/flags.make
src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.o: /work2/08335/jms26/frontera/jots/src/conduction_operator.cpp
src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.o: src/CMakeFiles/jotsLib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work2/08335/jms26/frontera/jots/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.o"
	cd /work2/08335/jms26/frontera/jots/build/src && /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.o -MF CMakeFiles/jotsLib.dir/conduction_operator.cpp.o.d -o CMakeFiles/jotsLib.dir/conduction_operator.cpp.o -c /work2/08335/jms26/frontera/jots/src/conduction_operator.cpp

src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/jotsLib.dir/conduction_operator.cpp.i"
	cd /work2/08335/jms26/frontera/jots/build/src && /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work2/08335/jms26/frontera/jots/src/conduction_operator.cpp > CMakeFiles/jotsLib.dir/conduction_operator.cpp.i

src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/jotsLib.dir/conduction_operator.cpp.s"
	cd /work2/08335/jms26/frontera/jots/build/src && /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work2/08335/jms26/frontera/jots/src/conduction_operator.cpp -o CMakeFiles/jotsLib.dir/conduction_operator.cpp.s

src/CMakeFiles/jotsLib.dir/config_file.cpp.o: src/CMakeFiles/jotsLib.dir/flags.make
src/CMakeFiles/jotsLib.dir/config_file.cpp.o: /work2/08335/jms26/frontera/jots/src/config_file.cpp
src/CMakeFiles/jotsLib.dir/config_file.cpp.o: src/CMakeFiles/jotsLib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work2/08335/jms26/frontera/jots/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/jotsLib.dir/config_file.cpp.o"
	cd /work2/08335/jms26/frontera/jots/build/src && /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/jotsLib.dir/config_file.cpp.o -MF CMakeFiles/jotsLib.dir/config_file.cpp.o.d -o CMakeFiles/jotsLib.dir/config_file.cpp.o -c /work2/08335/jms26/frontera/jots/src/config_file.cpp

src/CMakeFiles/jotsLib.dir/config_file.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/jotsLib.dir/config_file.cpp.i"
	cd /work2/08335/jms26/frontera/jots/build/src && /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work2/08335/jms26/frontera/jots/src/config_file.cpp > CMakeFiles/jotsLib.dir/config_file.cpp.i

src/CMakeFiles/jotsLib.dir/config_file.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/jotsLib.dir/config_file.cpp.s"
	cd /work2/08335/jms26/frontera/jots/build/src && /opt/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work2/08335/jms26/frontera/jots/src/config_file.cpp -o CMakeFiles/jotsLib.dir/config_file.cpp.s

# Object files for target jotsLib
jotsLib_OBJECTS = \
"CMakeFiles/jotsLib.dir/conduction_operator.cpp.o" \
"CMakeFiles/jotsLib.dir/config_file.cpp.o"

# External object files for target jotsLib
jotsLib_EXTERNAL_OBJECTS =

src/libjotsLib.a: src/CMakeFiles/jotsLib.dir/conduction_operator.cpp.o
src/libjotsLib.a: src/CMakeFiles/jotsLib.dir/config_file.cpp.o
src/libjotsLib.a: src/CMakeFiles/jotsLib.dir/build.make
src/libjotsLib.a: src/CMakeFiles/jotsLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work2/08335/jms26/frontera/jots/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libjotsLib.a"
	cd /work2/08335/jms26/frontera/jots/build/src && $(CMAKE_COMMAND) -P CMakeFiles/jotsLib.dir/cmake_clean_target.cmake
	cd /work2/08335/jms26/frontera/jots/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/jotsLib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/jotsLib.dir/build: src/libjotsLib.a
.PHONY : src/CMakeFiles/jotsLib.dir/build

src/CMakeFiles/jotsLib.dir/clean:
	cd /work2/08335/jms26/frontera/jots/build/src && $(CMAKE_COMMAND) -P CMakeFiles/jotsLib.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/jotsLib.dir/clean

src/CMakeFiles/jotsLib.dir/depend:
	cd /work2/08335/jms26/frontera/jots/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work2/08335/jms26/frontera/jots /work2/08335/jms26/frontera/jots/src /work2/08335/jms26/frontera/jots/build /work2/08335/jms26/frontera/jots/build/src /work2/08335/jms26/frontera/jots/build/src/CMakeFiles/jotsLib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/jotsLib.dir/depend

