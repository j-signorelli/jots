# Require same CMake as MFEM
cmake_minimum_required(VERSION 3.10)

# Define our project
project(jots VERSION 2.0 LANGUAGES CXX)

# Create interface library for compiler flag requirements (requiring C++11)
add_library(compiler_flag_requirements INTERFACE)
target_compile_features(compiler_flag_requirements INTERFACE cxx_std_11)

# Get preCICE
find_package(PRECICE REQUIRED CONFIG NAMES precice)
message(STATUS "Found preCICE config: TRUE")

# Get MFEM
find_package(MFEM REQUIRED CONFIG)
message(STATUS "Found MFEM config: TRUE")

# Get Boost headers
find_package(Boost REQUIRED CONFIG COMPONENTS headers)
message(STATUS "Found Boost headers: TRUE")

# Add src directory (creates jotsLib)
add_subdirectory(src)

# Add regression tests if specified
option(ENABLE_TESTING "Enable testing" OFF)
if (ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Create executable
add_executable(jots jots_main.cpp)
target_link_libraries(jots PUBLIC jotsLib)

message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")

# Install
install(TARGETS jots DESTINATION bin)
install(FILES "${PROJECT_SOURCE_DIR}/pw_mfem_mesh_fix.py" DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)