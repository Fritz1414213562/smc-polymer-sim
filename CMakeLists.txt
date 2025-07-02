cmake_minimum_required(VERSION 3.10)
enable_testing()
project(mdmc)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          17)

# for headder include with the path from project root
include_directories(${PROJECT_SOURCE_DIR})


# =============================================================================
# include CMake modules
include(CheckCXXCompilerFlag)

# =============================================================================
# Setting for OpenMM library

# search installed OpenMM lib, include and plugin path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(OpenMM REQUIRED)
message(STATUS "OpenMM include path : ${OPENMM_INCLUDE_DIR}")
message(STATUS "OpenMM library path : ${OPENMM_LIBRARY}")

# set option for openmm plugin directory
get_target_property(EXPECTED_OPENMM_PLUGIN_DIR OpenMM::OpenMM PLUGIN_DIRECTORY)
set(OPENMM_PLUGIN_DIR ${EXPECTED_OPENMM_PLUGIN_DIR} CACHE PATH "The plugin path for OpenMM")
message(STATUS "OpenMM plugin  path : ${OPENMM_PLUGIN_DIR}")

# =============================================================================
# Setting for toml11
set(TOML11_PRECOMPILE ON)

# =============================================================================
# warning flags

set(POLYMER_MDMC_WARNING_FLAGS "")
check_cxx_compiler_flag("-Wall" POLYMER_MDMC_COMPILER_HAS_WALL)
if(POLYMER_MDMC_COMPILER_HAS_WALL)
    set(POLYMER_MDMC_WARNING_FLAGS "${POLYMER_MDMC_WARNING_FLAGS} -Wall")
endif()

message(STATUS "warning flags are ${POLYMER_MDMC_WARNING_FLAGS}")

# =============================================================================
# optimization flags

set(POLYMER_MDMC_OPTIMIZATION_FLAGS "")
option(DEBUG "compile in debug mode" OFF)
if(DEBUG)
    check_cxx_compiler_flag("-O0" POLYMER_MDMC_COMPILER_HAS_O0)
    if(POLYMER_MDMC_COMPILER_HAS_O0)
        set(POLYMER_MDMC_OPTIMIZATION_FLAGS "${POLYMER_MDMC_OPTIMIZATION_FLAGS} -O0")
    endif()
    check_cxx_compiler_flag("-g" POLYMER_MDMC_COMPILER_HAS_G)
    if(POLYMER_MDMC_COMPILER_HAS_O0)
        set(POLYMER_MDMC_OPTIMIZATION_FLAGS "${POLYMER_MDMC_OPTIMIZATION_FLAGS} -g")
    endif()
else()
    check_cxx_compiler_flag("-Ofast" POLYMER_MDMC_COMPILER_HAS_OFAST)
    if(POLYMER_MDMC_COMPILER_HAS_OFAST)
        set(POLYMER_MDMC_OPTIMIZATION_FLAGS "${POLYMER_MDMC_OPTIMIZATION_FLAGS} -Ofast")
    endif()
endif()

message(STATUS "optimization flags are ${POLYMER_MDMC_OPTIMIZATION_FLAGS}")

# =============================================================================
# debugger flags

set(POLYMER_MDMC_DEBUGGER_FLAGS "")
option(GDB "add flag for gdb" OFF)
if(GDB)
    check_cxx_compiler_flag("-ggdb" POLYMER_MDMC_COMPILER_HAS_GGDB)
    if(POLYMER_MDMC_COMPILER_HAS_GGDB)
        set(POLYMER_MDMC_DEBUGGER_FLAGS "${POLYMER_MDMC_DEBUGGER_FLAGS} -ggdb")
    endif()
endif()

message(STATUS "debugger flags are ${POLYMER_MDMC_DEBUGGER_FLAGS}")

# =============================================================================
# compiler-dependent library configuration
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 8.0 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        link_libraries(stdc++fs)
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
    link_libraries(stdc++fs)
endif()

# =============================================================================
# Build the code

add_subdirectory(toml11)
add_subdirectory(fmt)
add_subdirectory(googletest)
add_subdirectory(src)
#add_subdirectory(test)
