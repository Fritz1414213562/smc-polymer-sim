# =============================================================================
# Setting for use OpenMM library

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
set(CMAKE_CXX_FLAGS "${POLYMER_MDMC_WARNING_FLAGS} ${POLYMER_MDMC_OPTIMIZATION_FLAGS} ${POLYMER_MDMC_DEBUGGER_FLAGS} -DOPENMM_PLUGIN_DIR=${OPENMM_PLUGIN_DIR}")
message(STATUS "CXX compiler flags are ${CMAKE_CXX_FLAGS}")

add_library(mdmclib STATIC
    IntegratorGenerator.cpp
    Observer.cpp
    Simulator.cpp
    SystemGenerator.cpp
    Topology.cpp
    )
target_link_libraries(mdmclib PRIVATE
    OpenMM::OpenMM
    toml11::toml11
    fmt
    )

add_executable(mdmc main.cpp)
set_target_properties(mdmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
target_include_directories(mdmc PRIVATE OpenMM::OpenMM)
target_link_libraries(mdmc PRIVATE
    OpenMM::OpenMM
    toml11::toml11
    fmt
    mdmclib
    )
add_subdirectory(util)
add_subdirectory(input)
add_subdirectory(ff)
