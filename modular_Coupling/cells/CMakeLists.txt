cmake_minimum_required(VERSION 3.19.3)
project(cells)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags for better optimization and debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})
find_package(OpenMP REQUIRED)

set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}")

# Include modular subdirectories
include_directories("${PROJECT_SOURCE_DIR}/config")
include_directories("${PROJECT_SOURCE_DIR}/agents")
include_directories("${PROJECT_SOURCE_DIR}/behaviors")
include_directories("${PROJECT_SOURCE_DIR}/io")
include_directories("${PROJECT_SOURCE_DIR}/environment")

# External bdm_adapter header-only path (built separately)
get_filename_component(REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(BDM_ADAPTER_DIR "${REPO_ROOT}/adapters/bdm-adapter" CACHE PATH "Path to bdm-adapter source")
include_directories("${BDM_ADAPTER_DIR}/src")

# Find preCICE (try pkg-config first, then fallback to find_library)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PRECICE QUIET libprecice)
endif()

if(NOT PRECICE_FOUND)
  find_library(PRECICE_LIBRARY precice)
  if(PRECICE_LIBRARY)
    set(PRECICE_LIBRARIES ${PRECICE_LIBRARY})
    set(PRECICE_FOUND TRUE)
  else()
    message(FATAL_ERROR "Could not find the preCICE library. Install preCICE or set PRECICE_LIBRARY to its location.")
  endif()
endif()

# Collect source files from all modular subdirectories
file(GLOB_RECURSE PROJECT_HEADERS "${PROJECT_SOURCE_DIR}/*.h")
file(GLOB_RECURSE PROJECT_SOURCES "${PROJECT_SOURCE_DIR}/*.cc")

# bdm-adapter is header-only, no library to link
# Just include the headers (already done above)

bdm_add_executable(${CMAKE_PROJECT_NAME}
                   HEADERS ${PROJECT_HEADERS}
                   SOURCES ${PROJECT_SOURCES}
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES} ${PRECICE_LIBRARIES} OpenMP::OpenMP_CXX)

enable_testing()
add_executable(boundary_geometry_test tests/boundary_geometry_test.cc)
target_include_directories(boundary_geometry_test PRIVATE "${BDM_ADAPTER_DIR}/src")
add_test(NAME boundary_geometry COMMAND boundary_geometry_test)

add_executable(rons_intracellular_model_test
               tests/rons_intracellular_model_test.cc)
target_include_directories(rons_intracellular_model_test
                           PRIVATE "${PROJECT_SOURCE_DIR}/behaviors")
add_test(NAME rons_intracellular_model COMMAND rons_intracellular_model_test)

# Add custom targets for different configurations
add_custom_target(cavity
    COMMAND ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}
    DEPENDS ${CMAKE_PROJECT_NAME}
    COMMENT "Running cavity coupling simulation"
)

add_custom_target(full
    COMMAND ${CMAKE_COMMAND} -E env BDM_SIMULATION_MODE=full ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}
    DEPENDS ${CMAKE_PROJECT_NAME}
    COMMENT "Running full coupling simulation with temperature"
)

add_custom_target(validate
    COMMAND ${CMAKE_COMMAND} -E env BDM_SIMULATION_MODE=validate ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}
    DEPENDS ${CMAKE_PROJECT_NAME}
    COMMENT "Validating coupling configuration"
)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Coupling Configuration ===")
message(STATUS "BioDynaMo found: ${BioDynaMo_FOUND}")
message(STATUS "preCICE found: ${PRECICE_FOUND}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
