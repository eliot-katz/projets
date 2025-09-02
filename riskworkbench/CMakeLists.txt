# CMakeLists.txt (racine)
cmake_minimum_required(VERSION 3.22)
project(rw LANGUAGES CXX)

# Option confort
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Si l'utilisateur n'a pas précisé -DCMAKE_BUILD_TYPE=...
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Bibliothèque centrale ----------------------------------------------------
add_library(rw_core)

target_sources(rw_core PRIVATE
  src/core/normal_rng.cpp
  src/core/stats.cpp
  src/models/gbm.cpp
  src/pricing/analytic_bs.cpp
  src/pricing/mc_pricer.cpp
)

# Publie les headers => #include <rw/...>
target_include_directories(rw_core PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

# C++17
target_compile_features(rw_core PUBLIC cxx_std_17)

# Warnings sympas
if (MSVC)
  target_compile_options(rw_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(rw_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Met tous les binaires dans build/bin (pratique pour lancer)
set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)

# --- CLIs ---------------------------------------------------------------------
add_executable(bs_ref app/cli/bs_ref/main.cpp)
target_link_libraries(bs_ref PRIVATE rw_core)
set_target_properties(bs_ref PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(mc_runner app/cli/mc_runner/main.cpp)
target_link_libraries(mc_runner PRIVATE rw_core)
set_target_properties(mc_runner PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
