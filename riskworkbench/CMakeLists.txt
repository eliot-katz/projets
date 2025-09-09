# CMakeLists.txt (racine)
cmake_minimum_required(VERSION 3.22)
project(rw LANGUAGES CXX)

# Option confort
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type par défaut
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Où déposer les binaires
set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)

# =======================
# Bibliothèque centrale
# =======================
add_library(rw_core)

target_sources(rw_core PRIVATE
  src/core/normal_rng.cpp
  src/core/stats.cpp
  src/models/gbm.cpp
  src/pricing/analytic_bs.cpp
  src/pricing/mc_pricer.cpp
  src/pricing/greeks.cpp
)

target_include_directories(rw_core PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(rw_core PUBLIC cxx_std_17)

if (MSVC)
  target_compile_options(rw_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(rw_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# =======================
# CLIs
# =======================
add_executable(bs_ref app/cli/bs_ref/main.cpp)
target_link_libraries(bs_ref PRIVATE rw_core)
set_target_properties(bs_ref PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(mc_runner app/cli/mc_runner/main.cpp)
target_link_libraries(mc_runner PRIVATE rw_core)
set_target_properties(mc_runner PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(vr_bench app/cli/vr_bench/main.cpp)
target_link_libraries(vr_bench PRIVATE rw_core)
set_target_properties(vr_bench PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(greeks_runner app/cli/greeks_runner/main.cpp)
target_link_libraries(greeks_runner PRIVATE rw_core)
set_target_properties(greeks_runner PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# =======================
# GUI (Qt Widgets)
# =======================

# Auto MOC/UIC/RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Trouver Qt (essaie Qt6 puis Qt5)
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

add_executable(rw_gui
  app/gui/src/main.cpp
  app/gui/src/MainWindow.cpp
  app/gui/src/McWorker.cpp
  app/gui/ui/MainWindow.ui       # <-- très important : référencer le .ui
)


# Pour que AUTOUIC trouve MainWindow.ui quand MainWindow.cpp inclut "ui_MainWindow.h"
set_property(TARGET rw_gui PROPERTY AUTOUIC_SEARCH_PATHS
  ${PROJECT_SOURCE_DIR}/app/gui/ui
)

# Includes “locaux” de l’app
target_include_directories(rw_gui PRIVATE
  ${PROJECT_SOURCE_DIR}/app/gui/src
)

# Standard, liens
target_compile_features(rw_gui PRIVATE cxx_std_17)
target_link_libraries(rw_gui PRIVATE
  rw_core
  Qt${QT_VERSION_MAJOR}::Widgets
)

set_target_properties(rw_gui PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

