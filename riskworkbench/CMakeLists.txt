# CMakeLists.txt (racine)
cmake_minimum_required(VERSION 3.22)
project(rw LANGUAGES CXX)

# Confort
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
  src/smile/smile.cpp
  src/pricing/analytic_bs.cpp
  src/pricing/mc_pricer.cpp
  src/pricing/greeks.cpp
  src/pricing/implied_vol.cpp  
  src/calibration/calib.cpp
  src/market/iv_surface.cpp  
  src/calibration/svi.cpp     
  src/calibration/calib_build_smile.cpp
  src/qc/qc.cpp
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
# IO marché (CSV)
# =======================
add_library(rw_io)

target_sources(rw_io PRIVATE
  src/io/market_csv.cpp         
  src/market/surface.cpp      
)

target_include_directories(rw_io PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(rw_io PUBLIC cxx_std_17)

if (MSVC)
  target_compile_options(rw_io PRIVATE /W4 /permissive-)
else()
  target_compile_options(rw_io PRIVATE -Wall -Wextra -Wpedantic)
endif()

# =======================
# CLIs (core)
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
# GUI (Qt Widgets + Charts)
# =======================
# Auto MOC/UIC/RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Trouver Qt (essaie Qt6 puis Qt5), une seule fois pour Widgets + Charts
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets Charts)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets Charts Concurrent)

add_executable(rw_gui
  app/gui/src/main.cpp
  app/gui/src/MainWindow.cpp
  app/gui/src/McWorker.cpp
  app/gui/ui/MainWindow.ui
  app/gui/src/CalibWorker.cpp
)

# Pour que AUTOUIC trouve MainWindow.ui
set_property(TARGET rw_gui PROPERTY AUTOUIC_SEARCH_PATHS
  ${PROJECT_SOURCE_DIR}/app/gui/ui
)

# Includes “locaux” de l’app
target_include_directories(rw_gui PRIVATE
  ${PROJECT_SOURCE_DIR}/app/gui/src
  ${PROJECT_SOURCE_DIR}/app/gui/widgets
)

target_compile_features(rw_gui PRIVATE cxx_std_17)

target_link_libraries(rw_gui PRIVATE
  rw_core
  rw_io
  Qt${QT_VERSION_MAJOR}::Widgets
  Qt${QT_VERSION_MAJOR}::Charts
  Qt${QT_VERSION_MAJOR}::Concurrent
)

set_target_properties(rw_gui PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
)

# Sanitizer pour la GUI
#if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  #target_compile_options(rw_gui PRIVATE -fsanitize=address -fno-omit-frame-pointer)
  #target_link_options(rw_gui PRIVATE -fsanitize=address)
#endif()

# =======================
# CLI test: lecture CSV
# =======================
add_executable(test_read_csv app/cli/market_csv/test_read_csv.cpp)
target_link_libraries(test_read_csv PRIVATE rw_io)
set_target_properties(test_read_csv PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# Sanitizer pour ce CLI (optionnel)
#if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  #target_compile_options(test_read_csv PRIVATE -fsanitize=address -fno-omit-frame-pointer)
  #target_link_options(test_read_csv PRIVATE -fsanitize=address)
#endif()

# =======================
# CLI: market_csv_info
# =======================
add_executable(market_csv_info app/cli/market_csv/market_csv_info.cpp)
target_link_libraries(market_csv_info PRIVATE rw_io)
set_target_properties(market_csv_info PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# =======================
# CLI test: IV solver
# =======================
add_executable(test_iv_solver app/cli/iv_solver/test_iv_solver.cpp)
target_link_libraries(test_iv_solver PRIVATE rw_core)
set_target_properties(test_iv_solver PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
#set_tests_properties(test_iv_solver PROPERTIES
#  ENVIRONMENT "ASAN_OPTIONS=abort_on_error=1:halt_on_error=1:print_stacktrace=1:symbolize=1")

#if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  #target_compile_options(test_iv_solver PRIVATE -fsanitize=address -fno-omit-frame-pointer)
  #target_link_options(test_iv_solver PRIVATE -fsanitize=address)
#endif()

add_executable(iv_from_csv app/cli/iv_solver/iv_from_csv.cpp)
target_link_libraries(iv_from_csv PRIVATE rw_core rw_io)
set_target_properties(iv_from_csv PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(iv_surface_info app/cli/iv_solver/iv_surface_info.cpp)
target_link_libraries(iv_surface_info PRIVATE rw_core rw_io)
set_target_properties(iv_surface_info PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(iv_export_and_qc app/cli/iv_solver/iv_export_and_qc.cpp)
target_link_libraries(iv_export_and_qc PRIVATE rw_core rw_io)
set_target_properties(iv_export_and_qc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# =======================
# CLI test: calib
# =======================

add_executable(calib_bs_global app/cli/calib/fit_bs_global.cpp)
target_link_libraries(calib_bs_global PRIVATE rw_core rw_io)
set_target_properties(calib_bs_global PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(make_synth_bs_market app/cli/calib/make_synth_bs_market.cpp)
target_link_libraries(make_synth_bs_market PRIVATE rw_core) # pas strictement nécessaire ici
set_target_properties(make_synth_bs_market PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(fit_svi_slices app/cli/svi/fit_svi.cpp)
target_link_libraries(fit_svi_slices PRIVATE rw_core rw_io)
set_target_properties(fit_svi_slices PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(smile_info app/cli/smile/smile_info.cpp)
target_link_libraries(smile_info PRIVATE rw_core rw_io)
set_target_properties(smile_info PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(qc_smile_fit app/cli/qc/qc_smile_fit.cpp)
target_link_libraries(qc_smile_fit PRIVATE rw_core rw_io)
set_target_properties(qc_smile_fit PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})



# =======================
# Tests CTest
# =======================
enable_testing()

add_test(NAME test_read_csv
  COMMAND ${BIN_DIR}/test_read_csv ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv
)

add_test(NAME market_csv_info_sample
  COMMAND ${BIN_DIR}/market_csv_info -f ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv -w -m
)

add_test(NAME test_iv_solver
  COMMAND ${BIN_DIR}/test_iv_solver
)

add_test(NAME iv_from_csv_sample
  COMMAND ${BIN_DIR}/iv_from_csv ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv)

add_test(NAME iv_surface_info_sample
  COMMAND ${BIN_DIR}/iv_surface_info ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv)

add_test(NAME iv_export_and_qc_sample
  COMMAND ${BIN_DIR}/iv_export_and_qc ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv)

add_test(NAME iv_surface_qc_thresholds
  COMMAND bash ${CMAKE_SOURCE_DIR}/tests/check_iv_surface.sh)
set_tests_properties(iv_surface_qc_thresholds PROPERTIES TIMEOUT 10)

add_test(NAME smile_build_and_check
  COMMAND ${BIN_DIR}/smile_info -f ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv --check)

add_test(NAME qc_smile_fit_sample
  COMMAND ${BIN_DIR}/qc_smile_fit -f ${CMAKE_SOURCE_DIR}/data/market_samples/sample_market.csv --thresh 0.5)
# ici seuil large (0.5%) pour l’exemple; ajuste selon tes données
