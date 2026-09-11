cmake_minimum_required(VERSION 3.16)
project(cns VERSION 1.0 LANGUAGES CXX
        DESCRIPTION "Celestial navigation: strapdown star fixes, dead reckoning")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED)

# ERFA: BSD-licensed fork of IAU SOFA.  Ubuntu/Debian: apt install liberfa-dev
find_path(ERFA_INCLUDE_DIR erfa.h)
find_library(ERFA_LIBRARY NAMES erfa)
if(NOT ERFA_INCLUDE_DIR OR NOT ERFA_LIBRARY)
  message(FATAL_ERROR "ERFA not found. Install with: sudo apt install liberfa-dev")
endif()

add_library(cns
  src/star_catalog.cpp
  src/star_catalog_data.cpp
  src/sky_model.cpp
  src/attitude.cpp
  src/orbit.cpp
  src/prior_net.cpp
  src/horizon.cpp
  src/imaging.cpp
  src/deadreckon.cpp
)
target_include_directories(cns PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${ERFA_INCLUDE_DIR}
)
target_link_libraries(cns PUBLIC Eigen3::Eigen ${ERFA_LIBRARY})
target_compile_options(cns PRIVATE -Wall -Wextra)

enable_testing()
add_executable(test_roundtrip test/test_roundtrip.cpp)
target_link_libraries(test_roundtrip PRIVATE cns)
add_test(NAME roundtrip COMMAND test_roundtrip)

add_executable(test_orbit test/test_orbit.cpp)
target_link_libraries(test_orbit PRIVATE cns)
add_test(NAME orbit COMMAND test_orbit)

add_executable(test_horizon test/test_horizon.cpp)
target_link_libraries(test_horizon PRIVATE cns)
add_test(NAME horizon COMMAND test_horizon)

add_executable(test_imaging test/test_imaging.cpp)
target_link_libraries(test_imaging PRIVATE cns)
add_test(NAME imaging COMMAND test_imaging)


# The deterministic end-to-end demonstration. Same numbers every run.
add_test(NAME scenario COMMAND transit_scenario)

# Milestone 3: analyse an ArduPilot SITL log.
#   python3 tools/dump_log.py flight.BIN > flight.csv
#   ./build/analyse_log flight.csv --imaging
add_executable(analyse_log tools/analyse_log.cpp)
target_link_libraries(analyse_log PRIVATE cns)

# Live MAVLink node -- needs the header-only MAVLink v2 library:
#   git clone --depth 1 https://github.com/mavlink/c_library_v2 $HOME/c_library_v2
#   cmake -S . -B build -DMAVLINK_DIR=$HOME/c_library_v2
find_path(MAVLINK_INCLUDE ardupilotmega/mavlink.h HINTS ${MAVLINK_DIR})
if(MAVLINK_INCLUDE)
  foreach(t celestial_node fake_sitl)
    if(${t} STREQUAL celestial_node)
      add_executable(${t} tools/celestial_node.cpp tools/mavlink_source.cpp tools/pipeline.cpp)
    else()
      add_executable(${t} tools/${t}.cpp)
    endif()
    target_link_libraries(${t} PRIVATE cns)
    target_include_directories(${t} PRIVATE ${MAVLINK_INCLUDE})
    target_compile_options(${t} PRIVATE -Wno-address-of-packed-member -Wno-cast-align)
  endforeach()
else()
  message(STATUS "MAVLink headers not found; skipping celestial_node and fake_sitl")
endif()

add_executable(transit_scenario examples/transit_scenario.cpp)
target_link_libraries(transit_scenario PRIVATE cns)

# OpenCV is OPTIONAL and used for exactly one thing: running an ONNX detection
# prior through cv::dnn, so a learned model does not drag in a second runtime.
# Everything else -- thresholding, connected components, centroiding, the
# matched filter, the mesh background -- is deliberately hand-written, because
# each is measured and documented and a library call would lose the reasoning.
# Without OpenCV the build is unchanged and PriorNet compiles to a stub.
# OpenCV 5 is the target: its DNN engine takes ONNX operator coverage from ~22%
# to over 80% and benchmarks competitively with ONNX Runtime on CPU. Only the
# stable subset is used here -- readNetFromONNX, blobFromImage, forward, resize
# -- which is unchanged from 4.x, so 4.x works too and is what CI has.
find_package(OpenCV QUIET COMPONENTS core imgproc dnn)
if(OpenCV_FOUND)
  message(STATUS "OpenCV ${OpenCV_VERSION} found: PriorNet enabled")
  if(OpenCV_VERSION VERSION_LESS 5.0)
    message(STATUS "  (OpenCV 5 preferred: new graph DNN engine, 80%+ ONNX)")
  endif()
  target_compile_definitions(cns PUBLIC CNS_WITH_OPENCV)
  target_include_directories(cns PUBLIC ${OpenCV_INCLUDE_DIRS})
  target_link_libraries(cns PUBLIC ${OpenCV_LIBS})
else()
  message(STATUS "OpenCV not found: PriorNet builds as a stub (this is fine)")
endif()

add_executable(gen_dataset tools/gen_dataset.cpp)
target_link_libraries(gen_dataset PRIVATE cns)

add_executable(bench tools/bench.cpp)
target_link_libraries(bench PRIVATE cns)

add_executable(visualise tools/visualise.cpp)
target_link_libraries(visualise PRIVATE cns)

