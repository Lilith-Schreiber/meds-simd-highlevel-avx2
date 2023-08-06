cmake_minimum_required(VERSION 3.10)

project(MEDS_simd)

find_package(OpenSSL REQUIRED)

option(DEBUG "Build with debug flags" OFF)
if(DEBUG)
  message(STATUS "Building with debug flags")
  set(CMAKE_C_FLAGS "-g -Wall -mavx512f")
  add_compile_definitions(DEBUG)
else()
  message(STATUS "Building with release flags")
  set(CMAKE_C_FLAGS "-O3 -Wall -mavx512f")
endif()

# set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

set(NIST_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/NIST)

set(LIBS
    bitstream
    fips202
    matrixmod
    measure
    meds
    randombytes
    seed
    util

    bitstream_batch
    vec_x16
    util_x16
    matrixmod_x16
    fips202_x8
    vec_x8)

# ##############################################################################
# Dump parameters
# ##############################################################################

# Create paramters header
execute_process(
  COMMAND bash "-c" "./params.py -p > ${PROJECT_INCLUDE_DIR}/params.h"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# List all parameter sets
execute_process(
  COMMAND ./params.py -l
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE PARAMS)
string(REPLACE "\n" ";" PARAMS ${PARAMS})

foreach(param ${PARAMS})
  message(STATUS ${param})

  # ############################################################################
  # Create configuration specified headers
  # ############################################################################

  set(API_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/${param})

  execute_process(COMMAND bash "-c" "mkdir -p ${API_INCLUDE_DIR}"
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  # Create api header
  execute_process(
    COMMAND bash "-c" "./params.py -a ${param} > ${API_INCLUDE_DIR}/api.h"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  set(INCLUDE_DIRECTORIES ${PROJECT_INCLUDE_DIR} ${NIST_INCLUDE_DIR}
                          ${API_INCLUDE_DIR})

  # ############################################################################
  # Build libraries
  # ############################################################################

  set(LIBRARIES)

  foreach(lib ${LIBS})
    set(LIBRARY ${lib}_${param})

    add_library(${LIBRARY} STATIC ${PROJECT_SOURCE_DIR}/${lib}.c)
    target_include_directories(${LIBRARY} PRIVATE ${INCLUDE_DIRECTORIES})
    target_compile_definitions(${LIBRARY} PRIVATE ${param})

    set_target_properties(${LIBRARY} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
                                                ${CMAKE_BINARY_DIR}/${param})

    list(APPEND LIBRARIES ${LIBRARY})
  endforeach()

  target_link_libraries(
    meds_${param}
    bitstream_${param}
    matrixmod_${param}
    measure_${param}
    randombytes_${param}
    seed_${param}
    util_${param}
    
    vec_x16_${param}
    util_x16_${param}
    matrixmod_x16_${param}
    bitstream_batch_${param})

  target_link_libraries(seed_${param} fips202_${param})
  target_link_libraries(util_${param} matrixmod_${param} fips202_${param})
  target_link_libraries(randombytes_${param} ${OPENSSL_LIBRARIES})

  target_link_libraries(matrixmod_x16_${param} vec_x16_${param})
  target_link_libraries(util_x16_${param} vec_x16_${param} matrixmod_x16_${param})
  target_link_libraries(fips202_x8_${param} vec_x8_${param})

  # ############################################################################
  # Build test
  # ############################################################################

  set(TEST test_${param})

  add_executable(${TEST} ${PROJECT_SOURCE_DIR}/test.c)
  target_include_directories(${TEST} PRIVATE ${INCLUDE_DIRECTORIES})
  target_link_libraries(${TEST} ${LIBRARIES})
  target_compile_definitions(${TEST} PRIVATE ${param})

  set_target_properties(
    ${TEST} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${param}
                       OUTPUT_NAME test)

  if(param STREQUAL "toy")
    set(RUN_test RUN)
  else()
    set(RUN_test RUN_${param})
  endif()
  add_custom_target(
    ${RUN_test}
    COMMAND ./${param}/test
    DEPENDS ${TEST})

  list(APPEND RUN_ALL_DEPENDS ${RUN_test})

  # ############################################################################
  # Build KAT_test
  # ############################################################################

  set(KAT kat_${param})

  add_executable(${KAT} ${PROJECT_SOURCE_DIR}/KAT_test.c)
  target_include_directories(${KAT} PRIVATE ${INCLUDE_DIRECTORIES})
  target_link_libraries(${KAT} ${LIBRARIES})
  target_compile_definitions(${KAT} PRIVATE ${param})

  set_target_properties(
    ${KAT}
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${param}
               RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${param}
               OUTPUT_NAME KAT_test)

  if(param STREQUAL "toy")
    set(RUN_KAT KAT)
  else()
    set(RUN_KAT KAT_${param})
  endif()
  add_custom_target(
    ${RUN_KAT}
    COMMAND ./${param}/KAT_test
    DEPENDS ${KAT})

  list(APPEND KAT_ALL_DEPENDS ${RUN_KAT})

  # ############################################################################
  # Build bench
  # ############################################################################

  set(BENCH bench_${param})

  add_executable(${BENCH} ${PROJECT_SOURCE_DIR}/bench.c)
  target_include_directories(${BENCH} PRIVATE ${INCLUDE_DIRECTORIES})
  target_link_libraries(${BENCH} ${LIBRARIES})
  target_compile_definitions(${BENCH} PRIVATE ${param})

  set_target_properties(
    ${BENCH}
    PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${param}
               RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${param}
               OUTPUT_NAME bench)

  if(param STREQUAL "toy")
    set(RUN_bench BENCH)
  else()
    set(RUN_bench BENCH_${param})
  endif()
  add_custom_target(
    ${RUN_bench}
    COMMAND bash "-c" "./${param}/bench" "|" "${CMAKE_SOURCE_DIR}/proc_bench.py"
            ">>" "bench.txt"
    DEPENDS ${BENCH})

  list(APPEND BENCH_ALL_DEPENDS ${RUN_bench})

endforeach()

add_custom_target(RUN_ALL DEPENDS ${RUN_ALL_DEPENDS})

add_custom_target(KAT_ALL DEPENDS ${KAT_ALL_DEPENDS})

add_custom_target(BENCH_ALL DEPENDS ${BENCH_ALL_DEPENDS})
