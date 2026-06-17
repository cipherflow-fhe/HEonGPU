# Copyright 2025-2026 Yanbin Li
# Licensed under the Apache License, Version 2.0, see LICENSE for details.
# SPDX-License-Identifier: Apache-2.0
# Developer: Yanbin Li

set(HEONGPU_FFT_SOURCES
    src/lib/fft.cu
)

add_library(heongpu_fft STATIC ${HEONGPU_FFT_SOURCES})
target_compile_features(heongpu_fft PUBLIC cxx_std_17)
target_include_directories(heongpu_fft
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
        $<INSTALL_INTERFACE:${INCLUDES_INSTALL_DIR}>
)
target_link_libraries(heongpu_fft PUBLIC fft)
set_target_properties(heongpu_fft PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CUDA_RUNTIME_LIBRARY Static
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

install(TARGETS heongpu_fft
  EXPORT ${HEonGPU_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/include/
  DESTINATION ${INCLUDES_INSTALL_DIR}
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.cuh"
    PATTERN "*.hpp"
)
