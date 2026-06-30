# HEonGPU vendors PhantomFHE for the optional PhantomNTT backend. Keep this
# wrapper unless PhantomFHE's CMake is made subproject-safe, whose standalone
# install rules assume PhantomFHE is the project source root.

option(PHANTOM_USE_CUDA_PTX "Use CUDA PTX Assembly in vendored PhantomFHE" ON)
if(PHANTOM_USE_CUDA_PTX)
    add_compile_definitions(PHANTOM_USE_CUDA_PTX)
endif()

set(PHANTOM_ROOT ${CMAKE_CURRENT_LIST_DIR}/phantom-fhe)
set(PHANTOM_INCLUDE_DIR ${PHANTOM_ROOT}/include)
file(GLOB_RECURSE PHANTOM_SOURCES CONFIGURE_DEPENDS
    ${PHANTOM_ROOT}/src/*.cu
)

add_library(Phantom SHARED ${PHANTOM_SOURCES})
target_compile_features(Phantom PUBLIC cxx_std_17 cuda_std_17)
target_include_directories(Phantom
    PUBLIC
        $<BUILD_INTERFACE:${PHANTOM_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include/phantom>
)
target_compile_options(Phantom
    PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--default-stream per-thread>
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-src-in-ptx>
)
set_target_properties(Phantom PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RUNTIME_LIBRARY Static
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
)

install(TARGETS Phantom
    EXPORT ${HEonGPU_TARGETS_EXPORT_NAME}
    RUNTIME DESTINATION ${RUNTIME_DESTINATION}
    LIBRARY DESTINATION ${LIBRARY_DESTINATION}
    ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
)

install(
    DIRECTORY ${PHANTOM_INCLUDE_DIR}/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/phantom
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.cuh"
        PATTERN "*.hpp"
)
