########################################################################
# PothosComms: DSP and communications blocks for Pothos.
# Builds against an installed Pothos core (find_package(Pothos)).
# The modules use only the public API of the core. Their third party
# libraries (json, xsimd) come from find_package: from the system,
# or from a package manager that gives a toolchain file.
########################################################################
cmake_minimum_required(VERSION 3.25)
project(PothosComms C CXX)

option(COMMS_BUILD_TESTS "Register the CTest tests of the comms blocks" ${PROJECT_IS_TOP_LEVEL})
option(COMMS_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (GCC and Clang)" OFF)

find_package(Pothos CONFIG REQUIRED)
include(PothosConfigSIMD) #pothos_multiarch, PothosGenerateSIMDSources

#the core loads modules from its own prefix (POTHOS_ROOT): install there
#when the user does not give CMAKE_INSTALL_PREFIX
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${POTHOS_ROOT}" CACHE PATH "Install prefix (default: the prefix of the Pothos core)" FORCE)
    message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX} (the prefix of the Pothos core)")
endif()

########################################################################
# Third party libraries
########################################################################
set(comms_missing_dependencies)

find_package(nlohmann_json 3.11 CONFIG)
if(NOT nlohmann_json_FOUND)
    list(APPEND comms_missing_dependencies "nlohmann_json (>= 3.11)")
endif()

find_package(xsimd 14.2 CONFIG)
if(NOT xsimd_FOUND)
    list(APPEND comms_missing_dependencies "xsimd (>= 14.2)")
endif()

if(comms_missing_dependencies)
    list(JOIN comms_missing_dependencies ", " comms_missing_text)
    message(FATAL_ERROR
        "PothosComms: missing build dependencies: ${comms_missing_text}.\n"
        "Install the system packages (for example: nlohmann-json, xsimd), "
        "or install them with a package manager and give its toolchain file "
        "(-DCMAKE_TOOLCHAIN_FILE=...) or CMAKE_PREFIX_PATH. Refer to README.md.")
endif()

add_library(comms_json INTERFACE)
target_link_libraries(comms_json INTERFACE nlohmann_json::nlohmann_json)

add_library(comms_xsimd INTERFACE)
target_link_libraries(comms_xsimd INTERFACE xsimd)

########################################################################
# Modules
########################################################################
#comms_warnings(<target>): the common warnings (GCC and Clang, as the core),
#errors with COMMS_WARNINGS_AS_ERRORS (CI)
function(comms_warnings target)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra)
        if(COMMS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

#pothos_comms_module(TARGET <name> SOURCES <sources...>
#    [LIBRARIES <libraries...>] [TESTS <self test names...>])
#builds and installs a comms module, TESTS are names under /comms/tests
#run one by one by the comms.selftests test
function(pothos_comms_module)
    cmake_parse_arguments(COMMS "" "TARGET" "SOURCES;LIBRARIES;TESTS" ${ARGN})
    POTHOS_MODULE_UTIL(
        TARGET ${COMMS_TARGET}
        SOURCES ${COMMS_SOURCES}
        LIBRARIES comms_json ${COMMS_LIBRARIES}
        DESTINATION comms
        ENABLE_DOCS)
    set_target_properties(${COMMS_TARGET} PROPERTIES CXX_EXTENSIONS OFF)
    #the repository root for common/
    target_include_directories(${COMMS_TARGET} PRIVATE ${PROJECT_SOURCE_DIR})
    comms_warnings(${COMMS_TARGET})
    foreach(test IN LISTS COMMS_TESTS)
        set_property(GLOBAL APPEND PROPERTY POTHOS_COMMS_TESTS /comms/tests/${test})
    endforeach()
endfunction()

#pothos_comms_static_library(<target> <sources...>)
#a static helper library for the modules (PIC: it goes into modules)
function(pothos_comms_static_library target)
    add_library(${target} STATIC ${ARGN})
    set_target_properties(${target} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_EXTENSIONS OFF)
    target_link_libraries(${target} PRIVATE Pothos::Pothos)
    comms_warnings(${target})
endfunction()

add_subdirectory(functions)
add_subdirectory(demod)
add_subdirectory(digital)
add_subdirectory(fft)
add_subdirectory(mac)
add_subdirectory(math)
add_subdirectory(utility)
add_subdirectory(waveform)
add_subdirectory(filter/design)
add_subdirectory(filter)
add_subdirectory(window)

########################################################################
# Tests: the modules are staged into a prefix inside the build tree and
# loaded with POTHOS_PLUGIN_PATH by the PothosUtil of the installed core.
# The tests use blocks of PothosBlocks: install it into the core prefix.
########################################################################
if(COMMS_BUILD_TESTS)
    enable_testing()
    set(COMMS_STAGE_DIR ${PROJECT_BINARY_DIR}/stage)
    set(COMMS_SELFTEST_TIMEOUT 900 CACHE STRING "Timeout in seconds for the self-tests")

    add_test(NAME comms.stage
        COMMAND ${CMAKE_COMMAND} --install ${PROJECT_BINARY_DIR} --prefix ${COMMS_STAGE_DIR} --config $<CONFIG>)
    set_tests_properties(comms.stage PROPERTIES FIXTURES_SETUP comms_stage)

    get_property(comms_tests GLOBAL PROPERTY POTHOS_COMMS_TESTS)
    add_test(NAME comms.selftests
        COMMAND ${CMAKE_COMMAND}
            -DPREFIX=${POTHOS_ROOT}
            -DTEST_HOME_DIR=${PROJECT_BINARY_DIR}/test-home
            -DBINDIR=bin
            -DLIBDIR=${CMAKE_INSTALL_LIBDIR}
            -DPLUGIN_PATH=${COMMS_STAGE_DIR}/${CMAKE_INSTALL_LIBDIR}/Pothos/modules${POTHOS_ABI_VERSION}
            "-DTEST_PATHS=${comms_tests}"
            -P ${PROJECT_SOURCE_DIR}/cmake/PothosRunSelfTests.cmake)
    set_tests_properties(comms.selftests PROPERTIES
        FIXTURES_REQUIRED comms_stage
        TIMEOUT ${COMMS_SELFTEST_TIMEOUT})
endif()
