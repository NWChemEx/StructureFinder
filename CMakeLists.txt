cmake_minimum_required(VERSION 3.14)
set(VERSION 1.0.0) #TODO: get from git
project(structurefinder VERSION "${VERSION}" LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
    nwx_cmake
    GIT_REPOSITORY https://github.com/NWChemEx-Project/NWXCMake
)
FetchContent_MakeAvailable(nwx_cmake)
list(APPEND CMAKE_MODULE_PATH "${nwx_cmake_SOURCE_DIR}/cmake")

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(get_cmaize)
include(nwx_versions)

### Options ###
option(BUILD_TESTING "Should we build the tests?" OFF)

## Build StructureFinder's dependencies ##
cpp_find_or_build_dependency(
    nwchemex
    URL github.com/NWChemEx-Project/NWChemEx
    PRIVATE TRUE
    VERSION master
    BUILD_TARGET nwchemex
    FIND_TARGET nwx::nwchemex
    CMAKE_ARGS BUILD_TESTING=OFF
)

if("${BUILD_TESTING}")
    include(CTest)
    include(nwx_pybind11)
    
# override
function(nwx_pybind11_tests npt_name npt_driver)
    if(NOT "${BUILD_PYBIND11_PYBINDINGS}")
        return()
    endif()

    set(_npt_options "")
    set(_npt_one_val "")
    set(_npt_lists SUBMODULES)
    cmake_parse_arguments(
        "_npt" "${_npt_options}" "${_npt_one_val}" "${_npt_lists}" ${ARGN}
    )

    if("${BUILD_TESTING}")
        include(CTest)
        find_package(Python COMPONENTS Interpreter REQUIRED)

        # Build the PYTHONPATH for the test
        set(_npt_py_path "PYTHONPATH=${NWX_MODULE_DIRECTORY}")
        set(_npt_py_path "${_npt_py_path}:${CMAKE_BINARY_DIR}")
        foreach(_npt_submod ${_npt_SUBMODULES})
            set(_npt_dep_dir "${CMAKE_BINARY_DIR}/_deps/${_npt_submod}-build")
            set(_npt_py_path "${_npt_py_path}:${_npt_dep_dir}")
        endforeach()
        if(NOT "${NWX_PYTHON_EXTERNALS}" STREQUAL "")
            set(_npt_py_path "${_npt_py_path}:${NWX_PYTHON_EXTERNALS}")
        endif()

        add_test(
            NAME "${npt_name}"
            COMMAND "${Python_EXECUTABLE}" "${npt_driver}"
        )

        set_tests_properties(
            "${npt_name}"
            PROPERTIES ENVIRONMENT "${_npt_py_path}"
        )
    endif()
endfunction()

#
    
    set(PYTHON_TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/tests/python")

    message("calling nwx_pybind11_tests")
    nwx_pybind11_tests(
        pyStructureFinder_unit_tests
        "${PYTHON_TEST_DIR}/unit_tests/test_structurefinder.py"
        SUBMODULES nwchemex simde chemist pluginplay parallelzone
    )
endif()

install(
    DIRECTORY "${python_src_directory}/structurefinder"
    DESTINATION "${NWX_MODULE_DIRECTORY}"
)
