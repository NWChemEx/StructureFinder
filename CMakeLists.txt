# Copyright 2024 NWChemEx Community
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.14)

## Set Project and Version
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" VERSION)
project(structurefinder VERSION "${VERSION}" LANGUAGES CXX)

## Get CMaize
include(cmake/get_cmaize.cmake)

## Paths ##
set(${PROJECT_NAME}_INC_DIR "${CMAKE_CURRENT_LIST_DIR}/include/${PROJECT_NAME}")
set(${PROJECT_NAME}_SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/src/${PROJECT_NAME}")
set(${PROJECT_NAME}_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")
set(${PROJECT_NAME}_PYTHON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/python")
set(${PROJECT_NAME}_PYTHON_TEST_DIR "${${PROJECT_NAME}_TESTS_DIR}/python")

## Options ##
cmaize_option_list(
    BUILD_TESTING OFF "Should the tests be built?"
    BUILD_PYBIND11_PYBINDINGS ON "Use Pybind11 to build Python bindings?"
)

## Find or build dependencies ##
cmaize_find_or_build_dependency(
    simde
    URL github.com/NWChemEx/SimDE
    VERSION master
    BUILD_TARGET simde
    FIND_TARGET nwx::simde
    CMAKE_ARGS BUILD_TESTING=OFF
               BUILD_PYBIND11_PYBINDINGS=${BUILD_PYBIND11_PYBINDINGS}
)

cmaize_find_or_build_dependency(
    pyberny
    PACKAGE_MANAGER pip
)

# Builds C++ API documentation
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("${${PROJECT_NAME}_INC_DIR}" "${${PROJECT_NAME}_SRC_DIR}")

## Add libraries ##
add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE simde pyberny)

## Build tests ##
if("${BUILD_TESTING}")
    include(CTest)
    include(nwx_pybind11)

    ## Add Tests ##
    nwx_pybind11_tests(
        py_${PROJECT_NAME}
        "${${PROJECT_NAME}_PYTHON_TEST_DIR}/unit_tests/test_${PROJECT_NAME}.py"
        SUBMODULES simde chemist pluginplay parallelzone
    )

endif()

## Add package ##
install(
    DIRECTORY "${${PROJECT_NAME}_PYTHON_DIR}/${PROJECT_NAME}"
    DESTINATION "${NWX_MODULE_DIRECTORY}"
)
