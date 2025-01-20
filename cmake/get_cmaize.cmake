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

function(get_cmaize)

    if("${CMAIZE_VERSION}" STREQUAL "")
        set(CMAIZE_VERSION ca0f41e6d42829f9afb3c02f68ba4030fff6fee7 )
    endif()

    # Store whether we are building tests or not, then turn off the tests
    if(BUILD_TESTING)
        set(build_testing_old "${BUILD_TESTING}")
    endif()
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

    # Download CMakePP and bring it into scope
    include(FetchContent)
    FetchContent_Declare(
        cmaize
        GIT_REPOSITORY https://github.com/CMakePP/CMaize
        GIT_TAG ${CMAIZE_VERSION}
    )
    FetchContent_MakeAvailable(cmaize)

    # Restore the previous value, if set
    # Unset otherwise
    if(build_testing_old)
        set(BUILD_TESTING "${build_testing_old}" CACHE BOOL "" FORCE)
    else()
        unset(BUILD_TESTING CACHE)
    endif()
endfunction()

# Call the function we just wrote to get CMaize
get_cmaize()

# Include CMaize
include(cmaize/cmaize)
