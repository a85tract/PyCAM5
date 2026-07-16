cmake_minimum_required(VERSION 3.20)

# Configure the dependency gates even on a login node without a Fortran
# compiler.  Fortran is enabled only when an actual standalone source build is
# requested.
project(PIAtmStandalonePhysics LANGUAGES NONE)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
include(CTest)

get_filename_component(TPHYS_CAM_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
set(
  TPHYS_STANDALONE_STATUS_FILE
  "${TPHYS_CAM_ROOT}/src/atmos_phys/standalone/standalone_status.json"
  CACHE FILEPATH
  "Phase-2 standalone ledger"
)
set(
  TPHYS_STANDALONE_PROGRESS_FILE
  "${TPHYS_CAM_ROOT}/src/atmos_phys/standalone/STANDALONE_PROGRESS.md"
  CACHE FILEPATH
  "Generated phase-2 progress Markdown"
)
set(
  TPHYS_STANDALONE_SCOPE
  "completed-only"
  CACHE STRING
  "Dependency audit scope: completed-only or all"
)
set_property(CACHE TPHYS_STANDALONE_SCOPE PROPERTY STRINGS completed-only all)
set(
  TPHYS_STANDALONE_MODE
  "both"
  CACHE STRING
  "Dependency audit mode: direct, transitive, or both"
)
set_property(CACHE TPHYS_STANDALONE_MODE PROPERTY STRINGS direct transitive both)

option(
  TPHYS_STANDALONE_BUILD_SOURCES
  "Compile the selected transitive source closure as a standalone library"
  OFF
)
set(
  TPHYS_STANDALONE_MODULE_INCLUDE_DIRS
  ""
  CACHE STRING
  "Semicolon-separated standalone dependency module directories; CAM EXEROOT paths are forbidden"
)
set(
  TPHYS_STANDALONE_EXTERNAL_LIBRARIES
  ""
  CACHE STRING
  "Semicolon-separated external CMake library targets or link items"
)

if(NOT TPHYS_STANDALONE_SCOPE MATCHES "^(completed-only|all)$")
  message(FATAL_ERROR "TPHYS_STANDALONE_SCOPE must be completed-only or all")
endif()
if(NOT TPHYS_STANDALONE_MODE MATCHES "^(direct|transitive|both)$")
  message(FATAL_ERROR "TPHYS_STANDALONE_MODE must be direct, transitive, or both")
endif()

set(_tracker "${TPHYS_CAM_ROOT}/tools/tphys_standalone_progress.py")
set(_dependency_checker "${TPHYS_CAM_ROOT}/tools/check_atmos_phys_standalone.py")
set(_python_command "${CMAKE_COMMAND}" -E env PYTHONDONTWRITEBYTECODE=1 "${Python3_EXECUTABLE}")

add_custom_target(
  check-standalone
  COMMAND
    ${_python_command} "${_tracker}" check
    --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
    --progress-file "${TPHYS_STANDALONE_PROGRESS_FILE}"
  COMMAND
    ${_python_command} "${_dependency_checker}"
    --repo-root "${TPHYS_CAM_ROOT}"
    --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
    --scope "${TPHYS_STANDALONE_SCOPE}"
    --mode "${TPHYS_STANDALONE_MODE}"
  USES_TERMINAL
  VERBATIM
  COMMENT "Checking phase-2 ledger and ${TPHYS_STANDALONE_SCOPE} standalone dependencies"
)

if(BUILD_TESTING)
  add_test(
    NAME tphys-standalone-ledger
    COMMAND
      ${_python_command} "${_tracker}" check
      --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
      --progress-file "${TPHYS_STANDALONE_PROGRESS_FILE}"
  )
  add_test(
    NAME tphys-standalone-dependencies
    COMMAND
      ${_python_command} "${_dependency_checker}"
      --repo-root "${TPHYS_CAM_ROOT}"
      --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
      --scope "${TPHYS_STANDALONE_SCOPE}"
      --mode "${TPHYS_STANDALONE_MODE}"
  )
endif()

add_custom_target(
  check-standalone-all
  COMMAND
    ${_python_command} "${_tracker}" check
    --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
    --progress-file "${TPHYS_STANDALONE_PROGRESS_FILE}"
  COMMAND
    ${_python_command} "${_dependency_checker}"
    --repo-root "${TPHYS_CAM_ROOT}"
    --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
    --scope all
    --mode both
  USES_TERMINAL
  VERBATIM
  COMMENT "Checking direct and transitive closure for all 36 standalone entries"
)

function(tphys_reject_cam_exeroot_module_dir raw_directory)
  if(raw_directory STREQUAL "")
    return()
  endif()
  get_filename_component(_directory "${raw_directory}" ABSOLUTE)
  string(REPLACE "\\" "/" _normalized "${_directory}")
  if(_normalized MATCHES "(^|/)(bld|Buildconf)(/|$)")
    message(
      FATAL_ERROR
      "CAM/CIME build module directory is forbidden: ${_directory}. "
      "Build every atmos_phys module in this standalone CMake tree."
    )
  endif()
  set(_cursor "${_directory}")
  foreach(_level RANGE 0 5)
    if(EXISTS "${_cursor}/cesm.exe")
      message(
        FATAL_ERROR
        "Module directory ${_directory} is inside CAM EXEROOT ${_cursor}; borrowing its .mod files is forbidden."
      )
    endif()
    get_filename_component(_parent "${_cursor}" DIRECTORY)
    if(_parent STREQUAL _cursor)
      break()
    endif()
    set(_cursor "${_parent}")
  endforeach()
endfunction()

foreach(_module_directory IN LISTS TPHYS_STANDALONE_MODULE_INCLUDE_DIRS)
  tphys_reject_cam_exeroot_module_dir("${_module_directory}")
endforeach()

if(TPHYS_STANDALONE_BUILD_SOURCES)
  if(TPHYS_STANDALONE_MODE STREQUAL "direct")
    message(FATAL_ERROR "A source build requires transitive or both dependency mode")
  endif()

  execute_process(
    COMMAND
      ${_python_command} "${_dependency_checker}"
      --repo-root "${TPHYS_CAM_ROOT}"
      --status-file "${TPHYS_STANDALONE_STATUS_FILE}"
      --scope "${TPHYS_STANDALONE_SCOPE}"
      --mode transitive
      --list-sources
    RESULT_VARIABLE _source_result
    OUTPUT_VARIABLE _source_output
    ERROR_VARIABLE _source_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _source_result EQUAL 0)
    message(FATAL_ERROR "Standalone dependency closure failed:\n${_source_error}")
  endif()

  if(_source_output STREQUAL "")
    message(STATUS "No ${TPHYS_STANDALONE_SCOPE} source entries are ready to build")
    add_custom_target(tphys-standalone-schemes)
  else()
    string(REPLACE "\n" ";" _standalone_sources "${_source_output}")
    enable_language(Fortran)
    add_library(tphys-standalone-schemes STATIC ${_standalone_sources})

    # This is the only module output directory for atmos_phys sources.  Raw
    # external include directories are accepted only after the EXEROOT guard
    # above; no CAM build tree is searched implicitly.
    set_target_properties(
      tphys-standalone-schemes
      PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
        POSITION_INDEPENDENT_CODE ON
    )
    target_include_directories(
      tphys-standalone-schemes
      PRIVATE
        "${CMAKE_CURRENT_BINARY_DIR}/modules"
        ${TPHYS_STANDALONE_MODULE_INCLUDE_DIRS}
    )
    if(TPHYS_STANDALONE_EXTERNAL_LIBRARIES)
      target_link_libraries(
        tphys-standalone-schemes
        PRIVATE
          ${TPHYS_STANDALONE_EXTERNAL_LIBRARIES}
      )
    endif()
    add_dependencies(tphys-standalone-schemes check-standalone)
  endif()
endif()
