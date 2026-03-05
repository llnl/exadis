cmake_minimum_required(VERSION 3.16)

if(NOT DEFINED MOBILITY_HEADERS)
  message(FATAL_ERROR "MOBILITY_HEADERS not set")
endif()
if(NOT DEFINED OUT_FILE)
  message(FATAL_ERROR "OUT_FILE not set")
endif()

set(includes "")
set(entries "")

foreach(h IN LISTS MOBILITY_HEADERS)
  # Compute an include path
  file(TO_CMAKE_PATH "${h}" h_norm)
  if(DEFINED ROOT_DIR)
    file(RELATIVE_PATH inc_path "${ROOT_DIR}" "${h_norm}")
  else()
    get_filename_component(inc_path "${h_norm}" NAME)
  endif()

  # Collect all EXADIS_MOBILITY(Type, Alias) matches per file
  file(READ "${h_norm}" content)
  string(REGEX MATCHALL "EXADIS_MOBILITY\\([^\n]*\\)" matches "${content}")

  if(matches)
    list(APPEND includes "#include \"${inc_path}\"")
  endif()
endforeach()

list(REMOVE_DUPLICATES includes)

# Write the output file
file(WRITE "${OUT_FILE}" "")
foreach(inc IN LISTS includes)
  file(APPEND "${OUT_FILE}" "${inc}\n")
endforeach()
