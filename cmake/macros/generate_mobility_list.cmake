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
    foreach(m IN LISTS matches)
      string(REGEX REPLACE "EXADIS_MOBILITY\\([ \t]*([^,]+)[ \t]*,[ \t]*([^\\)]+)[ \t]*\\)" "\\1" mname "${m}")
      string(REGEX REPLACE "EXADIS_MOBILITY\\([ \t]*([^,]+)[ \t]*,[ \t]*([^\\)]+)[ \t]*\\)" "\\2" malias "${m}")
      string(STRIP "${mname}" mname)
      string(STRIP "${malias}" malias)
      list(APPEND entries "    X(${mname}, ${malias})")
    endforeach()
  endif()
endforeach()

list(REMOVE_DUPLICATES includes)

# Write the output file
file(WRITE "${OUT_FILE}" "")
foreach(inc IN LISTS includes)
  file(APPEND "${OUT_FILE}" "${inc}\n")
endforeach()

file(APPEND "${OUT_FILE}" "\n#define EXADIS_MOBILITY_GLOBAL_LIST \\\n")

list(LENGTH entries num_entries)
if(num_entries EQUAL 0)
  file(APPEND "${OUT_FILE}" "    /* no mobilities discovered */\n")
else()
  math(EXPR last_index "${num_entries} - 1")
  foreach(idx RANGE 0 ${last_index})
    list(GET entries ${idx} entry)
    if(idx LESS last_index)
      set(suffix " \\")
    else()
      set(suffix "")
    endif()
    file(APPEND "${OUT_FILE}" "${entry}${suffix}\n")
  endforeach()
endif()
