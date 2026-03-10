
SET(SOURCES
)

target_sources(exadis PRIVATE ${SOURCES})

# Set directories to scan to register mobilities
set(MOBILITY_SEARCH_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Gather headers for build dependencies
set(MOBILITY_HEADERS "")
foreach(d IN LISTS MOBILITY_SEARCH_DIRS)
    file(GLOB tmp CONFIGURE_DEPENDS "${d}/mobility_*.h")
    list(APPEND MOBILITY_HEADERS ${tmp})
endforeach()
list(REMOVE_DUPLICATES MOBILITY_HEADERS)

# Generate global mobility header at build time using CMake script mode
set(MOBILITY_LIST_OUT "${CMAKE_CURRENT_BINARY_DIR}/mobility_global_list.h")
add_custom_command(
    OUTPUT "${MOBILITY_LIST_OUT}"
    COMMAND ${CMAKE_COMMAND}
            "-DMOBILITY_HEADERS=${MOBILITY_HEADERS}"
            "-DROOT_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
            "-DOUT_FILE=${MOBILITY_LIST_OUT}"
            -P ${EXADIS_MACROS_DIR}/generate_mobility_list.cmake
    DEPENDS ${MOBILITY_HEADERS}
    COMMENT "Generating mobility_global_list.h from mobility_*.h registrations"
    VERBATIM
)

add_custom_target(generate_mobility_list DEPENDS "${MOBILITY_LIST_OUT}")

add_dependencies(exadis generate_mobility_list)
