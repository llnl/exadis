
SET(SOURCES
)

target_sources(exadis PRIVATE ${SOURCES})

# Set directories to scan to register forces
set(FORCE_SEARCH_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Gather headers for build dependencies
set(FORCE_HEADERS "")
foreach(d IN LISTS FORCE_SEARCH_DIRS)
    file(GLOB tmp CONFIGURE_DEPENDS "${d}/force_*.h")
    list(APPEND FORCE_HEADERS ${tmp})
endforeach()
list(REMOVE_DUPLICATES FORCE_HEADERS)

# Generate global force header at build time using CMake script mode
set(FORCE_LIST_OUT "${CMAKE_CURRENT_BINARY_DIR}/force_global_list.h")
add_custom_command(
    OUTPUT "${FORCE_LIST_OUT}"
    COMMAND ${CMAKE_COMMAND}
            "-DFORCE_HEADERS=${FORCE_HEADERS}"
            "-DROOT_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
            "-DOUT_FILE=${FORCE_LIST_OUT}"
            -P ${EXADIS_MACROS_DIR}/generate_force_list.cmake
    DEPENDS ${FORCE_HEADERS}
    COMMENT "Generating force_global_list.h from force_*.h registrations"
    VERBATIM
)

add_custom_target(generate_force_list DEPENDS "${FORCE_LIST_OUT}")

add_dependencies(exadis generate_force_list)
