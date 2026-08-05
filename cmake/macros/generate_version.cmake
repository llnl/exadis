function(exadis_generate_version_header TARGET_NAME)
    # Define the output file path in the target's binary build folder
    set(VERSION_HDR_OUT "${CMAKE_CURRENT_BINARY_DIR}/exadis_version.h")

    # Generate the header contents dynamically using string contents
    file(WRITE "${VERSION_HDR_OUT}"
"#ifndef EXADIS_VERSION_H\n"
"#define EXADIS_VERSION_H\n\n"
"#define EXADIS_VERSION \"${PROJECT_VERSION}\"\n"
"#define EXADIS_VERSION_MAJOR  ${PROJECT_VERSION_MAJOR}\n"
"#define EXADIS_VERSION_MINOR  ${PROJECT_VERSION_MINOR}\n"
"#define EXADIS_VERSION_PATCH  ${PROJECT_VERSION_PATCH}\n\n"
"#endif // EXADIS_VERSION_H\n"
    )

    # Inform CMake that this file is dynamically created
    set_source_files_properties("${VERSION_HDR_OUT}" PROPERTIES GENERATED TRUE)

    # Attach the file safely to your target sources list
    target_sources(${TARGET_NAME} PRIVATE "${VERSION_HDR_OUT}")

    # Track it globally so our root install loop copies it automatically
    set_property(TARGET ${TARGET_NAME} APPEND PROPERTY 
        EXADIS_GENERATED_HEADERS "${VERSION_HDR_OUT}"
    )
endfunction()