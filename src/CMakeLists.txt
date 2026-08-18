cmake_minimum_required(VERSION 3.14)

set(SOURCES
    driver.cpp
    generate.cpp
    network.cpp
    oprec.cpp
    system.cpp
    utils.cpp
)

# Use SHARED library with SYCL
# Static libraries can cause device context corruption with SYCL when Kokkos::SharedSpace
if(Kokkos_ENABLE_SYCL AND NOT EXADIS_SYCL_USE_STATIC)
    add_library(exadis SHARED ${SOURCES})
    message(STATUS "Building ExaDiS as SHARED library for SYCL (set EXADIS_SYCL_USE_STATIC=ON to override)")
else()
    add_library(exadis ${SOURCES})
endif()

foreach(SUBDIR ${EXADIS_SUBDIRECTORIES})
    add_subdirectory(${SUBDIR})
endforeach()

# Generate version number header
include("${EXADIS_MACROS_DIR}/generate_version.cmake")
exadis_generate_version_header(exadis)

# Collect everything that was registered globally in the subdirectories
get_target_property(ALL_GENERATED_HEADERS exadis EXADIS_GENERATED_HEADERS)

if(ALL_GENERATED_HEADERS)
    # Apply the dynamic generation tag at the ROOT scope level
    set_source_files_properties(${ALL_GENERATED_HEADERS} PROPERTIES GENERATED TRUE)
    # Append the absolute paths to the main target
    target_sources(exadis PRIVATE ${ALL_GENERATED_HEADERS})
endif()

# Extracts exadis subdirectory names (e.g. "force_types") from the absolute paths
set(EXADIS_BUILD_INCLUDES "")
set(EXADIS_INSTALL_INCLUDES "")
foreach(SUBDIR ${EXADIS_SUBDIRECTORIES})
    get_filename_component(FOLDER_NAME "${SUBDIR}" NAME)
    list(APPEND EXADIS_BUILD_INCLUDES 
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${FOLDER_NAME}>"
    )
    list(APPEND EXADIS_INSTALL_INCLUDES 
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/exadis/${FOLDER_NAME}>"
    )
endforeach()

target_include_directories(exadis 
    PUBLIC
        # Paths used during the build of this repository
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        ${EXADIS_BUILD_INCLUDES}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/force_types>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mobility_types>
        # Path used by external users after they install exadis
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/exadis>
        ${EXADIS_INSTALL_INCLUDES}
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/kokkos>
    PRIVATE
        ${FFT_INC_DIR}
)

target_link_libraries(exadis
    PUBLIC
        kokkos
        ${EXADIS_LINK_OPTIONS}
    PRIVATE
        ${FFT_LIB}
)

# Install rules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install the compiled libraries and source headers
install(TARGETS exadis kokkos kokkoscore kokkoscontainers kokkosalgorithms kokkossimd LIBDL
    EXPORT exadisTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ./
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/exadis
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

# Install all dynamically generated headers from the subdirectories
if(ALL_GENERATED_HEADERS)
    foreach(GEN_HEADER ${ALL_GENERATED_HEADERS})
        # Extract the relative path from the build folder to preserve layout
        file(RELATIVE_PATH REL_PATH "${CMAKE_CURRENT_BINARY_DIR}" "${GEN_HEADER}")
        get_filename_component(REL_DIR "${REL_PATH}" DIRECTORY)
        
        # Install each generated header to its matching folder in the install tree
        install(FILES "${GEN_HEADER}"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/exadis/${REL_DIR}
        )
    endforeach()
endif()


# Export the CMake target configurations so one can link against it
set(EXADIS_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/exadis")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/exadisConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(EXPORT exadisTargets
    FILE exadisConfig.cmake
    NAMESPACE exadis::
    DESTINATION "${EXADIS_CMAKE_INSTALL_DIR}"
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/exadisConfigVersion.cmake"
    DESTINATION "${EXADIS_CMAKE_INSTALL_DIR}"
)
