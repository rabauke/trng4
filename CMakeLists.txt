cmake_minimum_required(VERSION 3.10)
project(trng4 VERSION 4.24 DESCRIPTION "pseudo random number generator library")

include(GNUInstallDirs)
enable_testing()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_LIBCPP_DEBUG=1")
endif ()

# Configure package config files.
set(pkg_config_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
set(project_pkg_config "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake")
set(version_pkg_config "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${version_pkg_config}" COMPATIBILITY SameMajorVersion)
configure_package_config_file("cmake/Config.cmake.in" "${project_pkg_config}" INSTALL_DESTINATION "${pkg_config_install_dir}")

add_subdirectory(trng)
add_subdirectory(examples)
add_subdirectory(tests)

# Install package config files.
install(FILES "${project_pkg_config}" "${version_pkg_config}" DESTINATION "${pkg_config_install_dir}")
install(EXPORT "${targets_export_name}" DESTINATION "${pkg_config_install_dir}")

find_program(RPMBUILD_FOUND rpmbuild)

if (CMAKE_SYSTEM_NAME MATCHES Linux)
    if (RPMBUILD_FOUND)
        set(CPACK_GENERATOR "TGZ;DEB;RPM")
    else ()
        set(CPACK_GENERATOR "TGZ;DEB")
    endif ()
else ()
    set(CPACK_GENERATOR "ZIP")
endif ()

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Heiko Bauke")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://www.numbercrunch.de/trng/")
INCLUDE(CPack)
