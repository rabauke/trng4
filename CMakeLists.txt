cmake_minimum_required(VERSION 3.21)
project(trng4 VERSION 4.24 DESCRIPTION "pseudo random number generator library")

include(GNUInstallDirs)
enable_testing()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

option(TRNG_ENABLE_TESTS "Enable/Disable the compilation of the TRNG tests" ON)
option(TRNG_ENABLE_EXAMPLES "Enable/Disable the compilation of the TRNG examples" ON)

if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wmaybe-uninitialized")
    # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-shadow-field-in-constructor -Wno-float-equal -Wno-sign-conversion -Wno-implicit-float-conversion -Wno-implicit-int-float-conversion -Wno-global-constructors -Wno-exit-time-destructors")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_LIBCPP_DEBUG=1")
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive- /W4 /WX")
endif ()

if (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
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
if (TRNG_ENABLE_TESTS)
    add_subdirectory(external/Catch2)
    add_subdirectory(tests)
endif ()
if (TRNG_ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif ()

# Install package config files.
install(EXPORT "${targets_export_name}" DESTINATION "${pkg_config_install_dir}" NAMESPACE trng4::)
install(FILES "${project_pkg_config}" "${version_pkg_config}" DESTINATION "${pkg_config_install_dir}")

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
include(CPack)
