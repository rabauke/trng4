find_package(Boost)

if (Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)
    include_directories(..)
    link_libraries(trng4)
    link_directories(${PROJECT_BINARY_DIR}/trng)

    add_executable(test_all test_all.cc test_engines.cc test_distributions.cc test_math.cc type_names.cc test_linear_algebra.cc)
    target_link_libraries(test_all PUBLIC Boost::unit_test_framework)
endif ()
