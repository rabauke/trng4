find_package(Boost)

if (Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)

    add_executable(test_all test_all.cc test_engines.cc test_distributions.cc test_math.cc test_int_math.cc type_names.cc test_linear_algebra.cc)
    target_link_libraries(test_all PUBLIC trng4 Boost::unit_test_framework)
    add_test(NAME run_test_all COMMAND test_all --log_level=test_suite)
endif ()
