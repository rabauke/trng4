find_package(Boost)

if (Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)

    add_executable(test_all test_all.cc test_engines.cc test_distributions.cc test_math.cc test_int_math.cc type_names.cc test_linear_algebra.cc)
    target_link_libraries(test_all PUBLIC trng4::trng4 Boost::unit_test_framework)
    if (WIN32 AND BUILD_SHARED_LIBS)
        add_custom_command(TARGET test_all POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_RUNTIME_DLLS:test_all>
                $<TARGET_FILE_DIR:test_all>
                COMMAND_EXPAND_LISTS)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
        target_compile_options(test_all PUBLIC "/bigobj")
    endif ()
    add_test(NAME run_test_all COMMAND test_all --log_level=test_suite CONFIGURATIONS Debug Release RelWithDebInfo)
    if (WIN32)
        set(exec_path_debug "PATH=${Boost_LIBRARY_DIR_DEBUG}$<SEMICOLON>$ENV{PATH}")
        set(exec_path_release "PATH=${Boost_LIBRARY_DIR_RELEASE}$<SEMICOLON>$ENV{PATH}")
        string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
        if (build_type STREQUAL debug)
            set_tests_properties(run_test_all PROPERTIES ENVIRONMENT "${exec_path_debug}")
        elseif (build_type STREQUAL release)
            set_tests_properties(run_test_all PROPERTIES ENVIRONMENT "${exec_path_release}")
        elseif (build_type STREQUAL relwithdebinfo)
            set_tests_properties(run_test_all PROPERTIES ENVIRONMENT "${exec_path_release}")
        endif ()
    endif ()
endif ()
