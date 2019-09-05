find_package(Boost)

if (Boost_FOUND)
    set(Boost_USE_STATIC_LIBS OFF)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(..)
    link_libraries(trng4)
    link_directories(${PROJECT_BINARY_DIR}/trng)

    add_executable(test_engines test_engines.cc)
    target_link_libraries(test_engines ${Boost_LIBRARIES})
endif ()
