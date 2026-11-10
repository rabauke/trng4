add_executable(test_all test_engines.cc test_distributions.cc test_special_functions.cc test_int_math.cc test_linear_algebra.cc test_uint128.cc)
target_link_libraries(test_all PRIVATE trng4::trng4 Catch2::Catch2WithMain)
if (WIN32 AND BUILD_SHARED_LIBS)
  add_custom_command(TARGET test_all POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
      $<TARGET_RUNTIME_DLLS:test_all>
      $<TARGET_FILE_DIR:test_all>
      COMMAND_EXPAND_LISTS)
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/external/Catch2/contrib ${CMAKE_MODULE_PATH})
include(CTest)
include(Catch)

catch_discover_tests(test_all)

if (CMAKE_CUDA_COMPILER)
  add_executable(test_cuda_engines test_cuda_engines.cu)
  target_link_libraries(test_cuda_engines PRIVATE trng4::trng4)
  set_target_properties(test_cuda_engines PROPERTIES
      CUDA_STANDARD 17
      CUDA_STANDARD_REQUIRED YES)
  add_test(NAME cuda_engines COMMAND test_cuda_engines)
  set_tests_properties(cuda_engines PROPERTIES SKIP_RETURN_CODE 77)
endif ()
