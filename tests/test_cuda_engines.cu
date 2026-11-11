// Copyright (c) 2000-2026, Heiko Bauke
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <cuda_runtime.h>

#include <trng/lcg64.hpp>
#include <trng/lcg64_count_shift.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/mrg2.hpp>
#include <trng/mrg3.hpp>
#include <trng/mrg3s.hpp>
#include <trng/mrg4.hpp>
#include <trng/mrg5.hpp>
#include <trng/mrg5s.hpp>
#include <trng/mwc128x.hpp>
#include <trng/mwc128x32.hpp>
#include <trng/mwc192x.hpp>
#include <trng/mwc64x32.hpp>
#include <trng/xoshiro256plus.hpp>
#include <trng/yarn2.hpp>
#include <trng/yarn3.hpp>
#include <trng/yarn3s.hpp>
#include <trng/yarn4.hpp>
#include <trng/yarn5.hpp>
#include <trng/yarn5s.hpp>


namespace {

  constexpr unsigned int thread_count{32};
  constexpr unsigned int sequence_length{256};


  template<typename Engine>
  __global__ void generate_sequences(Engine *engines, typename Engine::result_type *values) {
    const unsigned int thread{threadIdx.x};
    if (thread < thread_count) {
      for (unsigned int i{0}; i < sequence_length; ++i)
        values[thread * sequence_length + i] = engines[thread]();
    }
  }


  void check_cuda(const cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
      throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
  }


  template<typename Engine>
  void test_engine(const char *name) {
    static_assert(std::is_trivially_copyable<Engine>::value,
                  "CUDA tests require trivially copyable engine state");
    using result_type = typename Engine::result_type;

    std::array<Engine, thread_count> host_engines;
    for (auto &engine : host_engines)
      engine.seed(0ul);

    Engine engine;
    engine.seed(0ul);
    std::array<result_type, sequence_length> expected{};
    for (auto &value : expected)
      value = engine();

    Engine *device_engines{};
    result_type *device_values{};
    check_cuda(cudaMalloc(reinterpret_cast<void **>(&device_engines), sizeof(host_engines)),
               "cudaMalloc engine states");
    check_cuda(cudaMalloc(reinterpret_cast<void **>(&device_values),
                          thread_count * sequence_length * sizeof(result_type)),
               "cudaMalloc generated values");
    check_cuda(cudaMemcpy(device_engines, host_engines.data(), sizeof(host_engines),
                          cudaMemcpyHostToDevice),
               "cudaMemcpy engine states");

    generate_sequences<<<1, thread_count>>>(device_engines, device_values);
    check_cuda(cudaGetLastError(), "launch engine generation kernel");
    check_cuda(cudaDeviceSynchronize(), "synchronize engine generation kernel");

    std::array<result_type, thread_count * sequence_length> actual{};
    check_cuda(cudaMemcpy(actual.data(), device_values, sizeof(actual), cudaMemcpyDeviceToHost),
               "cudaMemcpy generated values");
    check_cuda(cudaFree(device_values), "cudaFree generated values");
    check_cuda(cudaFree(device_engines), "cudaFree engine states");

    for (unsigned int thread{0}; thread < thread_count; ++thread) {
      for (unsigned int i{0}; i < sequence_length; ++i) {
        if (actual[thread * sequence_length + i] != expected[i])
          throw std::runtime_error(std::string(name) + " differs at thread " +
                                   std::to_string(thread) + ", output " + std::to_string(i));
      }
    }
    std::cout << name << ": passed\n";
  }

}  // namespace


int main() {
  int device_count{0};
  const cudaError_t device_status{cudaGetDeviceCount(&device_count)};
  if (device_status == cudaErrorNoDevice or
      (device_status == cudaSuccess and device_count == 0)) {
    std::cerr << "CUDA engine tests skipped: no CUDA device is available\n";
    return 77;
  }

  try {
    check_cuda(device_status, "cudaGetDeviceCount");
    test_engine<trng::lcg64>("lcg64");
    test_engine<trng::lcg64_shift>("lcg64_shift");
    test_engine<trng::lcg64_count_shift>("lcg64_count_shift");
    test_engine<trng::mrg2>("mrg2");
    test_engine<trng::mrg3>("mrg3");
    test_engine<trng::mrg3s>("mrg3s");
    test_engine<trng::mrg4>("mrg4");
    test_engine<trng::mrg5>("mrg5");
    test_engine<trng::mrg5s>("mrg5s");
    test_engine<trng::yarn2>("yarn2");
    test_engine<trng::yarn3>("yarn3");
    test_engine<trng::yarn3s>("yarn3s");
    test_engine<trng::yarn4>("yarn4");
    test_engine<trng::yarn5>("yarn5");
    test_engine<trng::yarn5s>("yarn5s");
    test_engine<trng::xoshiro256plus>("xoshiro256plus");
    test_engine<trng::mwc64x32>("mwc64x32");
    test_engine<trng::mwc128x32>("mwc128x32");
    test_engine<trng::mwc128x>("mwc128x");
    test_engine<trng::mwc192x>("mwc192x");
  } catch (const std::exception &error) {
    std::cerr << "CUDA engine test failed: " << error.what() << '\n';
    return 1;
  }
}

// Excluded because their headers do not mark operator() as device-callable:
// lagfib2xor_521_64, lagfib4xor_521_32, lagfib2plus_521_32,
// lagfib4plus_521_64, mt19937, mt19937_64, and minstd.
// count128_lcg_shift is also excluded: its uint128 operations are host-only
// and produced a sequence that differed from the CPU reference.
