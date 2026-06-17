// Copyright 2025-2026 Yanbin Li
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Yanbin Li

#ifndef HEONGPU_FFT_CORE_H
#define HEONGPU_FFT_CORE_H

#include "gpufft/fft.cuh"

#include <unordered_map>
#include <vector>

namespace heongpu
{
namespace fft
{
    template <typename T>
    __host__ void GPU_Special_FFT(COMPLEX<T>* device_inout,
                                  COMPLEX<T>* root_of_unity_table,
                                  gpufft::fft_configuration<T> cfg,
                                  int batch_size);

    template <typename T> auto CreateForwardSpecialFFTKernel()
    {
        return std::unordered_map<int, std::vector<gpufft::KernelConfig>>{
            {1,
             {{1, 1, 1, 1, 512 * sizeof(COMPLEX<T>), 0, 0, 0, 1, false}}},
            {2,
             {{1, 1, 2, 1, 512 * sizeof(COMPLEX<T>), 1, 1, 0, 2, false}}},
            {3,
             {{1, 1, 4, 1, 512 * sizeof(COMPLEX<T>), 2, 2, 0, 3, false}}},
            {4,
             {{1, 1, 8, 1, 512 * sizeof(COMPLEX<T>), 3, 3, 0, 4, false}}},
            {5,
             {{1, 1, 16, 1, 512 * sizeof(COMPLEX<T>), 4, 4, 0, 5, false}}},
            {6,
             {{1, 1, 32, 1, 512 * sizeof(COMPLEX<T>), 5, 5, 0, 6, false}}},
            {7,
             {{1, 1, 64, 1, 512 * sizeof(COMPLEX<T>), 6, 6, 0, 7, false}}},
            {8,
             {{1, 1, 128, 1, 512 * sizeof(COMPLEX<T>), 7, 7, 0, 8, false}}},
            {9,
             {{1, 1, 256, 1, 512 * sizeof(COMPLEX<T>), 8, 8, 0, 9, false}}},
            {10,
             {{1, 2, 256, 1, 512 * sizeof(COMPLEX<T>), 8, 9, 1, 9, false},
              {2, 1, 256, 1, 512 * sizeof(COMPLEX<T>), 8, 0, 0, 1, true}}}};
    }

    template <typename T> auto CreateInverseSpecialFFTKernel()
    {
        return std::unordered_map<int, std::vector<gpufft::KernelConfig>>{
            {1,
             {{1, 1, 1, 1, 512 * sizeof(COMPLEX<T>), 0, 0, 0, 1, false}}},
            {2,
             {{1, 1, 2, 1, 512 * sizeof(COMPLEX<T>), 1, 0, 0, 2, false}}},
            {3,
             {{1, 1, 4, 1, 512 * sizeof(COMPLEX<T>), 2, 0, 0, 3, false}}},
            {4,
             {{1, 1, 8, 1, 512 * sizeof(COMPLEX<T>), 3, 0, 0, 4, false}}},
            {5,
             {{1, 1, 16, 1, 512 * sizeof(COMPLEX<T>), 4, 0, 0, 5, false}}},
            {6,
             {{1, 1, 32, 1, 512 * sizeof(COMPLEX<T>), 5, 0, 0, 6, false}}},
            {7,
             {{1, 1, 64, 1, 512 * sizeof(COMPLEX<T>), 6, 0, 0, 7, false}}},
            {8,
             {{1, 1, 128, 1, 512 * sizeof(COMPLEX<T>), 7, 0, 0, 8, false}}},
            {9,
             {{1, 1, 256, 1, 512 * sizeof(COMPLEX<T>), 8, 0, 0, 9, false}}},
            {10,
             {{2, 1, 256, 1, 512 * sizeof(COMPLEX<T>), 8, 0, 0, 1, true},
              {1, 2, 256, 1, 512 * sizeof(COMPLEX<T>), 8, 1, 0, 9, false}}}};
    }

} // namespace fft
} // namespace heongpu

#endif // HEONGPU_FFT_CORE_H
