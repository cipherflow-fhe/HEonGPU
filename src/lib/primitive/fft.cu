// Copyright 2025-2026 Yanbin Li
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Yanbin Li

#include <heongpu/primitive/fft.cuh>

#include <heongpu/fft/fft.cuh>

namespace heongpu
{
namespace primitive
{
    __host__ void special_fft(Complex64* device_inout,
                              Complex64* root_of_unity_table,
                              gpufft::fft_configuration<Float64> config,
                              int batch_size)
    {
        if (config.n_power <= 10)
        {
            fft::GPU_Special_FFT(device_inout, root_of_unity_table, config,
                                 batch_size);
            return;
        }

        gpufft::GPU_Special_FFT(device_inout, root_of_unity_table, config,
                                batch_size);
    }
} // namespace primitive
} // namespace heongpu
