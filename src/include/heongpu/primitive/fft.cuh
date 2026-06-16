// Copyright 2025-2026 Yanbin Li
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Yanbin Li

#ifndef HEONGPU_PRIMITIVE_FFT_H
#define HEONGPU_PRIMITIVE_FFT_H

#include "gpufft/fft.cuh"

namespace heongpu
{
namespace primitive
{
    __host__ void special_fft(Complex64* device_inout,
                              Complex64* root_of_unity_table,
                              gpufft::fft_configuration<Float64> config,
                              int batch_size);
} // namespace primitive
} // namespace heongpu

#endif // HEONGPU_PRIMITIVE_FFT_H
