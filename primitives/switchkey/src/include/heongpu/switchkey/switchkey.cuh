// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#ifndef HEONGPU_SWITCHKEY_PRIMITIVE_H
#define HEONGPU_SWITCHKEY_PRIMITIVE_H

#include <cuda_runtime.h>

#include <heongpu/util/util.cuh>

namespace heongpu
{
namespace switchkey
{
    void base_conversion_DtoQtilde_relin_leveled_modmajor(
        const Data64* ciphertext_coeff,
        const Data64* ciphertext_ntt,
        Data64* output,
        const Modulus64* modulus,
        const Data64* base_change_matrix_D_to_Qtilda,
        const Data64* Mi_inv_D_to_Qtilda,
        const Data64* prod_D_to_Qtilda,
        const int* I_j,
        const int* I_location,
        int n_power,
        int d,
        int current_Qtilda_size,
        int current_Q_size,
        int level,
        bool copy_excluded,
        cudaStream_t stream);

    void keyswitch_multiply_accumulate_leveled_method_II_modmajor(
        const Data64* input,
        const Data64* relinkey,
        Data64* output,
        const Modulus64* modulus,
        int first_rns_mod_count,
        int current_decomp_mod_count,
        int current_rns_mod_count,
        int iteration_count1,
        int iteration_count2,
        int level,
        int n_power,
        cudaStream_t stream);

} // namespace switchkey
} // namespace heongpu

#endif
