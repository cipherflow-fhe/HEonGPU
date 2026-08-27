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

    void bs_add_permute_fused(
        const Data64* input,
        const Data64* addend,
        Data64* output,
        const Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count,
        cudaStream_t stream);

    void gs_add_permute_acc_fused(
        const Data64* input,
        const Data64* addend,
        const Data64* accum,
        Data64* output,
        const Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count,
        cudaStream_t stream);

    void divide_round_lastq_extended_leveled_add_first(
        const Data64* input,
        const Data64* addend_first,
        Data64* output,
        const Modulus64* modulus,
        const Data64* half,
        const Data64* half_mod,
        const Data64* last_q_modinv,
        int n_power,
        int q_prime_size,
        int q_size,
        int first_q_prime_size,
        int first_q_size,
        int p_size,
        cudaStream_t stream);

} // namespace switchkey
} // namespace heongpu

#endif
