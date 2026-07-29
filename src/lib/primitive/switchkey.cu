// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#include <heongpu/kernel/switchkey.cuh>
#include <heongpu/primitive/ntt.cuh>
#include <heongpu/primitive/switchkey.cuh>
#include <heongpu/switchkey/switchkey.cuh>

namespace heongpu
{
namespace primitive
{
    void base_conversion_DtoQtilde_relin_leveled_ntt(
        Data64* ciphertext_coeff,
        Data64* ciphertext_ntt,
        Data64* output,
        Root64* roots,
        Modulus64* modulus,
        gpuntt::ntt_rns_configuration<Data64> cfg_ntt,
        Data64* base_change_matrix_D_to_Qtilda,
        Data64* Mi_inv_D_to_Qtilda,
        Data64* prod_D_to_Qtilda,
        int* I_j,
        int* I_location,
        int n_power,
        int d,
        int current_Qtilda_size,
        int current_Q_size,
        int first_Q_size,
        int level,
        int* mod_index,
        int* order,
        bool copy_excluded,
        const std::shared_ptr<PhantomNttTables>& tables,
        cudaStream_t stream)
    {
        if (!tables)
        {
            base_conversion_DtoQtilde_relin_leveled_kernel<<<
                dim3((1 << n_power) >> 8, d, 1), 256, 0, stream>>>(
                ciphertext_coeff, output, modulus, base_change_matrix_D_to_Qtilda,
                Mi_inv_D_to_Qtilda, prod_D_to_Qtilda, I_j, I_location, n_power,
                d, current_Qtilda_size, current_Q_size, level, mod_index);
            NTT_modulus_ordered_inplace(
                output, roots, modulus, cfg_ntt, d * current_Qtilda_size,
                current_Qtilda_size, order, tables);
            return;
        }

        switchkey::base_conversion_DtoQtilde_relin_leveled_modmajor(
            ciphertext_coeff, ciphertext_ntt, output, modulus,
            base_change_matrix_D_to_Qtilda, Mi_inv_D_to_Qtilda,
            prod_D_to_Qtilda, I_j, I_location, n_power, d,
            current_Qtilda_size, current_Q_size, level, copy_excluded, stream);

        NTT_modmajor_inplace(
            output, cfg_ntt, current_Qtilda_size, current_Q_size, first_Q_size,
            d, I_j, I_location, tables, copy_excluded);
    }

    void keyswitch_multiply_accumulate_leveled_method_II(
        Data64* input,
        const Data64* relinkey,
        Data64* output,
        Modulus64* modulus,
        int first_rns_mod_count,
        int current_decomp_mod_count,
        int current_rns_mod_count,
        int iteration_count1,
        int iteration_count2,
        int level,
        int n_power,
        const std::shared_ptr<PhantomNttTables>& tables,
        cudaStream_t stream)
    {
        if (!tables)
        {
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((1 << n_power) >> 8, current_rns_mod_count, 1), 256, 0,
                stream>>>(
                input, relinkey, output, modulus, first_rns_mod_count,
                current_decomp_mod_count, current_rns_mod_count,
                iteration_count1, iteration_count2, level, n_power);
            return;
        }

        switchkey::keyswitch_multiply_accumulate_leveled_method_II_modmajor(
            input, relinkey, output, modulus, first_rns_mod_count,
            current_decomp_mod_count, current_rns_mod_count, iteration_count1,
            iteration_count2, level, n_power, stream);
    }

} // namespace primitive
} // namespace heongpu
