// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#include <heongpu/kernel/switchkey.cuh>
#include <heongpu/primitive/ntt.cuh>
#include <heongpu/primitive/switchkey.cuh>
#include <heongpu/switchkey/switchkey.cuh>
#include <heongpu/util/util.cuh>

#include <cstdlib>

namespace heongpu
{
namespace primitive
{
namespace
{
    bool use_modmajor_keyswitch(const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            return false;
        }

        const char* enabled = std::getenv("HEONGPU_USE_MOD_KSWITCH");
        return !enabled || enabled[0] != '0';
    }

    bool use_bsgs_fusion(const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            return false;
        }

        const char* enabled = std::getenv("HEONGPU_USE_BSGS_FUSION");
        return !enabled || enabled[0] != '0';
    }

    bool use_keyswitch_part2(const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            return false;
        }

        const char* enabled = std::getenv("HEONGPU_USE_KSWITCH_P2");
        return !enabled || enabled[0] != '0';
    }

} // namespace

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
        if (!use_modmajor_keyswitch(tables))
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
        if (!use_modmajor_keyswitch(tables))
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

    void bs_add_permute_fused(
        Data64* input,
        Data64* addend,
        Data64* output,
        Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count,
        const std::shared_ptr<PhantomNttTables>& tables,
        cudaStream_t stream)
    {
        if (!use_bsgs_fusion(tables))
        {
            addition_pql_kernel<<<dim3((1 << n_power) >> 8, pql_count, 1),
                                  256, 0, stream>>>(
                input, addend, input, pq_modulus, n_power, pql_count);
            galois_permute_ntt_pql_kernel<<<
                dim3((1 << n_power) >> 8, pql_count, 2), 256, 0, stream>>>(
                input, output, galois_elt, n_power, pql_count);
            return;
        }

        switchkey::bs_add_permute_fused(
            input, addend, output, pq_modulus, galois_elt, n_power, pql_count,
            stream);
    }

    void gs_add_permute_acc_fused(
        Data64* input,
        Data64* addend,
        Data64* accum,
        Data64* scratch,
        Data64* output,
        Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count,
        const std::shared_ptr<PhantomNttTables>& tables,
        cudaStream_t stream)
    {
        if (!use_bsgs_fusion(tables))
        {
            addition_pql_kernel<<<dim3((1 << n_power) >> 8, pql_count, 1),
                                  256, 0, stream>>>(
                input, addend, input, pq_modulus, n_power, pql_count);
            galois_permute_ntt_pql_kernel<<<
                dim3((1 << n_power) >> 8, pql_count, 2), 256, 0, stream>>>(
                input, scratch, galois_elt, n_power, pql_count);
            addition_pql_kernel<<<dim3((1 << n_power) >> 8, pql_count, 2),
                                  256, 0, stream>>>(
                accum, scratch, output, pq_modulus, n_power, pql_count);
            return;
        }

        switchkey::gs_add_permute_acc_fused(
            input, addend, accum, output, pq_modulus, galois_elt, n_power,
            pql_count, stream);
    }

    void keyswitch_part2_fused_moddown_ntt(
        Data64* input,
        Data64* addend_first,
        Data64* scratch,
        Data64* output,
        Root64* roots,
        Modulus64* modulus,
        gpuntt::ntt_rns_configuration<Data64> cfg_ntt,
        Data64* half,
        Data64* half_mod,
        Data64* last_q_modinv,
        int n_power,
        int q_prime_size,
        int q_size,
        int first_q_prime_size,
        int first_q_size,
        int p_size,
        const std::shared_ptr<PhantomNttTables>& tables,
        cudaStream_t stream)
    {
        if (!use_keyswitch_part2(tables))
        {
            divide_round_lastq_extended_leveled_kernel<<<
                dim3((1 << n_power) >> 8, q_size, 2), 256, 0, stream>>>(
                input, output, modulus, half, half_mod, last_q_modinv,
                n_power, q_prime_size, q_size, first_q_prime_size,
                first_q_size, p_size);
            NTT_inplace(output, roots, modulus, cfg_ntt, 2 * q_size, q_size,
                        tables);

            NTT_inplace(addend_first, roots, modulus, cfg_ntt, q_size, q_size,
                        tables);
            addition_switchkey<<<dim3((1 << n_power) >> 8, q_size, 2), 256, 0,
                                 stream>>>(
                output, addend_first, output, modulus, n_power);
            return;
        }

        switchkey::divide_round_lastq_extended_leveled_add_first(
            input, addend_first, output, modulus, half, half_mod, last_q_modinv,
            n_power, q_prime_size, q_size, first_q_prime_size, first_q_size,
            p_size, stream);
        NTT_inplace(output, roots, modulus, cfg_ntt, 2 * q_size, q_size,
                    tables);
    }

} // namespace primitive
} // namespace heongpu
