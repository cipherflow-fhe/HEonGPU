// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#include <heongpu/switchkey/switchkey.cuh>

#include "gpuntt/common/modular_arith.cuh"

namespace heongpu
{
namespace
{
    template <bool CopyExcluded>
    __global__ void base_conversion_DtoQtilde_relin_leveled_modmajor_kernel(
        const Data64* __restrict__ ciphertext_coeff,
        const Data64* __restrict__ ciphertext_ntt,
        Data64* __restrict__ output,
        const Modulus64* __restrict__ modulus,
        const Data64* __restrict__ base_change_matrix_D_to_Qtilda,
        const Data64* __restrict__ Mi_inv_D_to_Qtilda,
        const Data64* __restrict__ prod_D_to_Qtilda,
        const int* __restrict__ I_j_,
        const int* __restrict__ I_location_, int n_power, int d,
        int current_Qtilda_size, int current_Q_size, int level)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int block_y = blockIdx.y;

        const int I_j = I_j_[block_y];
        int I_location = I_location_[block_y];

        int location = idx + (I_location << n_power);
        int matrix_index = I_location * current_Qtilda_size;

        Data64 partial[20];
        float r = 0.0F;
#pragma unroll
        for (int i = 0; i < I_j; i++)
        {
            Data64 temp = ciphertext_coeff[location + (i << n_power)];
            partial[i] =
                OPERATOR_GPU_64::mult(temp, Mi_inv_D_to_Qtilda[I_location + i],
                                      modulus[I_location + i]);
            r += static_cast<float>(partial[i]) /
                 static_cast<float>(modulus[I_location + i].value);
        }

        Data64 r_ = static_cast<Data64>(round(r));

        for (int i = 0; i < current_Qtilda_size; i++)
        {
            const int location_out =
                idx + (((i * d) + block_y) << n_power);

            if constexpr (CopyExcluded)
            {
                if (i >= I_location && i < I_location + I_j)
                {
                    output[location_out] =
                        ciphertext_ntt[idx + (i << n_power)];
                    continue;
                }
            }

            int mod_location = (i < current_Q_size) ? i : (i + level);

            Data64 temp = 0;
#pragma unroll
            for (int j = 0; j < I_j; j++)
            {
                Data64 mult = OPERATOR_GPU_64::reduce_forced(
                    partial[j], modulus[mod_location]);
                mult = OPERATOR_GPU_64::mult(
                    mult,
                    base_change_matrix_D_to_Qtilda[j + (i * I_j) +
                                                   matrix_index],
                    modulus[mod_location]);
                temp = OPERATOR_GPU_64::add(temp, mult, modulus[mod_location]);
            }

            Data64 r_mul = OPERATOR_GPU_64::mult(
                r_, prod_D_to_Qtilda[i + (block_y * current_Qtilda_size)],
                modulus[mod_location]);
            r_mul = OPERATOR_GPU_64::sub(temp, r_mul, modulus[mod_location]);
            output[location_out] = r_mul;
        }
    }

    __global__ void
    keyswitch_multiply_accumulate_leveled_method_II_modmajor_kernel(
        const Data64* __restrict__ input,
        const Data64* __restrict__ relinkey, Data64* __restrict__ output,
        const Modulus64* __restrict__ modulus, int first_rns_mod_count,
        int current_decomp_mod_count, int current_rns_mod_count,
        int iteration_count1, int iteration_count2, int level, int n_power)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int block_y = blockIdx.y;
        int key_index =
            (block_y < current_decomp_mod_count) ? block_y : (block_y + level);

        int key_offset1 = first_rns_mod_count << n_power;
        int key_offset2 = first_rns_mod_count << (n_power + 1);
        int index = idx + (block_y << n_power);
        int index2 = idx + (key_index << n_power);

        Modulus64 modulus_reg = modulus[key_index];
        Data64 ct_0_sum = 0;
        Data64 ct_1_sum = 0;

#pragma unroll
        for (int i = 0; i < iteration_count1; i++)
        {
            const int group0 = 4 * i;
            const std::size_t base =
                (static_cast<std::size_t>(block_y) *
                     (iteration_count1 * 4 + iteration_count2) +
                 group0)
                << n_power;

            Data64 in_piece1 = input[base + idx];
            Data64 rk0_1 =
                __ldg(&relinkey[index2 + key_offset2 * group0]);
            Data64 rk1_1 =
                __ldg(&relinkey[index2 + key_offset2 * group0 + key_offset1]);

            Data64 in_piece2 = input[base + (1 << n_power) + idx];
            Data64 rk0_2 =
                __ldg(&relinkey[index2 + key_offset2 * (group0 + 1)]);
            Data64 rk1_2 = __ldg(
                &relinkey[index2 + key_offset2 * (group0 + 1) + key_offset1]);

            Data64 in_piece3 = input[base + (2 << n_power) + idx];
            Data64 rk0_3 =
                __ldg(&relinkey[index2 + key_offset2 * (group0 + 2)]);
            Data64 rk1_3 = __ldg(
                &relinkey[index2 + key_offset2 * (group0 + 2) + key_offset1]);

            Data64 in_piece4 = input[base + (3 << n_power) + idx];
            Data64 rk0_4 =
                __ldg(&relinkey[index2 + key_offset2 * (group0 + 3)]);
            Data64 rk1_4 = __ldg(
                &relinkey[index2 + key_offset2 * (group0 + 3) + key_offset1]);

            ct_0_sum = OPERATOR_GPU_64::add(
                ct_0_sum, OPERATOR_GPU_64::mult(in_piece1, rk0_1, modulus_reg),
                modulus_reg);
            ct_1_sum = OPERATOR_GPU_64::add(
                ct_1_sum, OPERATOR_GPU_64::mult(in_piece1, rk1_1, modulus_reg),
                modulus_reg);
            ct_0_sum = OPERATOR_GPU_64::add(
                ct_0_sum, OPERATOR_GPU_64::mult(in_piece2, rk0_2, modulus_reg),
                modulus_reg);
            ct_1_sum = OPERATOR_GPU_64::add(
                ct_1_sum, OPERATOR_GPU_64::mult(in_piece2, rk1_2, modulus_reg),
                modulus_reg);
            ct_0_sum = OPERATOR_GPU_64::add(
                ct_0_sum, OPERATOR_GPU_64::mult(in_piece3, rk0_3, modulus_reg),
                modulus_reg);
            ct_1_sum = OPERATOR_GPU_64::add(
                ct_1_sum, OPERATOR_GPU_64::mult(in_piece3, rk1_3, modulus_reg),
                modulus_reg);
            ct_0_sum = OPERATOR_GPU_64::add(
                ct_0_sum, OPERATOR_GPU_64::mult(in_piece4, rk0_4, modulus_reg),
                modulus_reg);
            ct_1_sum = OPERATOR_GPU_64::add(
                ct_1_sum, OPERATOR_GPU_64::mult(in_piece4, rk1_4, modulus_reg),
                modulus_reg);
        }

        int loop_offset = iteration_count1 * 4;
#pragma unroll
        for (int i = loop_offset; i < loop_offset + iteration_count2; i++)
        {
            const std::size_t input_index =
                ((static_cast<std::size_t>(block_y) *
                      (iteration_count1 * 4 + iteration_count2) +
                  i)
                 << n_power) +
                idx;
            Data64 in_piece1 = input[input_index];
            Data64 rk0_1 = __ldg(&relinkey[index2 + (key_offset2 * i)]);
            Data64 rk1_1 =
                __ldg(&relinkey[index2 + (key_offset2 * i) + key_offset1]);

            ct_0_sum = OPERATOR_GPU_64::add(
                ct_0_sum, OPERATOR_GPU_64::mult(in_piece1, rk0_1, modulus_reg),
                modulus_reg);
            ct_1_sum = OPERATOR_GPU_64::add(
                ct_1_sum, OPERATOR_GPU_64::mult(in_piece1, rk1_1, modulus_reg),
                modulus_reg);
        }

        output[index] = ct_0_sum;
        output[index + (current_rns_mod_count << n_power)] = ct_1_sum;
    }

    __global__ void bs_add_permute_fused_kernel(
        const Data64* input,
        const Data64* addend,
        Data64* output,
        const Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int block_y = blockIdx.y;
        int block_z = blockIdx.z;

        int shift = 32 - n_power;
        int two_N = 2 << n_power;
        int br_j = __brev(idx) >> shift;
        int exp_j = 2 * br_j + 1;
        int new_exp =
            static_cast<int>((static_cast<long long>(galois_elt) * exp_j) %
                             two_N);
        int src_idx = __brev((new_exp - 1) >> 1) >> shift;

        int src_offset = src_idx + (block_y << n_power) +
                         ((pql_count << n_power) * block_z);
        int dst_offset =
            idx + (block_y << n_power) + ((pql_count << n_power) * block_z);

        Data64 value = input[src_offset];
        if (block_z == 0)
        {
            value = OPERATOR_GPU_64::add(
                value, addend[src_idx + (block_y << n_power)],
                pq_modulus[block_y]);
        }

        output[dst_offset] = value;
    }

    __global__ void gs_add_permute_acc_fused_kernel(
        const Data64* input,
        const Data64* addend,
        const Data64* accum,
        Data64* output,
        const Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int block_y = blockIdx.y;
        int block_z = blockIdx.z;

        int shift = 32 - n_power;
        int two_N = 2 << n_power;
        int br_j = __brev(idx) >> shift;
        int exp_j = 2 * br_j + 1;
        int new_exp =
            static_cast<int>((static_cast<long long>(galois_elt) * exp_j) %
                             two_N);
        int src_idx = __brev((new_exp - 1) >> 1) >> shift;

        int src_offset = src_idx + (block_y << n_power) +
                         ((pql_count << n_power) * block_z);
        int dst_offset =
            idx + (block_y << n_power) + ((pql_count << n_power) * block_z);

        Data64 value = input[src_offset];
        if (block_z == 0)
        {
            value = OPERATOR_GPU_64::add(
                value, addend[src_idx + (block_y << n_power)],
                pq_modulus[block_y]);
        }

        output[dst_offset] =
            OPERATOR_GPU_64::add(accum[dst_offset], value, pq_modulus[block_y]);
    }

    template <int PSize>
    __global__ void divide_round_lastq_extended_leveled_add_first_kernel(
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
        int p_size)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int block_y = blockIdx.y;
        int block_z = blockIdx.z;

        Data64 last_ct[15];
        const int active_p_size = PSize == 0 ? p_size : PSize;
#pragma unroll
        for (int i = 0; i < active_p_size; ++i)
        {
            last_ct[i] =
                input[idx + ((q_size + i) << n_power) +
                      ((q_prime_size << n_power) * block_z)];
        }

        Data64 input_ = input[idx + (block_y << n_power) +
                              ((q_prime_size << n_power) * block_z)];

        int location = 0;
#pragma unroll
        for (int i = 0; i < active_p_size; ++i)
        {
            Data64 last_ct_add_half = last_ct[active_p_size - 1 - i];
            last_ct_add_half =
                OPERATOR_GPU_64::add(last_ct_add_half, half[i],
                                     modulus[first_q_prime_size - 1 - i]);
            for (int j = 0; j < (active_p_size - 1 - i); ++j)
            {
                Data64 temp = OPERATOR_GPU_64::reduce_forced(
                    last_ct_add_half, modulus[first_q_size + j]);
                temp = OPERATOR_GPU_64::sub(
                    temp, half_mod[location + first_q_size + j],
                    modulus[first_q_size + j]);
                temp = OPERATOR_GPU_64::sub(last_ct[j], temp,
                                            modulus[first_q_size + j]);
                last_ct[j] = OPERATOR_GPU_64::mult(
                    temp, last_q_modinv[location + first_q_size + j],
                    modulus[first_q_size + j]);
            }

            Data64 temp = OPERATOR_GPU_64::reduce_forced(
                last_ct_add_half, modulus[block_y]);
            temp = OPERATOR_GPU_64::sub(temp, half_mod[location + block_y],
                                        modulus[block_y]);
            temp = OPERATOR_GPU_64::sub(input_, temp, modulus[block_y]);
            input_ = OPERATOR_GPU_64::mult(
                temp, last_q_modinv[location + block_y], modulus[block_y]);

            location += first_q_prime_size - 1 - i;
        }

        if (block_z == 0)
        {
            input_ = OPERATOR_GPU_64::add(
                input_, addend_first[idx + (block_y << n_power)],
                modulus[block_y]);
        }

        output[idx + (block_y << n_power) + ((q_size << n_power) * block_z)] =
            input_;
    }

} // namespace

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
        cudaStream_t stream)
    {
        if (copy_excluded)
        {
            base_conversion_DtoQtilde_relin_leveled_modmajor_kernel<true><<<
                dim3((1 << n_power) >> 8, d, 1), 256, 0, stream>>>(
                ciphertext_coeff, ciphertext_ntt, output, modulus,
                base_change_matrix_D_to_Qtilda, Mi_inv_D_to_Qtilda,
                prod_D_to_Qtilda, I_j, I_location, n_power, d,
                current_Qtilda_size, current_Q_size, level);
            return;
        }

        base_conversion_DtoQtilde_relin_leveled_modmajor_kernel<false><<<
            dim3((1 << n_power) >> 8, d, 1), 256, 0, stream>>>(
            ciphertext_coeff, nullptr, output, modulus,
            base_change_matrix_D_to_Qtilda, Mi_inv_D_to_Qtilda,
            prod_D_to_Qtilda, I_j, I_location, n_power, d,
            current_Qtilda_size, current_Q_size, level);
    }

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
        cudaStream_t stream)
    {
        keyswitch_multiply_accumulate_leveled_method_II_modmajor_kernel<<<
            dim3((1 << n_power) >> 8, current_rns_mod_count, 1), 256, 0,
            stream>>>(
            input, relinkey, output, modulus, first_rns_mod_count,
            current_decomp_mod_count, current_rns_mod_count, iteration_count1,
            iteration_count2, level, n_power);
    }

    void bs_add_permute_fused(
        const Data64* input,
        const Data64* addend,
        Data64* output,
        const Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count,
        cudaStream_t stream)
    {
        bs_add_permute_fused_kernel<<<
            dim3((1 << n_power) >> 8, pql_count, 2), 256, 0, stream>>>(
            input, addend, output, pq_modulus, galois_elt, n_power, pql_count);
    }

    void gs_add_permute_acc_fused(
        const Data64* input,
        const Data64* addend,
        const Data64* accum,
        Data64* output,
        const Modulus64* pq_modulus,
        int galois_elt,
        int n_power,
        int pql_count,
        cudaStream_t stream)
    {
        gs_add_permute_acc_fused_kernel<<<
            dim3((1 << n_power) >> 8, pql_count, 2), 256, 0, stream>>>(
            input, addend, accum, output, pq_modulus, galois_elt, n_power,
            pql_count);
    }

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
        cudaStream_t stream)
    {
        const dim3 grid((1 << n_power) >> 8, q_size, 2);
        switch (p_size)
        {
            case 1:
                divide_round_lastq_extended_leveled_add_first_kernel<1>
                    <<<grid, 256, 0, stream>>>(
                        input, addend_first, output, modulus, half, half_mod,
                        last_q_modinv, n_power, q_prime_size, q_size,
                        first_q_prime_size, first_q_size, p_size);
                return;
            case 2:
                divide_round_lastq_extended_leveled_add_first_kernel<2>
                    <<<grid, 256, 0, stream>>>(
                        input, addend_first, output, modulus, half, half_mod,
                        last_q_modinv, n_power, q_prime_size, q_size,
                        first_q_prime_size, first_q_size, p_size);
                return;
            case 3:
                divide_round_lastq_extended_leveled_add_first_kernel<3>
                    <<<grid, 256, 0, stream>>>(
                        input, addend_first, output, modulus, half, half_mod,
                        last_q_modinv, n_power, q_prime_size, q_size,
                        first_q_prime_size, first_q_size, p_size);
                return;
            case 4:
                divide_round_lastq_extended_leveled_add_first_kernel<4>
                    <<<grid, 256, 0, stream>>>(
                        input, addend_first, output, modulus, half, half_mod,
                        last_q_modinv, n_power, q_prime_size, q_size,
                        first_q_prime_size, first_q_size, p_size);
                return;
            case 5:
                divide_round_lastq_extended_leveled_add_first_kernel<5>
                    <<<grid, 256, 0, stream>>>(
                        input, addend_first, output, modulus, half, half_mod,
                        last_q_modinv, n_power, q_prime_size, q_size,
                        first_q_prime_size, first_q_size, p_size);
                return;
            default:
                divide_round_lastq_extended_leveled_add_first_kernel<0>
                    <<<grid, 256, 0, stream>>>(
                        input, addend_first, output, modulus, half, half_mod,
                        last_q_modinv, n_power, q_prime_size, q_size,
                        first_q_prime_size, first_q_size, p_size);
                return;
        }
    }

} // namespace switchkey
} // namespace heongpu
