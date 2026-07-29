// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#include <heongpu/ntt/ntt.cuh>

#include <vector>

#include "gpuntt/common/modular_arith.cuh"
#include "phantom_ntt_util.cuh"

namespace heongpu
{
namespace ntt
{
    struct PhantomNttTables
    {
        phantom_local::DNTTTable table;
    };

namespace
{
    using phantom_local::DNTTTable;
    using phantom_local::DModulus;

    using phantom_local::blockDimNTT;
    using phantom_local::gridDimNTT;
    using phantom_local::per_block_pad;
    using phantom_local::per_thread_sample_size;
    using phantom_local::sample_size;

    using phantom_local::csub_q;
    using phantom_local::ct_butterfly;
    using phantom_local::fntt4;
    using phantom_local::fntt8;
    using phantom_local::gs_butterfly;
    using phantom_local::intt4;
    using phantom_local::intt8;
    using phantom_local::multiply_and_reduce_shoup_lazy;

    enum class NTTMode
    {
        Inplace,
        ModulusOrdered,
        PolyOrdered
    };

    std::uint64_t compute_shoup(std::uint64_t operand, std::uint64_t modulus)
    {
        return static_cast<std::uint64_t>(
            (static_cast<__uint128_t>(operand) << 64) / modulus);
    }

    DModulus make_dmodulus(std::uint64_t modulus)
    {
        const __uint128_t quotient =
            (~static_cast<__uint128_t>(0) / modulus) + 1;
        return DModulus(
            modulus, static_cast<std::uint64_t>(quotient),
            static_cast<std::uint64_t>(quotient >> 64));
    }

    std::vector<std::uint64_t> make_shoup_table(
        const std::vector<std::uint64_t>& roots, std::uint64_t modulus)
    {
        std::vector<std::uint64_t> out(roots.size());
        for (std::size_t i = 0; i < roots.size(); ++i)
        {
            out[i] = compute_shoup(roots[i], modulus);
        }
        return out;
    }

    template <NTTMode Mode, bool Phase2>
    __device__ __forceinline__ void ntt_table_access(
        std::size_t logical_row_idx, std::size_t mod_count,
        std::size_t start_mod_idx, const int* order,
        std::size_t& data_row_idx, std::size_t& table_mod_idx)
    {
        const std::size_t local_row_idx = logical_row_idx % mod_count;
        const std::size_t data_row_base = logical_row_idx - local_row_idx;
        std::size_t local_mod_idx = local_row_idx;

        if constexpr (Phase2)
        {
            local_mod_idx = mod_count - 1 - local_row_idx;
        }

        data_row_idx = data_row_base + local_mod_idx;
        if constexpr (Mode == NTTMode::ModulusOrdered)
        {
            table_mod_idx = static_cast<std::size_t>(order[local_mod_idx]);
        }
        else
        {
            table_mod_idx = local_mod_idx + start_mod_idx;
        }
    }

    template <NTTMode Mode>
    __device__ __forceinline__ void intt_table_access(
        std::size_t logical_row_idx, std::size_t mod_count,
        std::size_t start_mod_idx, const int* order,
        std::size_t& data_row_idx, std::size_t& table_mod_idx)
    {
        const std::size_t local_mod_idx = logical_row_idx % mod_count;
        if constexpr (Mode == NTTMode::Inplace)
        {
            data_row_idx = logical_row_idx;
            table_mod_idx = local_mod_idx + start_mod_idx;
        }
        else if constexpr (Mode == NTTMode::ModulusOrdered)
        {
            data_row_idx = logical_row_idx;
            table_mod_idx = static_cast<std::size_t>(order[local_mod_idx]);
        }
        else
        {
            data_row_idx = static_cast<std::size_t>(order[logical_row_idx]);
            table_mod_idx = local_mod_idx + start_mod_idx;
        }
    }

    template <NTTMode Mode>
    __global__ void batched_inplace_fnwt_radix8_phase1(
        std::uint64_t* inout, const std::uint64_t* twiddles,
        const std::uint64_t* twiddles_shoup, const DModulus* modulus,
        std::size_t coeff_mod_size, std::size_t start_mod_idx,
        std::size_t poly_count, const int* order, std::size_t n,
        std::size_t n1, std::size_t pad)
    {
        extern __shared__ std::uint64_t buffer[];

        std::size_t pad_tid = threadIdx.x % pad;
        std::size_t pad_idx = threadIdx.x / pad;
        std::size_t group = n1 / 8;
        std::uint64_t samples[8];
        std::size_t t = n / 2;

        for (std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
             tid < (n / 8) * (coeff_mod_size * poly_count);
             tid += blockDim.x * gridDim.x)
        {
            std::size_t row_idx = tid / (n / 8);
            std::size_t data_row_idx = 0;
            std::size_t table_mod_idx = 0;
            ntt_table_access<Mode, false>(
                row_idx, coeff_mod_size, start_mod_idx, order, data_row_idx,
                table_mod_idx);
            std::size_t n_idx = tid % (n / 8);

            std::uint64_t* data_ptr = inout + data_row_idx * n;
            const std::uint64_t* psi = twiddles + table_mod_idx * n;
            const std::uint64_t* psi_shoup = twiddles_shoup + table_mod_idx * n;
            std::uint64_t modulus_value = modulus[table_mod_idx].value();
            std::size_t n_init = t / 4 / group * pad_idx + pad_tid + pad * (n_idx / (group * pad));

            for (std::size_t j = 0; j < 8; j++){
                samples[j] = *(data_ptr + n_init + t / 4 * j);
            }
            std::size_t tw_idx = 1;
            fntt8(samples, psi, psi_shoup, tw_idx, modulus_value);
            for (std::size_t j = 0; j < 8; j++){
                buffer[pad_tid * (n1 + pad) + pad_idx + group * j] = samples[j];
            }
            std::size_t remain_iters = 0;
            __syncthreads();
            for (std::size_t j = 8, k = group / 2; j < group + 1; j *= 8, k >>= 3){

                std::size_t m_idx2 = pad_idx / (k / 4);
                std::size_t t_idx2 = pad_idx % (k / 4);
                for (std::size_t l = 0; l < 8; l++){
                    samples[l] = buffer[(n1 + pad) * pad_tid + 2 * m_idx2 * k + t_idx2 + (k / 4) * l];
                }
                std::size_t tw_idx2 = j * tw_idx + m_idx2;
                fntt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
                for (std::size_t l = 0; l < 8; l++){
                    buffer[(n1 + pad) * pad_tid + 2 * m_idx2 * k + t_idx2 + (k / 4) * l] = samples[l];
                }
                if (j == group / 2) remain_iters = 1;
                if (j == group / 4) remain_iters = 2;
                __syncthreads();
            }

            if (group < 8) remain_iters = (group == 4) ? 2 : 1;
            for (std::size_t l = 0; l < 8; l++){
                samples[l] = buffer[(n1 + pad) * pad_tid + 8 * pad_idx + l];
            }
            if (remain_iters == 1){
                std::size_t tw_idx2 = 4 * group * tw_idx + 4 * pad_idx;
                ct_butterfly(samples[0], samples[1], psi[tw_idx2], psi_shoup[tw_idx2], modulus_value);
                ct_butterfly(samples[2], samples[3], psi[tw_idx2 + 1], psi_shoup[tw_idx2 + 1], modulus_value);
                ct_butterfly(samples[4], samples[5], psi[tw_idx2 + 2], psi_shoup[tw_idx2 + 2], modulus_value);
                ct_butterfly(samples[6], samples[7], psi[tw_idx2 + 3], psi_shoup[tw_idx2 + 3], modulus_value);
            }
            else if (remain_iters == 2){
                std::size_t tw_idx2 = 2 * group * tw_idx + 2 * pad_idx;
                fntt4(samples, psi, psi_shoup, tw_idx2, modulus_value);
                fntt4(samples + 4, psi, psi_shoup, tw_idx2 + 1, modulus_value);
            }
            for (std::size_t l = 0; l < 8; l++){
                buffer[(n1 + pad) * pad_tid + 8 * pad_idx + l] = samples[l];
            }

            __syncthreads();
            for (std::size_t j = 0; j < 8; j++){
                *(data_ptr + n_init + t / 4 * j) =
                    buffer[pad_tid * (n1 + pad) + pad_idx + group * j];
            }
        }
    }

    template <NTTMode Mode>
    __global__ void batched_inplace_fnwt_radix8_phase2(
        std::uint64_t* inout, const std::uint64_t* twiddles,
        const std::uint64_t* twiddles_shoup, const DModulus* modulus,
        std::size_t coeff_mod_size, std::size_t start_mod_idx,
        std::size_t poly_count, const int* order, std::size_t n,
        std::size_t n1,
        std::size_t n2)
    {
        extern __shared__ std::uint64_t buffer[];

        std::size_t group = n2 / 8;
        std::size_t set = threadIdx.x / group;
        std::uint64_t samples[8];
        std::size_t t = n2 / 2;

        for (std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
             tid < (n / 8) * (coeff_mod_size * poly_count);
             tid += blockDim.x * gridDim.x)
        {
            std::size_t row_idx = tid / (n / 8);
            std::size_t data_row_idx = 0;
            std::size_t table_mod_idx = 0;
            ntt_table_access<Mode, true>(
                row_idx, coeff_mod_size, start_mod_idx, order, data_row_idx,
                table_mod_idx);
            std::size_t n_idx = tid % (n / 8);
            std::size_t m_idx = n_idx / (t / 4);
            std::size_t t_idx = n_idx % (t / 4);

            std::uint64_t* data_ptr = inout + data_row_idx * n;
            std::uint64_t modulus_value = modulus[table_mod_idx].value();
            const std::uint64_t* psi = twiddles + n * table_mod_idx;
            const std::uint64_t* psi_shoup = twiddles_shoup + n * table_mod_idx;
            std::size_t n_init = 2 * m_idx * t + t_idx;
            for (std::size_t j = 0; j < 8; j++){
                samples[j] = *(data_ptr + n_init + t / 4 * j);
            }
            std::size_t tw_idx = n1 + m_idx;
            fntt8(samples, psi, psi_shoup, tw_idx, modulus_value);
            for (std::size_t j = 0; j < 8; j++){
                buffer[set * n2 + t_idx + t / 4 * j] = samples[j];
            }
            std::size_t tail = 0;
            __syncthreads();

            for (std::size_t j = 8, k = t / 8; j < t / 4 + 1; j *= 8, k >>= 3){
                std::size_t m_idx2 = t_idx / (k / 4);
                std::size_t t_idx2 = t_idx % (k / 4);
                for (std::size_t l = 0; l < 8; l++){
                    samples[l] = buffer[set * n2 + 2 * m_idx2 * k + t_idx2 + (k / 4) * l];
                }
                std::size_t tw_idx2 = j * tw_idx + m_idx2;
                fntt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
                for (std::size_t l = 0; l < 8; l++){
                    buffer[set * n2 + 2 * m_idx2 * k + t_idx2 + (k / 4) * l] = samples[l];
                }
                if (j == t / 8) tail = 1;
                if (j == t / 16) tail = 2;
                __syncthreads();
            }

            for (std::size_t l = 0; l < 8; l++){
                samples[l] = buffer[set * n2 + 8 * t_idx + l];
            }
            if (tail == 1){
                std::size_t tw_idx2 = t * tw_idx + 4 * t_idx;
                ct_butterfly(samples[0], samples[1], psi[tw_idx2], psi_shoup[tw_idx2], modulus_value);
                ct_butterfly(samples[2], samples[3], psi[tw_idx2 + 1], psi_shoup[tw_idx2 + 1], modulus_value);
                ct_butterfly(samples[4], samples[5], psi[tw_idx2 + 2], psi_shoup[tw_idx2 + 2], modulus_value);
                ct_butterfly(samples[6], samples[7], psi[tw_idx2 + 3], psi_shoup[tw_idx2 + 3], modulus_value);
            }
            else if (tail == 2){
                std::size_t tw_idx2 = (t / 2) * tw_idx + 2 * t_idx;
                fntt4(samples, psi, psi_shoup, tw_idx2, modulus_value);
                fntt4(samples + 4, psi, psi_shoup, tw_idx2 + 1, modulus_value);
            }
            for (std::size_t l = 0; l < 8; l++){
                buffer[set * n2 + 8 * t_idx + l] = samples[l];
            }
            __syncthreads();

            std::uint64_t modulus2 = modulus_value << 1;
            for (std::size_t j = 0; j < 8; j++){
                samples[j] = buffer[set * n2 + t_idx + t / 4 * j];
                csub_q(samples[j], modulus2);
                csub_q(samples[j], modulus_value);
            }
            for (std::size_t j = 0; j < 8; j++){
                *(data_ptr + n_init + t / 4 * j) = samples[j];
            }
        }
    }

    template <NTTMode Mode>
    __global__ void batched_inplace_inwt_radix8_phase1(
        const std::uint64_t* input, std::uint64_t* output,
        const std::uint64_t* itwiddles,
        const std::uint64_t* itwiddles_shoup, const DModulus* modulus,
        std::size_t coeff_mod_size, std::size_t start_mod_idx,
        std::size_t poly_count, const int* order, std::size_t n,
        std::size_t n1, std::size_t n2)
	    {
	        extern __shared__ std::uint64_t buffer[];

	        for (std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	             i < (n / 8) * (coeff_mod_size * poly_count);
	             i += blockDim.x * gridDim.x)
	        {
	            std::size_t group = n2 / 8;
	            std::size_t set = threadIdx.x / group;
	            std::uint64_t samples[8];
	            std::size_t t = n / 2 / n1;

	            std::size_t row_idx = i / (n / 8);
	            std::size_t data_row_idx = 0;
	            std::size_t twr_idx = 0;
                intt_table_access<Mode>(
                    row_idx, coeff_mod_size, start_mod_idx, order,
                    data_row_idx, twr_idx);
	            std::size_t n_idx = i % (n / 8);
	            std::size_t m_idx = n_idx / (t / 4);
	            std::size_t t_idx = n_idx % (t / 4);

	            const std::uint64_t* input_ptr = input + data_row_idx * n;
	            std::uint64_t* output_ptr = output + data_row_idx * n;
	            const std::uint64_t* psi = itwiddles + n * twr_idx;
	            const std::uint64_t* psi_shoup = itwiddles_shoup + n * twr_idx;
	            std::uint64_t modulus_value = modulus[twr_idx].value();
	            std::size_t n_init = 2 * m_idx * t + t_idx;

	            for (std::size_t j = 0; j < 8; j++){
		                buffer[set * n2 + t_idx + t / 4 * j] = *(input_ptr + n_init + t / 4 * j);
	            }
	            __syncthreads();

	            for (std::size_t l = 0; l < 8; l++){
	                samples[l] = buffer[set * n2 + 8 * t_idx + l];
	            }
	            std::size_t tw_idx = n1 + m_idx;
	            std::size_t tw_idx2 = (t / 4) * tw_idx + t_idx;
	            intt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
	            for (std::size_t l = 0; l < 8; l++){
	                buffer[set * n2 + 8 * t_idx + l] = samples[l];
	            }
	            std::size_t tail = 0;
	            __syncthreads();

	            for (std::size_t j = t / 32, k = 32; j > 0; j >>= 3, k *= 8){
	                std::size_t m_idx2 = t_idx / (k / 4);
	                std::size_t t_idx2 = t_idx % (k / 4);
	                for (std::size_t l = 0; l < 8; l++){
	                    samples[l] = buffer[set * n2 + 2 * m_idx2 * k + t_idx2 + (k / 4) * l];
	                }
	                tw_idx2 = j * tw_idx + m_idx2;
	                intt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
	                for (std::size_t l = 0; l < 8; l++){
	                    buffer[set * n2 + 2 * m_idx2 * k + t_idx2 + (k / 4) * l] = samples[l];
	                }
	                if (j == 2) tail = 1;
	                if (j == 4) tail = 2;
	                __syncthreads();
	            }

	            for (std::size_t j = 0; j < 8; j++){
	                samples[j] = buffer[set * n2 + t_idx + t / 4 * j];
	            }
	            if (tail == 1){
	                gs_butterfly(samples[0], samples[4], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	                gs_butterfly(samples[1], samples[5], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	                gs_butterfly(samples[2], samples[6], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	                gs_butterfly(samples[3], samples[7], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	            }
	            else if (tail == 2){
	                intt4(samples, psi, psi_shoup, tw_idx, modulus_value);
	                intt4(samples + 1, psi, psi_shoup, tw_idx, modulus_value);
	            }
	            for (std::size_t j = 0; j < 8; j++){
		                *(output_ptr + n_init + t / 4 * j) = samples[j];
	            }
	        }
	    }

	    template <NTTMode Mode>
	    __global__ void batched_inplace_inwt_radix8_phase2(
	        std::uint64_t* inout, const std::uint64_t* itwiddles,
	        const std::uint64_t* itwiddles_shoup,
	        const std::uint64_t* inv_degree_modulo,
	        const std::uint64_t* inv_degree_modulo_shoup,
	        const DModulus* modulus, std::size_t coeff_mod_size,
	        std::size_t start_mod_idx, std::size_t poly_count,
	        const int* order, std::size_t n, std::size_t n1,
	        std::size_t pad)
	    {
	        extern __shared__ std::uint64_t buffer[];

	        std::size_t pad_tid = threadIdx.x % pad;
	        std::size_t pad_idx = threadIdx.x / pad;
	        std::size_t group = n1 / 8;
	        std::uint64_t samples[8];
	        std::size_t t = n / 2;

	        for (std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	             i < (n / 8) * (coeff_mod_size * poly_count);
	             i += blockDim.x * gridDim.x)
	        {
	            std::size_t row_idx = i / (n / 8);
	            std::size_t data_row_idx = 0;
	            std::size_t twr_idx = 0;
                intt_table_access<Mode>(
                    row_idx, coeff_mod_size, start_mod_idx, order,
                    data_row_idx, twr_idx);
	            std::size_t n_idx = i % (n / 8);

	            std::uint64_t* data_ptr = inout + data_row_idx * n;
	            const std::uint64_t* psi = itwiddles + n * twr_idx;
	            const std::uint64_t* psi_shoup = itwiddles_shoup + n * twr_idx;
	            std::uint64_t modulus_value = modulus[twr_idx].value();
	            std::uint64_t inv_degree_mod = inv_degree_modulo[twr_idx];
	            std::uint64_t inv_degree_mod_shoup = inv_degree_modulo_shoup[twr_idx];
	            std::size_t n_init = 2 * t / group * pad_idx + pad_tid + pad * (n_idx / (group * pad));

	            for (std::size_t j = 0; j < 8; j++){
	                samples[j] = *(data_ptr + n_init + t / 4 / group * j);
	            }
	            std::size_t tw_idx = 1;
	            std::size_t tw_idx2 = group * tw_idx + pad_idx;
	            intt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
	            for (std::size_t j = 0; j < 8; j++){
	                buffer[pad_tid * (n1 + pad) + 8 * pad_idx + j] = samples[j];
	            }
	            std::size_t tail = 0;
	            __syncthreads();

	            for (std::size_t j = group / 8, k = 32; j > 0; j >>= 3, k *= 8){
	                std::size_t m_idx2 = pad_idx / (k / 4);
	                std::size_t t_idx2 = pad_idx % (k / 4);
	                for (std::size_t l = 0; l < 8; l++){
	                    samples[l] = buffer[(n1 + pad) * pad_tid + 2 * m_idx2 * k + t_idx2 + (k / 4) * l];
	                }
	                tw_idx2 = j * tw_idx + m_idx2;
	                intt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
	                for (std::size_t l = 0; l < 8; l++){
	                    buffer[(n1 + pad) * pad_tid + 2 * m_idx2 * k + t_idx2 + (k / 4) * l] = samples[l];
	                }
	                if (j == 2) tail = 1;
	                if (j == 4) tail = 2;
	                __syncthreads();
	            }
	            if (group < 8) tail = (group == 4) ? 2 : 1;

	            for (std::size_t l = 0; l < 8; l++){
	                samples[l] = buffer[pad_tid * (n1 + pad) + pad_idx + group * l];
	            }
	            if (tail == 1){
	                gs_butterfly(samples[0], samples[4], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	                gs_butterfly(samples[1], samples[5], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	                gs_butterfly(samples[2], samples[6], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	                gs_butterfly(samples[3], samples[7], psi[tw_idx], psi_shoup[tw_idx], modulus_value);
	            }
	            else if (tail == 2){
	                intt4(samples, psi, psi_shoup, tw_idx, modulus_value);
	                intt4(samples + 1, psi, psi_shoup, tw_idx, modulus_value);
	            }

	            for (std::size_t j = 0; j < 4; j++){
	                samples[j] = multiply_and_reduce_shoup_lazy(
	                    samples[j], inv_degree_mod, inv_degree_mod_shoup, modulus_value);
	            }

	            n_init = t / 4 / group * pad_idx + pad_tid + pad * (n_idx / (group * pad));
	            for (std::size_t j = 0; j < 8; j++){
	                csub_q(samples[j], modulus_value);
	                *(data_ptr + n_init + t / 4 * j) = samples[j];
	            }
	        }
	    }

    template <NTTMode Mode>
    void launch_ntt_radix8_batched(
        Data64* data, gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count, int start_mod_idx, const int* order,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        const int poly_count = batch_size / mod_count;
        const DNTTTable& ntt_tables = tables->table;
        const std::size_t n = ntt_tables.n();
        std::size_t phase1_sample_size = sample_size(n);
        const std::size_t phase2_sample_size = n / phase1_sample_size;
        constexpr std::size_t per_block_memory =
            blockDimNTT.x * per_thread_sample_size *
            sizeof(std::uint64_t);

        batched_inplace_fnwt_radix8_phase1<Mode><<<
            gridDimNTT, (phase1_sample_size / 8) * per_block_pad,
            (phase1_sample_size + per_block_pad + 1) * per_block_pad *
                sizeof(std::uint64_t),
            cfg.stream>>>(
            data, ntt_tables.twiddle(), ntt_tables.twiddle_shoup(),
            ntt_tables.modulus(), mod_count, start_mod_idx, poly_count,
            order, n, phase1_sample_size, per_block_pad);

        batched_inplace_fnwt_radix8_phase2<Mode><<<
            gridDimNTT, blockDimNTT, per_block_memory, cfg.stream>>>(
            data, ntt_tables.twiddle(), ntt_tables.twiddle_shoup(),
            ntt_tables.modulus(), mod_count, start_mod_idx, poly_count,
            order, n, phase1_sample_size, phase2_sample_size);
    }

    template <NTTMode Mode>
    void launch_intt_radix8_batched(
        const Data64* input, Data64* output,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count, int start_mod_idx, const int* order,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        const int poly_count = batch_size / mod_count;
        const DNTTTable& ntt_tables = tables->table;
        const std::size_t n = ntt_tables.n();
        std::size_t phase2_sample_size = sample_size(n);
        const std::size_t phase1_sample_size = n / phase2_sample_size;
        constexpr std::size_t per_block_memory =
            blockDimNTT.x * per_thread_sample_size * sizeof(std::uint64_t);

        batched_inplace_inwt_radix8_phase1<Mode><<<
            gridDimNTT, blockDimNTT, per_block_memory, cfg.stream>>>(
            input, output, ntt_tables.itwiddle(),
            ntt_tables.itwiddle_shoup(), ntt_tables.modulus(), mod_count,
            start_mod_idx, poly_count, order, n, phase1_sample_size,
            phase2_sample_size);

        batched_inplace_inwt_radix8_phase2<Mode><<<
            gridDimNTT, (phase1_sample_size / 8) * per_block_pad,
            (phase1_sample_size + per_block_pad + 1) * per_block_pad *
                sizeof(std::uint64_t),
            cfg.stream>>>(
            output, ntt_tables.itwiddle(), ntt_tables.itwiddle_shoup(),
            ntt_tables.n_inv_mod_q(), ntt_tables.n_inv_mod_q_shoup(),
            ntt_tables.modulus(), mod_count, start_mod_idx, poly_count,
            order, n, phase1_sample_size, per_block_pad);
    }

    __global__ void modmajor_exclude_ntt_phase1(
        Data64* data, const Data64* twiddles, const Data64* twiddles_shoup,
        const DModulus* ntt_modulus, const int* group_sizes,
        const int* group_locations, int current_qp, int current_q,
        int first_q_count, std::size_t n, std::size_t n1, std::size_t pad,
        int decomp_count, bool skip_excluded)
    {
        extern __shared__ Data64 buffer[];
        const std::size_t total_rows =
            static_cast<std::size_t>(current_qp) * decomp_count;

        for (std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
             tid < (n / 8) * total_rows; tid += blockDim.x * gridDim.x)
        {
            const int logical_row = static_cast<int>(tid / (n / 8));
            const int local_mod = logical_row / decomp_count;
            const int group_idx = logical_row - local_mod * decomp_count;
            const int group_size = group_sizes[group_idx];
            const int group_location = group_locations[group_idx];
            if (skip_excluded &&
                local_mod >= group_location &&
                local_mod < group_location + group_size)
            {
                continue;
            }

            const std::size_t pad_tid = threadIdx.x % pad;
            const std::size_t pad_idx = threadIdx.x / pad;
            const std::size_t group = n1 / 8;
            Data64 samples[8];
            const std::size_t t = n / 2;
            const std::size_t n_idx = tid % (n / 8);
            const std::size_t n_init =
                t / 4 / group * pad_idx + pad_tid +
                pad * (n_idx / (group * pad));
            Data64* data_ptr =
                data + static_cast<std::size_t>(logical_row) * n;

            const int table_idx =
                local_mod >= current_q ? first_q_count + local_mod - current_q
                                       : local_mod;
            const Data64* psi =
                twiddles + static_cast<std::size_t>(table_idx) * n;
            const Data64* psi_shoup =
                twiddles_shoup + static_cast<std::size_t>(table_idx) * n;
            const Data64 modulus_value = ntt_modulus[table_idx].value();

#pragma unroll
            for (std::size_t j = 0; j < 8; ++j)
            {
                samples[j] = data_ptr[n_init + t / 4 * j];
            }

            std::size_t tw_idx = 1;
            fntt8(samples, psi, psi_shoup, tw_idx, modulus_value);
#pragma unroll
            for (std::size_t j = 0; j < 8; ++j)
            {
                buffer[pad_tid * (n1 + pad) + pad_idx + group * j] =
                    samples[j];
            }
            std::size_t remain_iters = 0;
            __syncthreads();
#pragma unroll
            for (std::size_t j = 8, k = group / 2; j < group + 1;
                 j *= 8, k >>= 3)
            {
                const std::size_t m_idx2 = pad_idx / (k / 4);
                const std::size_t t_idx2 = pad_idx % (k / 4);
#pragma unroll
                for (std::size_t l = 0; l < 8; ++l)
                {
                    samples[l] =
                        buffer[(n1 + pad) * pad_tid + 2 * m_idx2 * k +
                               t_idx2 + (k / 4) * l];
                }
                const std::size_t tw_idx2 = j * tw_idx + m_idx2;
                fntt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
#pragma unroll
                for (std::size_t l = 0; l < 8; ++l)
                {
                    buffer[(n1 + pad) * pad_tid + 2 * m_idx2 * k + t_idx2 +
                           (k / 4) * l] = samples[l];
                }
                if (j == group / 2)
                {
                    remain_iters = 1;
                }
                if (j == group / 4)
                {
                    remain_iters = 2;
                }
                __syncthreads();
            }

            if (group < 8)
            {
                remain_iters = group == 4 ? 2 : 1;
            }
#pragma unroll
            for (std::size_t l = 0; l < 8; ++l)
            {
                samples[l] = buffer[(n1 + pad) * pad_tid + 8 * pad_idx + l];
            }
            if (remain_iters == 1)
            {
                const std::size_t tw_idx2 = 4 * group * tw_idx + 4 * pad_idx;
                ct_butterfly(samples[0], samples[1], psi[tw_idx2],
                             psi_shoup[tw_idx2], modulus_value);
                ct_butterfly(samples[2], samples[3], psi[tw_idx2 + 1],
                             psi_shoup[tw_idx2 + 1], modulus_value);
                ct_butterfly(samples[4], samples[5], psi[tw_idx2 + 2],
                             psi_shoup[tw_idx2 + 2], modulus_value);
                ct_butterfly(samples[6], samples[7], psi[tw_idx2 + 3],
                             psi_shoup[tw_idx2 + 3], modulus_value);
            }
            else if (remain_iters == 2)
            {
                const std::size_t tw_idx2 = 2 * group * tw_idx + 2 * pad_idx;
                fntt4(samples, psi, psi_shoup, tw_idx2, modulus_value);
                fntt4(samples + 4, psi, psi_shoup, tw_idx2 + 1,
                      modulus_value);
            }
#pragma unroll
            for (std::size_t l = 0; l < 8; ++l)
            {
                buffer[(n1 + pad) * pad_tid + 8 * pad_idx + l] = samples[l];
            }

            __syncthreads();
            for (std::size_t j = 0; j < 8; ++j)
            {
                data_ptr[n_init + t / 4 * j] =
                    buffer[pad_tid * (n1 + pad) + pad_idx + group * j];
            }
        }
    }

    __global__ void modmajor_exclude_ntt_phase2(
        Data64* data, const Data64* twiddles, const Data64* twiddles_shoup,
        const DModulus* modulus, const int* group_sizes,
        const int* group_locations, int current_qp, int current_q,
        int first_q_count, std::size_t n, std::size_t n1, std::size_t n2,
        int decomp_count, bool skip_excluded)
    {
        extern __shared__ Data64 buffer[];
        const std::size_t total_rows =
            static_cast<std::size_t>(current_qp) * decomp_count;

        for (std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
             tid < (n / 8) * total_rows; tid += blockDim.x * gridDim.x)
        {
            const int reverse_row =
                static_cast<int>(total_rows - 1 - tid / (n / 8));
            const int local_mod = reverse_row / decomp_count;
            const int group_idx = reverse_row - local_mod * decomp_count;
            const int group_size = group_sizes[group_idx];
            const int group_location = group_locations[group_idx];
            if (skip_excluded &&
                local_mod >= group_location &&
                local_mod < group_location + group_size)
            {
                continue;
            }

            const std::size_t group = n2 / 8;
            const std::size_t set = threadIdx.x / group;
            Data64 samples[8];
            const std::size_t t = n2 / 2;
            const int table_idx =
                local_mod >= current_q ? first_q_count + local_mod - current_q
                                       : local_mod;
            const std::size_t n_idx = tid % (n / 8);
            const std::size_t m_idx = n_idx / (t / 4);
            const std::size_t t_idx = n_idx % (t / 4);
            Data64* data_ptr = data + static_cast<std::size_t>(reverse_row) * n;
            const Data64* psi =
                twiddles + static_cast<std::size_t>(table_idx) * n;
            const Data64* psi_shoup =
                twiddles_shoup + static_cast<std::size_t>(table_idx) * n;
            const Data64 modulus_value = modulus[table_idx].value();
            const std::size_t n_init = 2 * m_idx * t + t_idx;
#pragma unroll
            for (std::size_t j = 0; j < 8; ++j)
            {
                samples[j] = *(data_ptr + n_init + t / 4 * j);
            }
            const std::size_t tw_idx = n1 + m_idx;
            fntt8(samples, psi, psi_shoup, tw_idx, modulus_value);
#pragma unroll
            for (std::size_t j = 0; j < 8; ++j)
            {
                buffer[set * n2 + t_idx + t / 4 * j] = samples[j];
            }
            std::size_t tail = 0;
            __syncthreads();

#pragma unroll
            for (std::size_t j = 8, k = t / 8; j < t / 4 + 1;
                 j *= 8, k >>= 3)
            {
                const std::size_t m_idx2 = t_idx / (k / 4);
                const std::size_t t_idx2 = t_idx % (k / 4);
#pragma unroll
                for (std::size_t l = 0; l < 8; ++l)
                {
                    samples[l] =
                        buffer[set * n2 + 2 * m_idx2 * k + t_idx2 +
                               (k / 4) * l];
                }
                const std::size_t tw_idx2 = j * tw_idx + m_idx2;
                fntt8(samples, psi, psi_shoup, tw_idx2, modulus_value);
#pragma unroll
                for (std::size_t l = 0; l < 8; ++l)
                {
                    buffer[set * n2 + 2 * m_idx2 * k + t_idx2 +
                           (k / 4) * l] = samples[l];
                }
                if (j == t / 8)
                {
                    tail = 1;
                }
                if (j == t / 16)
                {
                    tail = 2;
                }
                __syncthreads();
            }

#pragma unroll
            for (std::size_t l = 0; l < 8; ++l)
            {
                samples[l] = buffer[set * n2 + 8 * t_idx + l];
            }
            if (tail == 1)
            {
                const std::size_t tw_idx2 = t * tw_idx + 4 * t_idx;
                ct_butterfly(samples[0], samples[1], psi[tw_idx2],
                             psi_shoup[tw_idx2], modulus_value);
                ct_butterfly(samples[2], samples[3], psi[tw_idx2 + 1],
                             psi_shoup[tw_idx2 + 1], modulus_value);
                ct_butterfly(samples[4], samples[5], psi[tw_idx2 + 2],
                             psi_shoup[tw_idx2 + 2], modulus_value);
                ct_butterfly(samples[6], samples[7], psi[tw_idx2 + 3],
                             psi_shoup[tw_idx2 + 3], modulus_value);
            }
            else if (tail == 2)
            {
                const std::size_t tw_idx2 = (t / 2) * tw_idx + 2 * t_idx;
                fntt4(samples, psi, psi_shoup, tw_idx2, modulus_value);
                fntt4(samples + 4, psi, psi_shoup, tw_idx2 + 1,
                      modulus_value);
            }
#pragma unroll
            for (std::size_t l = 0; l < 8; ++l)
            {
                buffer[set * n2 + 8 * t_idx + l] = samples[l];
            }
            __syncthreads();

            const Data64 modulus2 = modulus_value << 1;
#pragma unroll
            for (std::size_t j = 0; j < 8; ++j)
            {
                samples[j] = buffer[set * n2 + t_idx + t / 4 * j];
                csub_q(samples[j], modulus2);
                csub_q(samples[j], modulus_value);
            }
#pragma unroll
            for (std::size_t j = 0; j < 8; ++j)
            {
                *(data_ptr + n_init + t / 4 * j) = samples[j];
            }
        }
    }

    void launch_ntt_modmajor_exclude(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int current_qp,
        int current_q,
        int first_q_count,
        int decomp_count,
        const int* group_sizes,
        const int* group_locations,
        const std::shared_ptr<PhantomNttTables>& tables,
        bool skip_excluded)
    {
        const DNTTTable& ntt_tables = tables->table;
        const std::size_t n = ntt_tables.n();
        const std::size_t phase1_sample_size = sample_size(n);
        const std::size_t phase2_sample_size = n / phase1_sample_size;
        const std::size_t phase1_shared =
            (phase1_sample_size + per_block_pad + 1) * per_block_pad *
            sizeof(Data64);
        const std::size_t phase2_shared =
            blockDimNTT.x * per_thread_sample_size * sizeof(Data64);
        constexpr unsigned int max_grid_x = 65535;
        const unsigned int grid_scale =
            decomp_count > 1 ? static_cast<unsigned int>(decomp_count) : 1U;
        const unsigned int grid_blocks =
            gridDimNTT.x > max_grid_x / grid_scale ? max_grid_x
                                                        : gridDimNTT.x *
                                                              grid_scale;
        const dim3 modmajor_grid(
            grid_blocks, gridDimNTT.y, gridDimNTT.z);

        modmajor_exclude_ntt_phase1<<<
            modmajor_grid, (phase1_sample_size / 8) * per_block_pad,
            phase1_shared, cfg.stream>>>(
            data, ntt_tables.twiddle(), ntt_tables.twiddle_shoup(),
            ntt_tables.modulus(), group_sizes, group_locations, current_qp,
            current_q, first_q_count, n, phase1_sample_size, per_block_pad,
            decomp_count, skip_excluded);

        modmajor_exclude_ntt_phase2<<<modmajor_grid, blockDimNTT,
                                      phase2_shared, cfg.stream>>>(
            data, ntt_tables.twiddle(), ntt_tables.twiddle_shoup(),
            ntt_tables.modulus(), group_sizes, group_locations, current_qp,
            current_q, first_q_count, n, phase1_sample_size,
            phase2_sample_size, decomp_count, skip_excluded);
    }

} // namespace

    std::shared_ptr<PhantomNttTables> make_phantom_ntt_tables_from_heongpu_roots(
        const std::vector<Modulus64>& moduli,
        const std::vector<Root64>& forward_roots,
        const std::vector<Root64>& inverse_roots,
        const std::vector<Ninverse64>& n_inverse,
        int n_power,
        cudaStream_t stream,
        bool forward_only)
    {
        auto out = std::make_shared<PhantomNttTables>();
        const std::size_t n = std::size_t{1} << n_power;
        out->table.init(n, moduli.size(), stream);

        for (std::size_t i = 0; i < moduli.size(); ++i){
            DModulus dmod = make_dmodulus(moduli[i].value);

            const std::size_t offset = i * n;
            std::vector<std::uint64_t> fwd(forward_roots.begin() + offset,
                                           forward_roots.begin() + offset + n);
            auto fwd_shoup = make_shoup_table(fwd, moduli[i].value);

            if (forward_only)   //generate forward table only for sparse NTT.
            {
                cudaMemcpyAsync(out->table.modulus() + i, &dmod,
                                sizeof(DModulus), cudaMemcpyHostToDevice,
                                stream);
                cudaMemcpyAsync(out->table.twiddle() + offset, fwd.data(),
                                n * sizeof(std::uint64_t),
                                cudaMemcpyHostToDevice, stream);
                cudaMemcpyAsync(out->table.twiddle_shoup() + offset,
                                fwd_shoup.data(),
                                n * sizeof(std::uint64_t),
                                cudaMemcpyHostToDevice, stream);
                continue;
            }

            std::vector<std::uint64_t> inv(
                inverse_roots.begin() + offset,
                inverse_roots.begin() + offset + n);

            Modulus64 heon_mod(moduli[i].value);
            if (n > 1)
                inv[1] = OPERATOR64::mult(inv[1], n_inverse[i], heon_mod);

            auto inv_shoup = make_shoup_table(inv, moduli[i].value);
            std::uint64_t n_inv_shoup =
                compute_shoup(n_inverse[i], moduli[i].value);

            out->table.set(&dmod, fwd.data(), fwd_shoup.data(), inv.data(), inv_shoup.data(), n_inverse[i], n_inv_shoup, i, stream);
        }

        return out;
    }

    void phantom_ntt_inplace(
        Data64* data, gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count,
        const std::shared_ptr<PhantomNttTables>& tables, bool batch)
    {
        if (!batch)
        {
            for (int poly = 0; poly < batch_size / mod_count; ++poly)
            {
                launch_ntt_radix8_batched<NTTMode::Inplace>(
                    data + std::size_t(poly) * mod_count * tables->table.n(),
                    cfg, mod_count, mod_count, 0, nullptr, tables);
            }
            return;
        }

        launch_ntt_radix8_batched<NTTMode::Inplace>(
            data, cfg, batch_size, mod_count, 0, nullptr, tables);
    }

    void phantom_ntt_modulus_ordered_inplace(
        Data64* data, gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count, const int* order,
        const std::shared_ptr<PhantomNttTables>& tables, bool batch)
    {
        if (!batch)
        {
            for (int poly = 0; poly < batch_size / mod_count; ++poly)
            {
                phantom_ntt_modulus_ordered_inplace(
                    data + std::size_t(poly) * mod_count * tables->table.n(),
                    cfg, mod_count, mod_count, order, tables, true);
            }
            return;
        }

        launch_ntt_radix8_batched<NTTMode::ModulusOrdered>(
            data, cfg, batch_size, mod_count, 0, order, tables);
    }

    void phantom_ntt_modmajor_inplace(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int current_qp,
        int current_q,
        int first_q_count,
        int decomp_count,
        const int* group_sizes,
        const int* group_locations,
        const std::shared_ptr<PhantomNttTables>& tables,
        bool skip_excluded)
    {
        launch_ntt_modmajor_exclude(
            data, cfg, current_qp, current_q, first_q_count, decomp_count,
            group_sizes, group_locations, tables, skip_excluded);
    }

    void phantom_intt_inplace(
        Data64* data, gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count,
        const std::shared_ptr<PhantomNttTables>& tables, bool batch)
    {
        if (!batch)
        {
            for (int poly = 0; poly < batch_size / mod_count; ++poly)
            {
                launch_intt_radix8_batched<NTTMode::Inplace>(
                    data + std::size_t(poly) * mod_count * tables->table.n(),
                    data + std::size_t(poly) * mod_count * tables->table.n(),
                    cfg, mod_count, mod_count, 0, nullptr, tables);
            }
            return;
        }

        launch_intt_radix8_batched<NTTMode::Inplace>(
            data, data, cfg, batch_size, mod_count, 0, nullptr, tables);
    }

    void phantom_intt(
        const Data64* input, Data64* output,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count,
        const std::shared_ptr<PhantomNttTables>& tables, bool batch)
    {
        if (!batch)
        {
            for (int poly = 0; poly < batch_size / mod_count; ++poly)
            {
                launch_intt_radix8_batched<NTTMode::Inplace>(
                    input + std::size_t(poly) * mod_count * tables->table.n(),
                    output +
                        std::size_t(poly) * mod_count * tables->table.n(),
                    cfg, mod_count, mod_count, 0, nullptr, tables);
            }
            return;
        }

        launch_intt_radix8_batched<NTTMode::Inplace>(
            input, output, cfg, batch_size, mod_count, 0, nullptr, tables);
    }

    void phantom_intt_modulus_ordered_inplace(
        Data64* data, gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count, const int* order,
        const std::shared_ptr<PhantomNttTables>& tables, bool batch)
    {
        if (!batch)
        {
            for (int poly = 0; poly < batch_size / mod_count; ++poly)
            {
                phantom_intt_modulus_ordered_inplace(
                    data + std::size_t(poly) * mod_count * tables->table.n(),
                    cfg, mod_count, mod_count, order, tables, true);
            }
            return;
        }

        launch_intt_radix8_batched<NTTMode::ModulusOrdered>(
            data, data, cfg, batch_size, mod_count, 0, order, tables);
    }

    void phantom_intt_poly_ordered_inplace_batched(
        Data64* data, gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size, int mod_count, const int* order, int start_mod_idx,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        launch_intt_radix8_batched<NTTMode::PolyOrdered>(
            data, data, cfg, batch_size, mod_count, start_mod_idx, order,
            tables);
    }

} // namespace ntt
} // namespace heongpu
