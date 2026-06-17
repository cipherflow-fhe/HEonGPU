// Copyright 2025-2026 Yanbin Li
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Yanbin Li

#include <heongpu/fft/fft.cuh>

#include <stdexcept>

namespace heongpu
{
namespace fft
{
    using gpufft::CudaException;

    template <typename T>
    __global__ void Small_Special_InverseCore(COMPLEX<T>* polynomial,
                                              COMPLEX<T>* root_of_unity_table,
                                              int shared_index, int logm,
                                              int outer_iteration_count,
                                              int N_power,
                                              COMPLEX<T> n_inverse,
                                              bool not_last_kernel)
    {
        const int idx_x = threadIdx.x;
        const int idx_y = threadIdx.y;
        const int block_x = blockIdx.x;
        const int block_y = blockIdx.y;
        const int block_z = blockIdx.z;

        extern __shared__ char shared_memory_typed[];
        COMPLEX<T>* shared_memory =
            reinterpret_cast<COMPLEX<T>*>(shared_memory_typed);

        gpufft::location_t offset =
            static_cast<gpufft::location_t>(1) << (N_power - logm - 1);
        int t_ = shared_index; // @company CipherFlow: allow n_power < 9.
        int loops = outer_iteration_count;
        gpufft::location_t m =
            static_cast<gpufft::location_t>(1) << (N_power - logm - 1);

        gpufft::location_t global_addresss =
            idx_x +
            static_cast<gpufft::location_t>(
                idx_y * (offset / (1 << (outer_iteration_count - 1)))) +
            static_cast<gpufft::location_t>(blockDim.x * block_x) +
            static_cast<gpufft::location_t>(2 * block_y * offset) +
            (static_cast<gpufft::location_t>(block_z) << N_power);
        gpufft::location_t omega_addresss =
            idx_x +
            static_cast<gpufft::location_t>(
                idx_y * (offset / (1 << (outer_iteration_count - 1)))) +
            static_cast<gpufft::location_t>(blockDim.x * block_x) +
            static_cast<gpufft::location_t>(block_y * offset);

        gpufft::location_t shared_addresss = idx_x + (idx_y * blockDim.x);

        shared_memory[shared_addresss] = polynomial[global_addresss];
        shared_memory[shared_addresss + (blockDim.x * blockDim.y)] =
            polynomial[global_addresss + offset];

        int t = 1 << t_;
        int in_shared_address =
            ((shared_addresss >> t_) << t_) + shared_addresss;

        gpufft::location_t current_root_index;
#pragma unroll
        for (int lp = 0; lp < loops; lp++)
        {
            current_root_index = m + (omega_addresss & (m - 1));

            gpufft::GentlemanSandeUnit(
                shared_memory[in_shared_address],
                shared_memory[in_shared_address + t],
                root_of_unity_table[current_root_index]);

            t = t >> 1;
            t_ -= 1;
            m >>= 1;

            if (lp + 1 < loops)
            {
                in_shared_address =
                    ((shared_addresss >> t_) << t_) + shared_addresss;
                __syncthreads();
            }
        }
        __syncthreads();

        if (not_last_kernel)
        {
            polynomial[global_addresss] = shared_memory[shared_addresss];
            polynomial[global_addresss + offset] =
                shared_memory[shared_addresss + (blockDim.x * blockDim.y)];
        }
        else
        {
            polynomial[global_addresss] =
                shared_memory[shared_addresss] * n_inverse;
            polynomial[global_addresss + offset] =
                shared_memory[shared_addresss + (blockDim.x * blockDim.y)] *
                n_inverse;
        }
    }

    template <typename T>
    __host__ void GPU_Special_FFT(COMPLEX<T>* device_inout,
                                  COMPLEX<T>* root_of_unity_table,
                                  gpufft::fft_configuration<T> cfg,
                                  int batch_size)
    {
        if (cfg.n_power < 1 || cfg.n_power > 10)
        {
            throw std::invalid_argument(
                "heongpu::fft small Special FFT currently implements n_power=1..10");
        }

        switch (cfg.fft_type)
        {
            case gpufft::FORWARD:
            {
                auto kernel_parameters = CreateForwardSpecialFFTKernel<T>();

                for (const auto& params : kernel_parameters[cfg.n_power])
                {
                    gpufft::Special_ForwardCore<<<
                        dim3(params.griddim_x, params.griddim_y, batch_size),
                        dim3(params.blockdim_x, params.blockdim_y),
                        params.shared_memory, cfg.stream>>>(
                        device_inout, root_of_unity_table, params.shared_index,
                        params.logm, params.k, params.outer_iteration_count,
                        cfg.n_power);
                    GPUFFT_CUDA_CHECK(cudaGetLastError());
                }
                break;
            }
            case gpufft::INVERSE:
            {
                auto kernel_parameters = CreateInverseSpecialFFTKernel<T>();

                for (int i = 0; i < kernel_parameters[cfg.n_power].size() - 1;
                     i++)
                {
                    auto& params = kernel_parameters[cfg.n_power][i];
                    if (cfg.n_power < 9) // @company CipherFlow
                    {
                        Small_Special_InverseCore<<<
                            dim3(params.griddim_x, params.griddim_y,
                                 batch_size),
                            dim3(params.blockdim_x, params.blockdim_y),
                            params.shared_memory, cfg.stream>>>(
                            device_inout, root_of_unity_table,
                            params.shared_index, params.logm,
                            params.outer_iteration_count, cfg.n_power,
                            cfg.mod_inverse, true);
                    }
                    else
                    {
                        gpufft::Special_InverseCore<<<
                            dim3(params.griddim_x, params.griddim_y,
                                 batch_size),
                            dim3(params.blockdim_x, params.blockdim_y),
                            params.shared_memory, cfg.stream>>>(
                            device_inout, root_of_unity_table,
                            params.shared_index, params.logm,
                            params.outer_iteration_count, cfg.n_power,
                            cfg.mod_inverse, true);
                    }
                    GPUFFT_CUDA_CHECK(cudaGetLastError());
                }

                auto& params =
                    kernel_parameters[cfg.n_power]
                                     [kernel_parameters[cfg.n_power].size() - 1];
                if (cfg.n_power < 9) // @company CipherFlow
                {
                    Small_Special_InverseCore<<<
                        dim3(params.griddim_x, params.griddim_y, batch_size),
                        dim3(params.blockdim_x, params.blockdim_y),
                        params.shared_memory, cfg.stream>>>(
                        device_inout, root_of_unity_table, params.shared_index,
                        params.logm, params.outer_iteration_count, cfg.n_power,
                        cfg.mod_inverse, false);
                }
                else
                {
                    gpufft::Special_InverseCore<<<
                        dim3(params.griddim_x, params.griddim_y, batch_size),
                        dim3(params.blockdim_x, params.blockdim_y),
                        params.shared_memory, cfg.stream>>>(
                        device_inout, root_of_unity_table, params.shared_index,
                        params.logm, params.outer_iteration_count, cfg.n_power,
                        cfg.mod_inverse, false);
                }
                GPUFFT_CUDA_CHECK(cudaGetLastError());
                break;
            }
            default:
                break;
        }
    }

    template __host__ void GPU_Special_FFT(
        Complex32* device_inout, Complex32* root_of_unity_table,
        gpufft::fft_configuration<Float32> cfg, int batch_size);

    template __host__ void GPU_Special_FFT(
        Complex64* device_inout, Complex64* root_of_unity_table,
        gpufft::fft_configuration<Float64> cfg, int batch_size);
} // namespace fft
} // namespace heongpu
