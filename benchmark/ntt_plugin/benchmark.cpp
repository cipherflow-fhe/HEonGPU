// Copyright 2024-2026 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include "common.hpp"

#include <heongpu/ntt/ntt.cuh>
#include <heongpu/primitive/ntt.cuh>
#include <heongpu/primitive/switchkey.cuh>
#include <heongpu/util/util.cuh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using ntt_plugin::build_key_modulus;
using ntt_plugin::make_context;
using ntt_plugin::make_bootstrapping_config;
using ntt_plugin::make_parameters;
using ntt_plugin::make_ntt_order;
using ntt_plugin::make_bsgs_input;
using ntt_plugin::ntt_surface_title;
using ntt_plugin::set_optimization_env;
using ntt_plugin::sort_parameters_by_n;
using ntt_plugin::NttSurface;
using ntt_plugin::ParameterSet;
using ntt_plugin::Scheme;
using heongpu::CudaException;

constexpr int kKernelPolyCounts[] = {1, 2, 4, 8, 16, 32};
constexpr int kWarmupCount = 2;
constexpr int kDefaultRepeatCount = 10;
constexpr int kDefaultBootstrapRepeatCount = 20;

struct OperationTiming
{
    float relinearize = 0.0F;
    float rotate = 0.0F;
    float bootstrapping = 0.0F;
};

struct BenchmarkResult
{
    std::string label;
    size_t poly_modulus_degree = 0;
    OperationTiming timing;
};

struct OptimizationRoute
{
    std::string name;
    bool use_phantom_ntt = false;
    bool use_mod_keyswitch = false;
    bool use_keyswitch_part2 = true;
    bool use_bsgs_fusion = true;
};

struct BootstrapRouteResult
{
    std::string name;
    BenchmarkResult result;
};

struct NttKernelResult
{
    std::string label;
    size_t poly_modulus_degree = 0;
    int poly_count = 0;
    float gpuntt_ms = 0.0F;
    float phantom_ms = 0.0F;
    float phantom_no_batch_ms = 0.0F;
};

template <typename Func>
void time_call(cudaEvent_t start_time, cudaEvent_t stop_time, float& total,
               Func&& func)
{
    cudaEventRecord(start_time);
    func();
    cudaEventRecord(stop_time);

    cudaEventSynchronize(stop_time);

    float elapsed = 0.0F;
    cudaEventElapsedTime(&elapsed, start_time, stop_time);
    total += elapsed;
}

template <typename Func>
float average_kernel_time(cudaStream_t stream, cudaEvent_t start_time,
                          cudaEvent_t stop_time, int warmup_count,
                          int repeat_count, Func&& func)
{
    for (int trial = 0; trial < warmup_count; ++trial)
    {
        func();
    }
    cudaStreamSynchronize(stream);

    float total = 0.0F;
    for (int trial = 0; trial < repeat_count; ++trial)
    {
        cudaEventRecord(start_time, stream);
        func();
        cudaEventRecord(stop_time, stream);
        cudaEventSynchronize(stop_time);

        float elapsed = 0.0F;
        cudaEventElapsedTime(&elapsed, start_time, stop_time);
        total += elapsed;
    }
    return total / repeat_count;
}

template <typename PrepareFunc, typename TimedFunc>
float average_prepared_kernel_time(cudaStream_t stream, cudaEvent_t start_time,
                                   cudaEvent_t stop_time, int warmup_count,
                                   int repeat_count, PrepareFunc&& prepare,
                                   TimedFunc&& func)
{
    for (int trial = 0; trial < warmup_count; ++trial)
    {
        prepare();
        cudaEventRecord(start_time, stream);
        func();
        cudaEventRecord(stop_time, stream);
        cudaEventSynchronize(stop_time);
    }
    cudaStreamSynchronize(stream);

    float total = 0.0F;
    for (int trial = 0; trial < repeat_count; ++trial)
    {
        prepare();
        cudaEventRecord(start_time, stream);
        func();
        cudaEventRecord(stop_time, stream);
        cudaEventSynchronize(stop_time);

        float elapsed = 0.0F;
        cudaEventElapsedTime(&elapsed, start_time, stop_time);
        total += elapsed;
    }
    return total / repeat_count;
}

BenchmarkResult run_operator_case(const ParameterSet& parameter,
                                  bool use_original_route)
{
    const bool use_phantom_ntt = !use_original_route;
    set_optimization_env(!use_original_route, !use_original_route,
                         !use_original_route);
    heongpu::HEContext<Scheme> context =
        make_context(parameter, use_phantom_ntt);

    heongpu::HEKeyGenerator<Scheme> keygen(context);
    heongpu::Secretkey<Scheme> secret_key(context);
    keygen.generate_secret_key(secret_key);

    heongpu::Publickey<Scheme> public_key(context);
    keygen.generate_public_key(public_key, secret_key);

    heongpu::Relinkey<Scheme> relin_key(context);
    keygen.generate_relin_key(relin_key, secret_key);

    std::vector<int> custom_key_index = {1};
    heongpu::Galoiskey<Scheme> galois_key(context, custom_key_index);
    keygen.generate_galois_key(galois_key, secret_key);

    heongpu::HEEncoder<Scheme> encoder(context);
    heongpu::HEEncryptor<Scheme> encryptor(context, public_key);
    heongpu::HEArithmeticOperator<Scheme> operators(context, encoder);

    const int row_size = parameter.poly_modulus_degree / 2;
    heongpu::HostVector<double> message(row_size, 1);

    OperationTiming timing;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    auto run_trial = [&](OperationTiming* measured_timing) {
        heongpu::Plaintext<Scheme> plaintext(context);
        encoder.encode(plaintext, message, parameter.scale);

        heongpu::Ciphertext<Scheme> c1(context);
        encryptor.encrypt(c1, plaintext);

        {
            heongpu::Ciphertext<Scheme> operand(context);
            operators.multiply(c1, c1, operand);
            if (measured_timing)
            {
                time_call(start_time, stop_time, measured_timing->relinearize,
                          [&]() {
                              operators.relinearize_inplace(operand,
                                                            relin_key);
                          });
            }
            else
            {
                operators.relinearize_inplace(operand, relin_key);
            }
        }

        {
            heongpu::Ciphertext<Scheme> operand = c1;
            if (measured_timing)
            {
                time_call(start_time, stop_time, measured_timing->rotate,
                          [&]() {
                              operators.rotate_rows(operand, operand,
                                                    galois_key, 1);
                          });
            }
            else
            {
                operators.rotate_rows(operand, operand, galois_key, 1);
            }
        }

        cudaDeviceSynchronize();
    };

    for (int trial = 0; trial < kWarmupCount; trial++)
    {
        run_trial(nullptr);
    }

    for (int trial = 0; trial < 100; trial++)
    {
        run_trial(&timing);
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);

    timing.relinearize /= 100;
    timing.rotate /= 100;

    return {parameter.label, parameter.poly_modulus_degree, timing};
}

BenchmarkResult run_bootstrapping_case(const ParameterSet& parameter,
                                       const OptimizationRoute& route,
                                       int repeat_count)
{
    set_optimization_env(route.use_mod_keyswitch, route.use_keyswitch_part2,
                         route.use_bsgs_fusion);
    heongpu::HEContext<Scheme> context =
        make_context(parameter, route.use_phantom_ntt, true);

    constexpr int secret_weight = 192;
    constexpr int ephemeral_secret_weight = 32;
    heongpu::HEKeyGenerator<Scheme> keygen(context);
    heongpu::Secretkey<Scheme> secret_key(context, secret_weight);
    keygen.generate_secret_key_v2(secret_key);

    heongpu::Publickey<Scheme> public_key(context);
    keygen.generate_public_key(public_key, secret_key);

    heongpu::Relinkey<Scheme> relin_key(context);
    keygen.generate_relin_key(relin_key, secret_key);

    heongpu::Secretkey<Scheme> sparse_secret_key(context,
                                                 ephemeral_secret_weight);
    keygen.generate_secret_key_v2(sparse_secret_key);

    heongpu::Switchkey<Scheme> swk_dense_to_sparse(context);
    keygen.generate_switch_key(swk_dense_to_sparse, sparse_secret_key,
                               secret_key);

    heongpu::Switchkey<Scheme> swk_sparse_to_dense(context);
    keygen.generate_switch_key(swk_sparse_to_dense, secret_key,
                               sparse_secret_key);

    heongpu::HEEncoder<Scheme> encoder(context);
    heongpu::HEEncryptor<Scheme> encryptor(context, public_key);
    heongpu::HEArithmeticOperator<Scheme> operators(context, encoder);

    const int slot_count = static_cast<int>(parameter.poly_modulus_degree / 2);
    std::vector<Complex64> message(slot_count, Complex64(0.2, 0.4));
    heongpu::Plaintext<Scheme> plaintext(context);
    encoder.encode(plaintext, message, parameter.scale);

    heongpu::Ciphertext<Scheme> boot_input_base(context);
    encryptor.encrypt(boot_input_base, plaintext);

    const auto boot_config =
        make_bootstrapping_config(context->get_key_modulus()[0].value);
    operators.generate_bootstrapping_params_v2(parameter.scale, boot_config);

    std::vector<int> key_index = operators.bootstrapping_key_indexs();
    heongpu::Galoiskey<Scheme> galois_key(context, key_index);
    keygen.generate_galois_key(galois_key, secret_key);

    for (int level = 1; level < static_cast<int>(parameter.q.size()); ++level)
    {
        operators.mod_drop_inplace(boot_input_base);
    }

    OperationTiming timing;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    for (int trial = 0; trial < repeat_count; ++trial)
    {
        heongpu::Ciphertext<Scheme> boot_input = boot_input_base;
        heongpu::Ciphertext<Scheme> boot_output(context);
        time_call(start_time, stop_time, timing.bootstrapping, [&]() {
            boot_output = operators.regular_bootstrapping_v2(
                boot_input, galois_key, relin_key, &swk_dense_to_sparse,
                &swk_sparse_to_dense);
        });
        cudaDeviceSynchronize();
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);

    timing.bootstrapping /= repeat_count;
    return {parameter.label, parameter.poly_modulus_degree, timing};
}

template <typename GpunttFunc, typename PhantomFunc, typename PhantomNoBatchFunc>
void measure_kernel_paths(cudaStream_t stream, cudaEvent_t start_time,
                          cudaEvent_t stop_time, int repeat_count,
                          std::size_t order_seed, GpunttFunc&& run_gpuntt,
                          PhantomFunc&& run_phantom,
                          PhantomNoBatchFunc&& run_phantom_no_batch,
                          float& gpuntt_ms, float& phantom_ms,
                          float& phantom_no_batch_ms)
{
    auto measure_gpuntt = [&]() {
        return average_kernel_time(
            stream, start_time, stop_time, kWarmupCount,
            repeat_count, run_gpuntt);
    };
    auto measure_phantom = [&]() {
        return average_kernel_time(
            stream, start_time, stop_time, kWarmupCount,
            repeat_count, run_phantom);
    };
    auto measure_phantom_no_batch = [&]() {
        return average_kernel_time(
            stream, start_time, stop_time, kWarmupCount,
            repeat_count, run_phantom_no_batch);
    };

    enum { GpunttPath, PhantomPath, PhantomNoBatchPath };
    const int path_orders[] = {
        GpunttPath, PhantomPath, PhantomNoBatchPath,
        PhantomPath, PhantomNoBatchPath, GpunttPath,
        PhantomNoBatchPath, GpunttPath, PhantomPath};
    const int* path_order = &path_orders[(order_seed % 3) * 3];

    for (int i = 0; i < 3; ++i)
    {
        switch (path_order[i])
        {
            case GpunttPath:
                gpuntt_ms = measure_gpuntt();
                break;
            case PhantomPath:
                phantom_ms = measure_phantom();
                break;
            default:
                phantom_no_batch_ms = measure_phantom_no_batch();
                break;
        }
    }
}

std::vector<NttKernelResult> run_ntt_kernel_comparison(
    const std::vector<ParameterSet>& parameters, int repeat_count,
    NttSurface mode)
{
    std::vector<NttKernelResult> results;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    for (const auto& parameter : parameters)
    {
        if (parameter.poly_modulus_degree < 4096 ||
            parameter.poly_modulus_degree > 131072)
        {
            continue;
        }

        const auto moduli = build_key_modulus(parameter);
        const int total_mod_count = static_cast<int>(moduli.size());
        const int mod_count =
            mode == NttSurface::PolyOrderedInverse ? 1 : total_mod_count;
        const int start_mod_idx =
            mode == NttSurface::PolyOrderedInverse
                ? std::max(0, total_mod_count - 1)
                : 0;
        const int n_power =
            static_cast<int>(std::log2(parameter.poly_modulus_degree));
        const auto roots_base = heongpu::generate_primitive_root_of_unity(
            parameter.poly_modulus_degree, moduli);
        const auto forward_roots =
            heongpu::generate_ntt_table(roots_base, moduli, n_power);
        const auto inverse_roots =
            heongpu::generate_intt_table(roots_base, moduli, n_power);
        const auto n_inverse = heongpu::generate_n_inverse(
            parameter.poly_modulus_degree, moduli);

        heongpu::DeviceVector<Modulus64> device_moduli(moduli);
        heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);
        heongpu::DeviceVector<Root64> device_inverse_roots(inverse_roots);
        heongpu::DeviceVector<Ninverse64> device_n_inverse(n_inverse);
        auto phantom_tables =
            heongpu::ntt::make_phantom_ntt_tables_from_heongpu_roots(
                moduli, forward_roots, inverse_roots, n_inverse, n_power,
                cudaStreamLegacy, mode == NttSurface::ModulusOrderedForward);

        for (int poly_count : kKernelPolyCounts)
        {
            const int batch_size = mod_count * poly_count;
            const std::size_t element_count =
                parameter.poly_modulus_degree *
                static_cast<std::size_t>(batch_size);
            std::size_t free_memory = 0;
            std::size_t total_memory = 0;
            cudaMemGetInfo(&free_memory, &total_memory);
            const std::size_t buffer_count = 4;
            if ((buffer_count * element_count * sizeof(Data64)) >
                static_cast<std::size_t>(free_memory * 0.70))
            {
                continue;
            }

            std::vector<int> order =
                make_ntt_order(mode, mod_count, batch_size);
            std::vector<int> device_order_values = order.empty()
                                                       ? std::vector<int>{0}
                                                       : order;
            heongpu::DeviceVector<int> device_order(device_order_values);

            cudaStream_t stream;
            cudaStreamCreate(&stream);

            gpuntt::ntt_rns_configuration<Data64> cfg = {
                .n_power = n_power,
                .ntt_type = mode == NttSurface::ForwardInplace ||
                            mode == NttSurface::ModulusOrderedForward
                        ? gpuntt::FORWARD
                        : gpuntt::INVERSE,
                .ntt_layout = gpuntt::PerPolynomial,
                .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
                .zero_padding = false,
                .mod_inverse = device_n_inverse.data() + start_mod_idx,
                .stream = stream};

            float gpuntt_ms = 0.0F;
            float phantom_ms = 0.0F;
            float phantom_no_batch_ms = 0.0F;

            heongpu::DeviceVector<Data64> buffer0(element_count);
            heongpu::DeviceVector<Data64> buffer1(element_count);
            heongpu::DeviceVector<Data64> buffer2(element_count);
            heongpu::DeviceVector<Data64> buffer3(element_count);
            cudaMemsetAsync(buffer0.data(), 0,
                            element_count * sizeof(Data64), stream);
            cudaMemsetAsync(buffer1.data(), 0,
                            element_count * sizeof(Data64), stream);
            cudaMemsetAsync(buffer2.data(), 0,
                            element_count * sizeof(Data64), stream);
            cudaMemsetAsync(buffer3.data(), 0,
                            element_count * sizeof(Data64), stream);
            cudaStreamSynchronize(stream);

            if (mode == NttSurface::InverseOutOfPlace)
            {
                Data64* input_data = buffer0.data();
                Data64* gpuntt_output = buffer1.data();
                Data64* phantom_output = buffer2.data();
                Data64* phantom_no_batch_output = buffer3.data();

                auto run_gpuntt = [&]() {
                    gpuntt::GPU_INTT(
                        input_data, gpuntt_output,
                        device_inverse_roots.data(), device_moduli.data(), cfg,
                        batch_size, mod_count);
                };
                auto run_phantom = [&]() {
                    heongpu::ntt::phantom_intt(
                        input_data, phantom_output, cfg, batch_size, mod_count,
                        phantom_tables, true);
                };

                auto run_phantom_no_batch = [&]() {
                    heongpu::ntt::phantom_intt(
                        input_data, phantom_no_batch_output, cfg, batch_size,
                        mod_count, phantom_tables, false);
                };
                measure_kernel_paths(
                    stream, start_time, stop_time, repeat_count, results.size(),
                    run_gpuntt, run_phantom, run_phantom_no_batch, gpuntt_ms,
                    phantom_ms, phantom_no_batch_ms);
            }
            else
            {
                Data64* gpuntt_data = buffer0.data();
                Data64* phantom_data = buffer1.data();
                Data64* phantom_no_batch_data = buffer2.data();

                auto run_gpuntt = [&]() {
                    if (mode == NttSurface::ForwardInplace)
                    {
                        gpuntt::GPU_NTT_Inplace(
                            gpuntt_data, device_forward_roots.data(),
                            device_moduli.data(), cfg, batch_size, mod_count);
                    }
                    else if (mode == NttSurface::InverseInplace)
                    {
                        gpuntt::GPU_INTT_Inplace(
                            gpuntt_data, device_inverse_roots.data(),
                            device_moduli.data(), cfg, batch_size, mod_count);
                    }
                    else if (mode == NttSurface::ModulusOrderedForward ||
                             mode == NttSurface::ModulusOrderedInverse)
                    {
                        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
                            gpuntt_data,
                            mode == NttSurface::ModulusOrderedForward
                                ? device_forward_roots.data()
                                : device_inverse_roots.data(),
                            device_moduli.data(), cfg, batch_size, mod_count,
                            device_order.data());
                    }
                    else
                    {
                        gpuntt::GPU_NTT_Poly_Ordered_Inplace(
                            gpuntt_data,
                            device_inverse_roots.data() +
                                (static_cast<std::size_t>(start_mod_idx)
                                 << n_power),
                            device_moduli.data() + start_mod_idx, cfg,
                            batch_size, mod_count, device_order.data());
                    }
                };
                auto run_phantom = [&]() {
                    if (mode == NttSurface::ForwardInplace)
                    {
                        heongpu::ntt::phantom_ntt_inplace(
                            phantom_data, cfg, batch_size, mod_count,
                            phantom_tables, true);
                    }
                    else if (mode == NttSurface::InverseInplace)
                    {
                        heongpu::ntt::phantom_intt_inplace(
                            phantom_data, cfg, batch_size, mod_count,
                            phantom_tables, true);
                    }
                    else if (mode == NttSurface::ModulusOrderedForward)
                    {
                        heongpu::ntt::phantom_ntt_modulus_ordered_inplace(
                            phantom_data, cfg, batch_size, mod_count,
                            device_order.data(), phantom_tables, true);
                    }
                    else if (mode == NttSurface::ModulusOrderedInverse)
                    {
                        heongpu::ntt::phantom_intt_modulus_ordered_inplace(
                            phantom_data, cfg, batch_size, mod_count,
                            device_order.data(), phantom_tables, true);
                    }
                    else
                    {
                        heongpu::ntt::phantom_intt_poly_ordered_inplace_batched(
                            phantom_data, cfg, batch_size, mod_count,
                            device_order.data(), start_mod_idx,
                            phantom_tables);
                    }
                };

                auto run_phantom_no_batch = [&]() {
                    if (mode == NttSurface::ForwardInplace)
                    {
                        heongpu::ntt::phantom_ntt_inplace(
                            phantom_no_batch_data, cfg, batch_size, mod_count,
                            phantom_tables, false);
                    }
                    else if (mode == NttSurface::InverseInplace)
                    {
                        heongpu::ntt::phantom_intt_inplace(
                            phantom_no_batch_data, cfg, batch_size, mod_count,
                            phantom_tables, false);
                    }
                    else if (mode == NttSurface::ModulusOrderedForward)
                    {
                        heongpu::ntt::phantom_ntt_modulus_ordered_inplace(
                            phantom_no_batch_data, cfg, batch_size, mod_count,
                            device_order.data(), phantom_tables, false);
                    }
                    else if (mode == NttSurface::ModulusOrderedInverse)
                    {
                        heongpu::ntt::phantom_intt_modulus_ordered_inplace(
                            phantom_no_batch_data, cfg, batch_size, mod_count,
                            device_order.data(), phantom_tables, false);
                    }
                    else
                    {
                        // Poly ordered didn't implement original phantom;
                        // Always batched due to tiny workload;
                        // So mimic one row per launch through the implemented batched primitive.
                        for (int row = 0; row < batch_size; ++row)
                        {
                            heongpu::ntt::phantom_intt_poly_ordered_inplace_batched(
                                phantom_no_batch_data, cfg, 1, 1,
                                device_order.data() + row, start_mod_idx,
                                phantom_tables);
                        }
                    }
                };

                measure_kernel_paths(
                    stream, start_time, stop_time, repeat_count, results.size(),
                    run_gpuntt, run_phantom, run_phantom_no_batch, gpuntt_ms,
                    phantom_ms, phantom_no_batch_ms);
            }

            cudaStreamDestroy(stream);

            results.push_back(
                {parameter.label, parameter.poly_modulus_degree, poly_count,
                 gpuntt_ms, phantom_ms, phantom_no_batch_ms});
        }
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);
    return results;
}

std::vector<NttKernelResult> run_bsgs_fusion_microbenchmark(
    const std::vector<ParameterSet>& parameters, int repeat_count)
{
    std::vector<NttKernelResult> results;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    for (const auto& parameter : parameters)
    {
        if (parameter.poly_modulus_degree > MAX_POLY_DEGREE)
        {
            continue;
        }

        const auto moduli = build_key_modulus(parameter);
        const int limb_count = static_cast<int>(moduli.size());
        const int n_power =
            static_cast<int>(std::log2(parameter.poly_modulus_degree));
        const std::size_t n = parameter.poly_modulus_degree;
        const std::size_t limb_elements =
            n * static_cast<std::size_t>(limb_count);
        const std::size_t ct_elements = 2 * limb_elements;
        constexpr int galois_elt = 5;
        const auto roots_base =
            heongpu::generate_primitive_root_of_unity(n, moduli);
        const auto forward_roots =
            heongpu::generate_ntt_table(roots_base, moduli, n_power);
        const auto inverse_roots =
            heongpu::generate_intt_table(roots_base, moduli, n_power);
        const auto n_inverse = heongpu::generate_n_inverse(n, moduli);
        const auto tables =
            heongpu::ntt::make_phantom_ntt_tables_from_heongpu_roots(
                moduli, forward_roots, inverse_roots, n_inverse, n_power, 0,
                true);

        const auto input =
            make_bsgs_input(n, 2 * limb_count, limb_count, moduli, 17);
        const auto addend =
            make_bsgs_input(n, limb_count, limb_count, moduli, 31);
        const auto accum =
            make_bsgs_input(n, 2 * limb_count, limb_count, moduli, 43);

        cudaStream_t stream;
        cudaStreamCreate(&stream);

        heongpu::DeviceVector<Modulus64> device_moduli(moduli);
        heongpu::DeviceVector<Data64> device_addend(addend);
        heongpu::DeviceVector<Data64> device_accum(accum);
        heongpu::DeviceVector<Data64> baby_input_unfused(input);
        heongpu::DeviceVector<Data64> baby_input_fused(input);
        heongpu::DeviceVector<Data64> giant_input_unfused(input);
        heongpu::DeviceVector<Data64> giant_input_fused(input);
        heongpu::DeviceVector<Data64> output_unfused(ct_elements);
        heongpu::DeviceVector<Data64> output_fused(ct_elements);
        heongpu::DeviceVector<Data64> scratch(ct_elements);
        cudaStreamSynchronize(stream);

        setenv("HEONGPU_USE_BSGS_FUSION", "0", 1);
        const float baby_original_ms = average_kernel_time(
            stream, start_time, stop_time, kWarmupCount, repeat_count,
            [&]() {
                heongpu::primitive::bs_add_permute_fused(
                    baby_input_unfused.data(), device_addend.data(),
                    output_unfused.data(), device_moduli.data(), galois_elt,
                    n_power, limb_count, tables, stream);
            });

        setenv("HEONGPU_USE_BSGS_FUSION", "1", 1);
        const float baby_optimized_ms = average_kernel_time(
            stream, start_time, stop_time, kWarmupCount, repeat_count,
            [&]() {
                heongpu::primitive::bs_add_permute_fused(
                    baby_input_fused.data(), device_addend.data(),
                    output_fused.data(), device_moduli.data(), galois_elt,
                    n_power, limb_count, tables, stream);
            });

        setenv("HEONGPU_USE_BSGS_FUSION", "0", 1);
        const float giant_original_ms = average_kernel_time(
            stream, start_time, stop_time, kWarmupCount, repeat_count,
            [&]() {
                heongpu::primitive::gs_add_permute_acc_fused(
                    giant_input_unfused.data(), device_addend.data(),
                    device_accum.data(), scratch.data(), output_unfused.data(),
                    device_moduli.data(), galois_elt, n_power, limb_count,
                    tables, stream);
            });

        setenv("HEONGPU_USE_BSGS_FUSION", "1", 1);
        const float giant_optimized_ms = average_kernel_time(
            stream, start_time, stop_time, kWarmupCount, repeat_count,
            [&]() {
                heongpu::primitive::gs_add_permute_acc_fused(
                    giant_input_fused.data(), device_addend.data(),
                    device_accum.data(), scratch.data(), output_fused.data(),
                    device_moduli.data(), galois_elt, n_power, limb_count,
                    tables, stream);
            });

        cudaStreamDestroy(stream);
        results.push_back({parameter.label, parameter.poly_modulus_degree,
                           0, baby_original_ms + giant_original_ms,
                           baby_optimized_ms + giant_optimized_ms, 0.0F});
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);
    return results;
}

std::vector<NttKernelResult> run_keyswitch_part2_microbenchmark(
    const std::vector<ParameterSet>& parameters,
    int repeat_count)
{
    std::vector<NttKernelResult> results;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    for (const auto& parameter : parameters)
    {
        if (parameter.poly_modulus_degree > MAX_POLY_DEGREE)
        {
            continue;
        }

        const int q_size = !parameter.q.empty()
                               ? static_cast<int>(parameter.q.size())
                               : static_cast<int>(parameter.q_bit_sizes.size());
        const int p_size = !parameter.p.empty()
                               ? static_cast<int>(parameter.p.size())
                               : static_cast<int>(parameter.p_bit_sizes.size());
        if (q_size == 0 || p_size == 0)
        {
            continue;
        }

        const int q_prime_size = q_size + p_size;
        const int n_power =
            static_cast<int>(std::log2(parameter.poly_modulus_degree));
        const std::size_t n = parameter.poly_modulus_degree;
        const std::size_t q_elements = n * static_cast<std::size_t>(q_size);
        const std::size_t output_elements = 2 * q_elements;

        const auto moduli = build_key_modulus(parameter);
        const auto half = heongpu::calculate_half(moduli, p_size);
        const auto half_mod =
            heongpu::calculate_half_mod(moduli, half, q_prime_size, p_size);
        const auto last_q_modinv =
            heongpu::calculate_last_q_modinv(moduli, q_prime_size, p_size);
        const auto roots_base =
            heongpu::generate_primitive_root_of_unity(n, moduli);
        const auto forward_roots =
            heongpu::generate_ntt_table(roots_base, moduli, n_power);
        const auto inverse_roots =
            heongpu::generate_intt_table(roots_base, moduli, n_power);
        const auto n_inverse = heongpu::generate_n_inverse(n, moduli);
        const auto tables =
            heongpu::ntt::make_phantom_ntt_tables_from_heongpu_roots(
                moduli, forward_roots, inverse_roots, n_inverse, n_power, 0,
                true);

        const auto input =
            make_bsgs_input(n, 2 * q_prime_size, q_prime_size, moduli, 71);
        const auto addend =
            make_bsgs_input(n, q_size, q_size, moduli, 97);

        cudaStream_t stream;
        cudaStreamCreate(&stream);

        heongpu::DeviceVector<Modulus64> device_moduli(moduli);
        heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);
        heongpu::DeviceVector<Data64> device_half(half);
        heongpu::DeviceVector<Data64> device_half_mod(half_mod);
        heongpu::DeviceVector<Data64> device_last_q_modinv(last_q_modinv);
        heongpu::DeviceVector<Data64> device_input(input);
        heongpu::DeviceVector<Data64> device_addend(addend);

        heongpu::DeviceVector<Data64> original_addend(q_elements);
        heongpu::DeviceVector<Data64> original_output(output_elements);
        heongpu::DeviceVector<Data64> optimized_output(output_elements);

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        const auto no_prepare = []() {};
        const auto prepare_original = [&]() {
            HEONGPU_CUDA_CHECK(cudaMemcpyAsync(
                original_addend.data(), device_addend.data(),
                q_elements * sizeof(Data64), cudaMemcpyDeviceToDevice,
                stream));
        };
        const auto run_original = [&]() {
            heongpu::primitive::keyswitch_part2_fused_moddown_ntt(
                device_input.data(), original_addend.data(),
                nullptr, original_output.data(),
                device_forward_roots.data(), device_moduli.data(), cfg_ntt,
                device_half.data(), device_half_mod.data(),
                device_last_q_modinv.data(), n_power, q_prime_size, q_size,
                q_prime_size, q_size, p_size, nullptr, stream);
        };
        const auto run_optimized = [&]() {
            heongpu::primitive::keyswitch_part2_fused_moddown_ntt(
                device_input.data(), device_addend.data(),
                nullptr, optimized_output.data(),
                device_forward_roots.data(),
                device_moduli.data(), cfg_ntt, device_half.data(),
                device_half_mod.data(), device_last_q_modinv.data(), n_power,
                q_prime_size, q_size, q_prime_size, q_size, p_size, tables,
                stream);
        };

        const float original_ms = average_prepared_kernel_time(
            stream, start_time, stop_time, kWarmupCount, repeat_count,
            prepare_original, run_original);
        const float optimized_ms = average_prepared_kernel_time(
            stream, start_time, stop_time, kWarmupCount, repeat_count,
            no_prepare, run_optimized);

        cudaStreamDestroy(stream);
        results.push_back({parameter.label, parameter.poly_modulus_degree,
                           p_size, original_ms, optimized_ms, 0.0F});
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);
    return results;
}

float speedup(float gpuntt_ms, float phantom_ms)
{
    return phantom_ms > 0.0F ? gpuntt_ms / phantom_ms : 0.0F;
}

int repeat_count_from_env(const char* name, int default_value)
{
    if (const char* value = std::getenv(name))
    {
        return std::max(1, std::atoi(value));
    }
    return default_value;
}

void print_ntt_kernel_table(std::vector<NttKernelResult> results,
                            NttSurface mode)
{
    std::sort(results.begin(), results.end(),
              [](const NttKernelResult& lhs, const NttKernelResult& rhs) {
                  if (lhs.poly_modulus_degree != rhs.poly_modulus_degree)
                  {
                      return lhs.poly_modulus_degree < rhs.poly_modulus_degree;
                  }
                  if (lhs.label != rhs.label)
                  {
                      return lhs.label < rhs.label;
                  }
                  return lhs.poly_count < rhs.poly_count;
              });

    const char* title = ntt_surface_title(mode);

    std::cout << "\n=================== Kernel microbenchmark: " << title
              << " speedup ===================" << std::endl;
    std::cout << "Speedup = GPUNTT ms / PhantomNTT ms; Interpretation: >1.00 PhantomNTT faster, <1.00 GPUNTT faster." << std::endl;

    std::vector<int> poly_counts;
    for (const auto& result : results)
    {
        if (std::find(poly_counts.begin(), poly_counts.end(),
                      result.poly_count) == poly_counts.end())
        {
            poly_counts.push_back(result.poly_count);
        }
    }
    std::sort(poly_counts.begin(), poly_counts.end());

    std::cout << std::left << std::setw(58) << "parameter / N / mode";
    for (int poly_count : poly_counts)
    {
        std::cout << std::right << std::setw(8) << poly_count;
    }
    std::cout << std::endl;
    std::cout << std::string(58 + (8 * poly_counts.size()), '-')
              << std::endl;

    for (size_t i = 0; i < results.size();)
    {
        const std::string current_label = results[i].label;
        const size_t current_n = results[i].poly_modulus_degree;
        std::vector<NttKernelResult> row_results;

        while (i < results.size() && results[i].label == current_label &&
               results[i].poly_modulus_degree == current_n)
        {
            row_results.push_back(results[i]);
            ++i;
        }

        const std::string row_prefix =
            current_label + " / N=" + std::to_string(current_n);

        auto print_row = [&](const std::string& mode, bool no_batch) {
            std::cout << std::left << std::setw(58)
                      << (row_prefix + " / " + mode);
            for (int poly_count : poly_counts)
            {
                const auto result = std::find_if(
                    row_results.begin(), row_results.end(),
                    [&](const NttKernelResult& row_result) {
                        return row_result.poly_count == poly_count;
                    });
                if (result == row_results.end())
                {
                    std::cout << std::right << std::setw(8) << "-";
                }
                else
                {
                    const float phantom_ms = no_batch
                                                 ? result->phantom_no_batch_ms
                                                 : result->phantom_ms;
                    std::cout << std::right << std::setw(8) << std::fixed
                              << std::setprecision(2)
                              << speedup(result->gpuntt_ms, phantom_ms)
                              << std::defaultfloat;
                }
            }
            std::cout << std::endl;
        };

        print_row("phantom_batch", false);
        print_row("phantom_no_batch", true);
    }

}

void print_bsgs_fusion_table(std::vector<NttKernelResult> results)
{
    std::sort(results.begin(), results.end(),
              [](const NttKernelResult& lhs, const NttKernelResult& rhs) {
                  if (lhs.poly_modulus_degree != rhs.poly_modulus_degree)
                  {
                      return lhs.poly_modulus_degree < rhs.poly_modulus_degree;
                  }
                  return lhs.label < rhs.label;
              });

    std::cout << "\n=================== BSGS fusion microbenchmark ==================="
              << std::endl;
    std::cout << "Speedup = original BSGS ms / optimized BSGS ms."
              << std::endl;
    std::cout << std::left << std::setw(30) << "parameter / N"
              << std::right << std::setw(18) << "speedup"
              << std::setw(22) << "original/opt ms"
              << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (const auto& result : results)
    {
        std::cout << std::left << std::setw(30)
                  << (result.label + " / " +
                      std::to_string(result.poly_modulus_degree))
                  << std::setw(18) << std::fixed << std::setprecision(2)
                  << speedup(result.gpuntt_ms, result.phantom_ms)
                  << std::setprecision(4)
                  << std::setw(12) << result.gpuntt_ms << "/"
                  << std::setw(9) << result.phantom_ms
                  << std::defaultfloat << std::endl;
    }
}

void print_keyswitch_optimization_table(std::vector<NttKernelResult> results)
{
    std::sort(results.begin(), results.end(),
              [](const NttKernelResult& lhs,
                 const NttKernelResult& rhs) {
                  if (lhs.poly_modulus_degree != rhs.poly_modulus_degree)
                  {
                      return lhs.poly_modulus_degree < rhs.poly_modulus_degree;
                  }
                  return lhs.label < rhs.label;
              });

    std::cout << "\n=================== KeySwitch optimization microbenchmark ==================="
              << std::endl;
    std::cout << "Baseline: original GPUNTT/no-table divide_round, "
                 "NTT(mod-down output), NTT(addend), add."
              << std::endl;
    std::cout << "KeySwitch optimization: P-specialized fused "
                 "divide_round/add_first plus PhantomNTT NTT."
              << std::endl;
    std::cout << "Speedup = original GPUNTT/no-table keyswitch ms / "
                 "optimized keyswitch ms."
              << std::endl;
    std::cout << std::left << std::setw(30) << "parameter / N"
              << std::right << std::setw(6) << "P"
              << std::right << std::setw(16) << "speedup"
              << std::endl;
    std::cout << std::string(52, '-') << std::endl;

    for (const auto& result : results)
    {
        std::cout << std::left << std::setw(30)
                  << (result.label + " / " +
                      std::to_string(result.poly_modulus_degree))
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(6) << result.poly_count
                  << std::setw(16)
                  << speedup(result.gpuntt_ms, result.phantom_ms)
                  << std::defaultfloat << std::endl;
    }

}

void print_operation_summary_table(
    const std::vector<BenchmarkResult>& gpuntt_results,
    const std::vector<BenchmarkResult>& route_results,
    const std::string& title = "HEOperator operation speedup summary")
{
    std::cout << "\n=================== " << title
              << " ===================" << std::endl;
    std::cout << "Measured HEOperator paths that call primitive::NTT/INTT wrappers:" << std::endl;
    std::cout << "  relin_II: HEArithmeticOperator::relinearize -> "
                 "HEOperator::relinearize_external_product_method2_inplace_ckks" << std::endl;
    std::cout << "  rotate_col_II: GPU rotate_col -> HEArithmeticOperator::rotate_rows -> "
                 "HEOperator::apply_galois_ckks_method_II" << std::endl;
    std::cout << "Optimized route: PhantomNTT + KeySwitch optimization + "
                 "BSGS fusion."
              << std::endl;
    std::cout << "Interpretation: GPUNTT ms / Phantom plugin, "
                 ">1.00 Phantom plugin faster, <1.00 GPUNTT faster."
              << std::endl;
    std::cout << std::left << std::setw(34) << "parameter / N"
              << std::right << std::setw(12) << "relin_II"
              << std::right << std::setw(16) << "rotate_col_II"
              << std::endl;
    std::cout << std::string(62, '-')
              << std::endl;

    for (size_t i = 0; i < gpuntt_results.size(); ++i)
    {
        const auto& g = gpuntt_results[i];
        const auto& p = route_results[i];
        std::cout << std::left << std::setw(34)
                  << (g.label + " / " +
                      std::to_string(g.poly_modulus_degree))
                  << std::right << std::setw(14)
                  << std::fixed << std::setprecision(2)
                  << speedup(g.timing.relinearize, p.timing.relinearize)
                  << std::right << std::setw(12)
                  << speedup(g.timing.rotate, p.timing.rotate)
                  << std::defaultfloat << std::endl;
    }
}

void print_bootstrap_summary_table(
    const std::vector<BootstrapRouteResult>& results)
{
    std::cout << "\n=================== regular_bootstrapping_v2 bootstrap summary ===================" << std::endl;
    std::cout << "PhantomNTT: CKKS NTT/INTT primitive replacement." << std::endl;
    std::cout << "KeySwitch optimization: mod-major base conversion, "
                 "excluded-limb NTT, key multiplication, and fused "
                 "divide_round/add_first before the final NTT."
              << std::endl;
    std::cout << "BSGS fusion: fuse add, NTT-domain permutation, and "
                 "accumulation inside BSGS matrix path."
              << std::endl;
    if (results.empty())
    {
        return;
    }

    const auto& baseline = results.front().result;
    std::cout << "Parameter: " << baseline.label << " / N="
              << baseline.poly_modulus_degree << std::endl;
    std::cout << std::left << std::setw(26) << "route"
              << std::right << std::setw(16) << "absolute ms"
              << std::setw(18) << "speedup"
              << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (const auto& result : results)
    {
        std::cout << std::left << std::setw(26) << result.name
                  << std::right << std::fixed << std::setprecision(3)
                  << std::setw(16) << result.result.timing.bootstrapping
                  << std::setprecision(2)
                  << std::setw(18)
                  << speedup(baseline.timing.bootstrapping,
                             result.result.timing.bootstrapping)
                  << std::defaultfloat << std::endl;
    }
}

int main()
{
    const int repeat_count =
        repeat_count_from_env("REPEAT",
                              kDefaultRepeatCount);
    const int bootstrap_repeat_count =
        repeat_count_from_env("REPEAT_BOOT",
                              kDefaultBootstrapRepeatCount);

    auto parameters = make_parameters();
    sort_parameters_by_n(parameters);

    std::vector<ParameterSet> operation_parameters;
    std::copy_if(parameters.begin(), parameters.end(),
                 std::back_inserter(operation_parameters),
                 [](const ParameterSet& parameter) {
                     return parameter.poly_modulus_degree <= MAX_POLY_DEGREE;
                 });

    std::cout << "CKKS NTT backend benchmark" << std::endl;
    std::cout << "Repeat count: " << repeat_count << std::endl;
    std::cout << "Bootstrap repeat count: " << bootstrap_repeat_count
              << std::endl;

    // Kernel microbenchmarks.
    std::cout << "\nOptimization: GPUNTT vs PhantomNTT kernel - Inplace, Modulus Ordered and Poly Ordered" << std::endl;

    const NttSurface kernel_modes[] = {
        NttSurface::ForwardInplace,
        NttSurface::InverseInplace,
        NttSurface::InverseOutOfPlace,
        NttSurface::ModulusOrderedForward,
        NttSurface::ModulusOrderedInverse,
        NttSurface::PolyOrderedInverse};

    for (NttSurface mode : kernel_modes)
    {
        const auto results =
            run_ntt_kernel_comparison(parameters, repeat_count, mode);
        print_ntt_kernel_table(results, mode);
    }

    const auto bsgs_results =
        run_bsgs_fusion_microbenchmark(parameters, repeat_count);
    print_bsgs_fusion_table(bsgs_results);

    const auto keyswitch_part2_results =
        run_keyswitch_part2_microbenchmark(parameters, repeat_count);
    print_keyswitch_optimization_table(keyswitch_part2_results);

    // HEOperator benchmarks.
    std::vector<BenchmarkResult> gpuntt_results;
    std::vector<BenchmarkResult> phantom_results;
    gpuntt_results.reserve(operation_parameters.size());
    phantom_results.reserve(operation_parameters.size());

    for (size_t i = 0; i < operation_parameters.size(); ++i)
    {
        const auto& parameter = operation_parameters[i];
        BenchmarkResult gpuntt;
        BenchmarkResult phantom;

        if ((i % 2) == 0)
        {
            gpuntt = run_operator_case(parameter, true);
            phantom = run_operator_case(parameter, false);
        }
        else
        {
            phantom = run_operator_case(parameter, false);
            gpuntt = run_operator_case(parameter, true);
        }

        gpuntt_results.push_back(gpuntt);
        phantom_results.push_back(phantom);
    }

    print_operation_summary_table(gpuntt_results, phantom_results,
                                  "HEOperator operation speedup summary");

    // Bootstrapping benchmark.
    const auto boot_parameter =
        std::find_if(operation_parameters.begin(), operation_parameters.end(),
                     [](const ParameterSet& parameter) {
                         return parameter.label == "N16QP1546H192H32";
                     });
    if (boot_parameter != operation_parameters.end())
    {
        const std::vector<OptimizationRoute> bootstrap_routes = {
            {"Original", false, false, false, false},
            {"KeySwitch optimization", true, true, true, false},
            {"Full optimized", true, true, true, true}};

        std::vector<BootstrapRouteResult> bootstrap_results;
        bootstrap_results.reserve(bootstrap_routes.size());
        for (const auto& route : bootstrap_routes)
        {
            bootstrap_results.push_back(
                {route.name,
                 run_bootstrapping_case(*boot_parameter, route,
                                        bootstrap_repeat_count)});
        }

        print_bootstrap_summary_table(bootstrap_results);
    }
    return EXIT_SUCCESS;
}
