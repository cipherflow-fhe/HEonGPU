// Copyright 2024-2026 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include <heongpu/heongpu.hpp>
#include <heongpu/ntt/ntt.cuh>
#include <heongpu/primitive/ntt.cuh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

constexpr auto Scheme = heongpu::Scheme::CKKS;
constexpr int kKernelPolyCounts[] = {1, 2, 4, 8, 16, 32};
constexpr int kKernelWarmupCount = 2;

struct ParameterSet
{
    std::string label;
    size_t poly_modulus_degree;
    std::vector<int> q_bit_sizes;
    std::vector<int> p_bit_sizes;
    std::vector<Data64> q;
    std::vector<Data64> p;
    double scale;
};

struct OperationTiming
{
    float relinearize = 0.0F;
    float rotate = 0.0F;
    float rescale = 0.0F;
    float bootstrapping = 0.0F;
};

struct BenchmarkResult
{
    std::string label;
    size_t poly_modulus_degree = 0;
    OperationTiming timing;
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

enum class NttKernelMode
{
    ForwardInplace,
    InverseInplace,
    InverseOutOfPlace
};

std::vector<ParameterSet> make_parameters()
{
    return {
        {.label = "PN12QP109",
         .poly_modulus_degree = 4096,
         .q = {0x200000e001ULL, 0x100006001ULL},
         .p = {0x3ffffea001ULL},
         .scale = std::ldexp(1.0, 32)},
        {.label = "PN13QP218",
         .poly_modulus_degree = 8192,
         .q = {0x1fffec001ULL, 0x3fff4001ULL, 0x3ffe8001ULL,
               0x40020001ULL, 0x40038001ULL, 0x3ffc0001ULL},
         .p = {0x800004001ULL},
         .scale = std::ldexp(1.0, 30)},
        {.label = "PN14QP438",
         .poly_modulus_degree = 16384,
         .q = {0x200000008001ULL, 0x400018001ULL, 0x3fffd0001ULL,
               0x400060001ULL, 0x400068001ULL, 0x3fff90001ULL,
               0x400080001ULL, 0x4000a8001ULL, 0x400108001ULL,
               0x3ffeb8001ULL},
         .p = {0x7fffffd8001ULL, 0x7fffffc8001ULL},
         .scale = std::ldexp(1.0, 34)},
        {.label = "PN15QP880",
         .poly_modulus_degree = 32768,
         .q = {0x4000000120001ULL, 0x10000140001ULL, 0xffffe80001ULL,
               0x10000290001ULL, 0xffffc40001ULL, 0x100003e0001ULL,
               0x10000470001ULL, 0x100004b0001ULL, 0xffffb20001ULL,
               0x10000500001ULL, 0x10000650001ULL, 0xffff940001ULL,
               0xffff8a0001ULL, 0xffff820001ULL, 0xffff780001ULL,
               0x10000890001ULL, 0xffff750001ULL, 0x10000960001ULL},
         .p = {0x40000001b0001ULL, 0x3ffffffdf0001ULL,
               0x4000000270001ULL},
         .scale = std::ldexp(1.0, 40)},
        {.label = "PN16QP1761",
         .poly_modulus_degree = 65536,
         .q = {0x80000000080001ULL, 0x2000000a0001ULL,
               0x2000000e0001ULL, 0x1fffffc20001ULL, 0x200000440001ULL,
               0x200000500001ULL, 0x200000620001ULL, 0x1fffff980001ULL,
               0x2000006a0001ULL, 0x1fffff7e0001ULL, 0x200000860001ULL,
               0x200000a60001ULL, 0x200000aa0001ULL, 0x200000b20001ULL,
               0x200000c80001ULL, 0x1fffff360001ULL, 0x200000e20001ULL,
               0x1fffff060001ULL, 0x200000fe0001ULL, 0x1ffffede0001ULL,
               0x1ffffeca0001ULL, 0x1ffffeb40001ULL, 0x200001520001ULL,
               0x1ffffe760001ULL, 0x2000019a0001ULL, 0x1ffffe640001ULL,
               0x200001a00001ULL, 0x1ffffe520001ULL, 0x200001e80001ULL,
               0x1ffffe0c0001ULL, 0x1ffffdee0001ULL, 0x200002480001ULL,
               0x1ffffdb60001ULL, 0x200002560001ULL},
         .p = {0x80000000440001ULL, 0x7fffffffba0001ULL,
               0x80000000500001ULL, 0x7fffffffaa0001ULL},
         .scale = std::ldexp(1.0, 45)},
        {.label = "N16QP1546H192H32",
         .poly_modulus_degree = 65536,
         .q = {0x10000000006e0001ULL, 0x10000140001ULL,
               0xffffe80001ULL, 0xffffc40001ULL, 0x100003e0001ULL,
               0xffffb20001ULL, 0x10000500001ULL, 0xffff940001ULL,
               0xffff8a0001ULL, 0xffff820001ULL, 0x7fffe60001ULL,
               0x7fffe40001ULL, 0x7fffe00001ULL, 0xfffffffff840001ULL,
               0x1000000000860001ULL, 0xfffffffff6a0001ULL,
               0x1000000000980001ULL, 0xfffffffff5a0001ULL,
               0x1000000000b00001ULL, 0x1000000000ce0001ULL,
               0xfffffffff2a0001ULL, 0x100000000060001ULL,
               0xfffffffff00001ULL, 0xffffffffd80001ULL,
               0x1000000002a0001ULL},
         .p = {0x1fffffffffe00001ULL, 0x1fffffffffc80001ULL,
               0x1fffffffffb40001ULL, 0x1fffffffff500001ULL,
               0x1fffffffff420001ULL},
         .scale = std::ldexp(1.0, 40)},
        {.label = "PhantomN17QP3570",
         .poly_modulus_degree = 131072,
         .q_bit_sizes = [] {
             std::vector<int> bits(54, 55);
             std::fill(bits.begin(), bits.begin() + 5, 60);
             return bits;
         }(),
         .p_bit_sizes = [] {
             std::vector<int> bits(10, 55);
             std::fill(bits.begin(), bits.begin() + 5, 60);
             return bits;
         }(),
         .scale = std::ldexp(1.0, 50)}};
}

std::vector<Modulus64> build_key_modulus(const ParameterSet& parameter)
{
    if (!parameter.q_bit_sizes.empty() || !parameter.p_bit_sizes.empty())
    {
        std::vector<int> bit_sizes = parameter.q_bit_sizes;
        bit_sizes.insert(bit_sizes.end(), parameter.p_bit_sizes.begin(),
                         parameter.p_bit_sizes.end());
        return heongpu::generate_primes(parameter.poly_modulus_degree,
                                        bit_sizes);
    }

    std::vector<Modulus64> moduli;
    moduli.reserve(parameter.q.size() + parameter.p.size());
    for (Data64 value : parameter.q)
        moduli.emplace_back(value);

    for (Data64 value : parameter.p)
        moduli.emplace_back(value);

    return moduli;
}

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
    for (int trial = 0; trial < warmup_count; ++trial){
        func();
    }
    cudaStreamSynchronize(stream);

    float total = 0.0F;
    for (int trial = 0; trial < repeat_count; ++trial){
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
                                  bool use_phantom_ntt, int repeat_count)
{
    setenv("HEONGPU_USE_PHANTOM_NTT", use_phantom_ntt ? "1" : "0", 1);

    heongpu::HEContext<Scheme> context = heongpu::GenHEContext<Scheme>(heongpu::sec_level_type::none);
    context->set_poly_modulus_degree(parameter.poly_modulus_degree);
    context->set_coeff_modulus_values(parameter.q, parameter.p);
    context->generate();

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

    for (int trial = 0; trial < repeat_count; trial++)
    {
        heongpu::Plaintext<Scheme> plaintext(context);
        encoder.encode(plaintext, message, parameter.scale);

        heongpu::Ciphertext<Scheme> c1(context);
        encryptor.encrypt(c1, plaintext);

        {
            heongpu::Ciphertext<Scheme> operand(context);
            operators.multiply(c1, c1, operand);
            time_call(start_time, stop_time, timing.relinearize, [&]() {
                operators.relinearize_inplace(operand, relin_key);
            });
        }

        {
            heongpu::Ciphertext<Scheme> operand(context);
            operators.multiply(c1, c1, operand);
            operators.relinearize_inplace(operand, relin_key);
            time_call(start_time, stop_time, timing.rescale,
                      [&]() { operators.rescale_inplace(operand); });
        }

        {
            heongpu::Ciphertext<Scheme> operand = c1;
            time_call(start_time, stop_time, timing.rotate, [&]() {
                operators.rotate_rows(operand, operand, galois_key, 1);
            });
        }

        cudaDeviceSynchronize();
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);

    timing.relinearize /= repeat_count;
    timing.rotate /= repeat_count;
    timing.rescale /= repeat_count;

    return {parameter.label, parameter.poly_modulus_degree, timing};
}

BenchmarkResult run_bootstrapping_case(const ParameterSet& parameter,
                                       bool use_phantom_ntt, int repeat_count)
{
    setenv("HEONGPU_USE_PHANTOM_NTT", use_phantom_ntt ? "1" : "0", 1);

    heongpu::HEContext<Scheme> context =
        heongpu::GenHEContext<Scheme>(heongpu::sec_level_type::none);
    context->set_poly_modulus_degree(parameter.poly_modulus_degree);
    context->set_slot_count(
        static_cast<int>(parameter.poly_modulus_degree / 2));
    context->set_coeff_modulus_values(parameter.q, parameter.p);
    context->generate();

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

    heongpu::EvalModConfig eval_mod_config(
        context->get_key_modulus()[0].value, 20, 256.0, 16, 30, 3, 0,
        std::ldexp(1.0, 60));

    heongpu::BootstrappingConfigV2 boot_config(
        heongpu::EncodingMatrixConfig(
            heongpu::LinearTransformType::SLOTS_TO_COEFFS, 12, 2.0, 3),
        eval_mod_config,
        heongpu::EncodingMatrixConfig(
            heongpu::LinearTransformType::COEFFS_TO_SLOTS, 24, 2.0, 4));

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

std::vector<NttKernelResult> run_ntt_kernel_comparison(
    const std::vector<ParameterSet>& parameters, int repeat_count,
    NttKernelMode mode)
{
    std::vector<NttKernelResult> results;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    auto measure_paths = [&](cudaStream_t stream, auto&& run_gpuntt,
                             auto&& run_phantom,
                             auto&& run_phantom_no_batch,
                             std::size_t order_seed,
                             float& gpuntt_ms, float& phantom_ms,
                             float& phantom_no_batch_ms) {
        auto measure_gpuntt = [&]() {
            return average_kernel_time(
                stream, start_time, stop_time, kKernelWarmupCount,
                repeat_count, run_gpuntt);
        };
        auto measure_phantom = [&]() {
            return average_kernel_time(
                stream, start_time, stop_time, kKernelWarmupCount,
                repeat_count, run_phantom);
        };
        auto measure_phantom_no_batch = [&]() {
            return average_kernel_time(
                stream, start_time, stop_time, kKernelWarmupCount,
                repeat_count, run_phantom_no_batch);
        };

        enum { GpunttPath, PhantomPath, PhantomNoBatchPath };
        const int three_path_orders[] = {
            GpunttPath, PhantomPath, PhantomNoBatchPath,
            PhantomPath, PhantomNoBatchPath, GpunttPath,
            PhantomNoBatchPath, GpunttPath, PhantomPath};
        const int* path_order = &three_path_orders[(order_seed % 3) * 3];

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
    };

    for (const auto& parameter : parameters)
    {
        if (parameter.poly_modulus_degree < 8192 ||
            parameter.poly_modulus_degree > 131072){
            continue;
        }

        const auto moduli = build_key_modulus(parameter);
        const int mod_count = static_cast<int>(moduli.size());
        const int n_power =
            static_cast<int>(std::log2(parameter.poly_modulus_degree));
        const auto roots_base = heongpu::generate_primitive_root_of_unity(
            parameter.poly_modulus_degree, moduli);
        const auto forward_roots = heongpu::generate_ntt_table(roots_base, moduli, n_power);
        const auto inverse_roots = heongpu::generate_intt_table(roots_base, moduli, n_power);
        const auto n_inverse = heongpu::generate_n_inverse(parameter.poly_modulus_degree, moduli);

        heongpu::DeviceVector<Modulus64> device_moduli(moduli);
        heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);
        heongpu::DeviceVector<Root64> device_inverse_roots(inverse_roots);
        heongpu::DeviceVector<Ninverse64> device_n_inverse(n_inverse);
        auto phantom_tables =
            heongpu::primitive::make_phantom_ntt_tables_from_heongpu_roots(
                moduli, forward_roots, inverse_roots, n_inverse, n_power, cudaStreamLegacy);

        for (int poly_count : kKernelPolyCounts)
        {
            const int batch_size = mod_count * poly_count;
            const std::size_t element_count = parameter.poly_modulus_degree * static_cast<std::size_t>(batch_size);
            std::size_t free_memory = 0;
            std::size_t total_memory = 0;
            cudaMemGetInfo(&free_memory, &total_memory);
            const std::size_t buffer_count = 4;
            if ((buffer_count * element_count * sizeof(Data64)) > static_cast<std::size_t>(free_memory * 0.70)){
                continue;
            }

            cudaStream_t stream;
            cudaStreamCreate(&stream);

            gpuntt::ntt_rns_configuration<Data64> cfg = {
                .n_power = n_power,
                .ntt_type = mode == NttKernelMode::ForwardInplace
                                ? gpuntt::FORWARD
                                : gpuntt::INVERSE,
                .ntt_layout = gpuntt::PerPolynomial,
                .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
                .zero_padding = false,
                .mod_inverse = device_n_inverse.data(),
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

            if (mode == NttKernelMode::InverseOutOfPlace)
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
                    heongpu::ntt::phantom_intt_batched(
                        input_data, phantom_output, cfg, batch_size, mod_count,
                        phantom_tables);
                };

                auto run_phantom_no_batch = [&]() {
                    heongpu::ntt::phantom_intt(
                        input_data, phantom_no_batch_output, cfg, batch_size,
                        mod_count, phantom_tables);
                };
                measure_paths(
                    stream, run_gpuntt, run_phantom, run_phantom_no_batch,
                    results.size(), gpuntt_ms, phantom_ms, phantom_no_batch_ms);
            }
            else
            {
                Data64* gpuntt_data = buffer0.data();
                Data64* phantom_data = buffer1.data();
                Data64* phantom_no_batch_data = buffer2.data();

                auto run_gpuntt = [&]() {
                    if (mode == NttKernelMode::InverseInplace)
                    {
                        gpuntt::GPU_INTT_Inplace(
                            gpuntt_data, device_inverse_roots.data(),
                            device_moduli.data(), cfg, batch_size, mod_count);
                    }
                    else{
                        gpuntt::GPU_NTT_Inplace(
                            gpuntt_data, device_forward_roots.data(),
                            device_moduli.data(), cfg, batch_size, mod_count);
                    }
                };
                auto run_phantom = [&]() {
                    if (mode == NttKernelMode::InverseInplace)
                    {
                        heongpu::ntt::phantom_intt_inplace_batched(
                            phantom_data, cfg, batch_size, mod_count,
                            phantom_tables);
                    }
                    else{
                        heongpu::ntt::phantom_ntt_inplace_batched(
                            phantom_data, cfg, batch_size, mod_count,
                            phantom_tables);
                    }
                };

                auto run_phantom_no_batch = [&]() {
                    if (mode == NttKernelMode::InverseInplace)
                    {
                        heongpu::ntt::phantom_intt_inplace(
                            phantom_no_batch_data, cfg, batch_size, mod_count,
                            phantom_tables);
                    }
                    else
                    {
                        heongpu::ntt::phantom_ntt_inplace(
                            phantom_no_batch_data, cfg, batch_size, mod_count,
                            phantom_tables);
                    }
                };

                measure_paths(
                    stream, run_gpuntt, run_phantom, run_phantom_no_batch,
                    results.size(), gpuntt_ms, phantom_ms, phantom_no_batch_ms);
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

float speedup(float gpuntt_ms, float phantom_ms){
    return phantom_ms > 0.0F ? gpuntt_ms / phantom_ms : 0.0F;
}

void sort_parameters_by_n(std::vector<ParameterSet>& parameters)
{
    std::sort(parameters.begin(), parameters.end(),
              [](const ParameterSet& lhs, const ParameterSet& rhs) {
                  if (lhs.poly_modulus_degree != rhs.poly_modulus_degree)
                  {
                      return lhs.poly_modulus_degree < rhs.poly_modulus_degree;
                  }
                  return lhs.label < rhs.label;
              });
}

void print_ntt_kernel_table(std::vector<NttKernelResult> results,
                            const std::string& transform)
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

    std::cout << "\n=================== Kernel microbenchmark: " << transform << " NTT speedup ===================" << std::endl;
    std::cout << "Comparison: GPUNTT RNS " << transform << " vs HEonGPU PhantomNTT " << transform << std::endl;
    std::cout << "Speedup = GPUNTT ms / PhantomNTT ms; Interpretation: >1.00 PhantomNTT faster, <1.00 GPUNTT faster." << std::endl;
    std::cout << "Modes: phantom_batch - batch poly operation; phantom_no_batch - one poly at a time." << std::endl;

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

void print_operation_summary_table(
    const std::vector<BenchmarkResult>& gpuntt_results,
    const std::vector<BenchmarkResult>& phantom_results)
{
    std::cout << "\n=================== HEOperator operation speedup summary ==================="<< std::endl;
    std::cout << "Measured HEOperator paths that call primitive::NTT/INTT wrappers:" << std::endl;
    std::cout << "  relin_II: HEArithmeticOperator::relinearize -> "
                 "HEOperator::relinearize_external_product_method2_inplace_ckks" << std::endl;
    std::cout << "  rotate_col_II: GPU rotate_col -> HEArithmeticOperator::rotate_rows -> "
                 "HEOperator::apply_galois_ckks_method_II" << std::endl;
    std::cout << "  rescale_leveled: HEArithmeticOperator::rescale -> HEOperator::rescale_inplace_ckks_leveled" << std::endl;
    std::cout << "  regular_bootstrapping_v2: HEArithmeticOperator::regular_bootstrapping_v2" << std::endl;
    std::cout << "Values are GPUNTT ms / PhantomNTT ms." << std::endl;
    std::cout << "Interpretation: >1.00 PhantomNTT faster, <1.00 GPUNTT faster." << std::endl;
    std::cout << std::left << std::setw(34) << "parameter / N"
              << std::right << std::setw(12) << "relin_II"
              << std::right << std::setw(16) << "rotate_col_II"
              << std::right << std::setw(16) << "rescale_leveled"
              << std::right << std::setw(24) << "regular_bootstrap_v2"
              << std::endl;
    std::cout << std::string(102, '-')
              << std::endl;

    for (size_t i = 0; i < gpuntt_results.size(); ++i)
    {
        const auto& g = gpuntt_results[i];
        const auto& p = phantom_results[i];
        std::cout << std::left << std::setw(34)
                  << (g.label + " / " +
                      std::to_string(g.poly_modulus_degree))
                  << std::right << std::setw(14)
                  << std::fixed << std::setprecision(2)
                  << speedup(g.timing.relinearize, p.timing.relinearize)
                  << std::right << std::setw(12)
                  << speedup(g.timing.rotate, p.timing.rotate)
                  << std::right << std::setw(10)
                  << speedup(g.timing.rescale, p.timing.rescale)
                  << std::right << std::setw(24);
        if (g.timing.bootstrapping > 0.0F && p.timing.bootstrapping > 0.0F)
        {
            std::cout << speedup(g.timing.bootstrapping,
                                 p.timing.bootstrapping);
        }
        std::cout << std::defaultfloat << std::endl;
    }
}

void run_phantomntt_per_poly_ntt_benchmark(
    const std::string& optimization_name, int repeat_count)
{
    auto parameters = make_parameters();
    sort_parameters_by_n(parameters);

    std::vector<ParameterSet> operation_parameters;
    std::copy_if(parameters.begin(), parameters.end(),
                 std::back_inserter(operation_parameters),
                 [](const ParameterSet& parameter) {
                     return parameter.poly_modulus_degree <= MAX_POLY_DEGREE;
                 });

    std::cout << "\nOptimization: " << optimization_name << std::endl;

    const auto forward_ntt_results =
        run_ntt_kernel_comparison(parameters, repeat_count, NttKernelMode::ForwardInplace);
    print_ntt_kernel_table(forward_ntt_results, "forward_inplace");

    const auto inverse_ntt_results =
        run_ntt_kernel_comparison(parameters, repeat_count, NttKernelMode::InverseInplace);
    print_ntt_kernel_table(inverse_ntt_results, "inverse_inplace");

    const auto inverse_outofplace_results =
        run_ntt_kernel_comparison(parameters, repeat_count, NttKernelMode::InverseOutOfPlace);
    print_ntt_kernel_table(inverse_outofplace_results, "inverse_outofplace");

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
            gpuntt = run_operator_case(parameter, false, repeat_count);
            phantom = run_operator_case(parameter, true, repeat_count);
        }
        else
        {
            phantom = run_operator_case(parameter, true, repeat_count);
            gpuntt = run_operator_case(parameter, false, repeat_count);
        }

        gpuntt_results.push_back(gpuntt);
        phantom_results.push_back(phantom);
    }

    const auto boot_parameter =
        std::find_if(operation_parameters.begin(), operation_parameters.end(),
                     [](const ParameterSet& parameter) {
                         return parameter.label == "N16QP1546H192H32";
                     });
    if (boot_parameter != operation_parameters.end())
    {
        BenchmarkResult gpuntt_boot =
            run_bootstrapping_case(*boot_parameter, false, repeat_count);
        BenchmarkResult phantom_boot =
            run_bootstrapping_case(*boot_parameter, true, repeat_count);

        for (size_t i = 0; i < gpuntt_results.size(); ++i)
        {
            if (gpuntt_results[i].label == boot_parameter->label)
            {
                gpuntt_results[i].timing.bootstrapping =
                    gpuntt_boot.timing.bootstrapping;
                phantom_results[i].timing.bootstrapping =
                    phantom_boot.timing.bootstrapping;
                break;
            }
        }
    }

    print_operation_summary_table(gpuntt_results, phantom_results);
}

int main()
{
    int repeat_count = 10;
    if (const char* repeat_env = std::getenv("HEONGPU_CKKS_NTT_BENCH_REPEAT"))
    {
        repeat_count = std::max(1, std::atoi(repeat_env));
    }

    std::cout << "CKKS NTT backend benchmark" << std::endl;
    std::cout << "Repeat count: " << repeat_count << std::endl;

    run_phantomntt_per_poly_ntt_benchmark(
        "PhantomNTT inplace NTT - batch switch", repeat_count);
    return EXIT_SUCCESS;
}
