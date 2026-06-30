// Copyright 2024-2026 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include <heongpu/heongpu.hpp>
#include <heongpu/primitive/ntt.cuh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
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
    std::vector<int> log_q_bits;
    std::vector<int> log_p_bits;
    std::vector<Data64> q_values;
    std::vector<Data64> p_values;
    double scale;
};

struct OperationTiming
{
    float relinearize = 0.0F;
    float rotate = 0.0F;
    float rescale = 0.0F;
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
    std::string transform;
    size_t poly_modulus_degree = 0;
    int poly_count = 0;
    int mod_count = 0;
    int batch_size = 0;
    float gpuntt_ms = 0.0F;
    float phantom_ms = 0.0F;
    float phantom_no_batch_ms = 0.0F;
};

std::vector<ParameterSet> make_parameters()
{
    return {
        {"PN12QP109", 4096, {}, {},
         {0x200000e001ULL, 0x100006001ULL}, {0x3ffffea001ULL},
         pow(2.0, 32)},
        {"PN13QP218", 8192, {}, {},
         {0x1fffec001ULL, 0x3fff4001ULL, 0x3ffe8001ULL, 0x40020001ULL,
          0x40038001ULL, 0x3ffc0001ULL},
         {0x800004001ULL}, pow(2.0, 30)},
        {"TestPN14QP438", 16384, {}, {},
         {0x100000000060001ULL, 0x80000000068001ULL, 0x80000000080001ULL,
          0x3fffffffef8001ULL, 0x40000000120001ULL,
          0x3fffffffeb8001ULL},
         {0x80000000130001ULL, 0x7fffffffe90001ULL}, pow(2.0, 34)},
        {"TestPN15QP880", 32768, {}, {},
         {0x7ffffffffe70001ULL, 0x7ffffffffe10001ULL,
          0x7ffffffffcc0001ULL, 0x400000000270001ULL,
          0x400000000350001ULL, 0x400000000360001ULL,
          0x3ffffffffc10001ULL, 0x3ffffffffbe0001ULL,
          0x3ffffffffbd0001ULL, 0x4000000004d0001ULL,
          0x400000000570001ULL, 0x400000000660001ULL},
         {0xffffffffffc0001ULL, 0x10000000001d0001ULL,
          0x10000000006e0001ULL},
         pow(2.0, 40)},
        {"TestPN16QP240", 65536, {60, 60, 60}, {60}, {}, {},
         pow(2.0, 40)},
        {"TestPN17QP360", 131072, {60, 60, 60, 60}, {60, 60}, {},
         {}, pow(2.0, 40)}};
}

void configure_context(heongpu::HEContext<Scheme>& context,
                       const ParameterSet& parameter)
{
    context->set_poly_modulus_degree(parameter.poly_modulus_degree);
    if (!parameter.q_values.empty()){
        context->set_coeff_modulus_values(parameter.q_values, parameter.p_values);
    }
    else{
        context->set_coeff_modulus_bit_sizes(parameter.log_q_bits, parameter.log_p_bits);
    }
    context->generate();
}

std::vector<Modulus64> build_key_modulus(const ParameterSet& parameter)
{
    if (!parameter.q_values.empty()){
        std::vector<Modulus64> moduli;
        moduli.reserve(parameter.q_values.size() + parameter.p_values.size());
        for (Data64 value : parameter.q_values)
            moduli.emplace_back(value);
        
        for (Data64 value : parameter.p_values)
            moduli.emplace_back(value);
        
        return moduli;
    }

    std::vector<int> bit_sizes = parameter.log_q_bits;
    bit_sizes.insert(bit_sizes.end(), parameter.log_p_bits.begin(), parameter.log_p_bits.end());
    return heongpu::generate_primes(parameter.poly_modulus_degree, bit_sizes);
}

template <typename Func>
void time_call(cudaEvent_t start_time, cudaEvent_t stop_time, float& total,
               bool record, Func&& func)
{
    cudaEventRecord(start_time);
    func();
    cudaEventRecord(stop_time);

    cudaEventSynchronize(stop_time);
    if (!record){
        return;
    }

    float elapsed = 0.0F;
    cudaEventElapsedTime(&elapsed, start_time, stop_time);
    total += elapsed;
}

template <typename Func>
float average_cuda_time(cudaEvent_t start_time, cudaEvent_t stop_time,
                        int warmup_count, int repeat_count, Func&& func)
{
    for (int trial = 0; trial < warmup_count; ++trial){
        func();
    }
    cudaDeviceSynchronize();

    float total = 0.0F;
    for (int trial = 0; trial < repeat_count; trial++){
        time_call(start_time, stop_time, total, true, func);
    }
    return total / repeat_count;
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

BenchmarkResult run_case(const ParameterSet& parameter, bool use_phantom_ntt,
                         int repeat_count)
{
    setenv("HEONGPU_USE_PHANTOM_NTT", use_phantom_ntt ? "1" : "0", 1);

    heongpu::HEContext<Scheme> context = heongpu::GenHEContext<Scheme>(heongpu::sec_level_type::none);
    configure_context(context, parameter);

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
    heongpu::HEDecryptor<Scheme> decryptor(context, secret_key);
    heongpu::HEArithmeticOperator<Scheme> operators(context, encoder);

    const int row_size = parameter.poly_modulus_degree / 2;
    heongpu::HostVector<double> message(row_size, 1);

    OperationTiming timing;
    float unreported_timing = 0.0F;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

    for (int trial = 0; trial < repeat_count; trial++)
    {
        const bool record = trial >= 0;
        heongpu::Plaintext<Scheme> plaintext(context);

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            encoder.encode(plaintext, message, parameter.scale);
        });

        heongpu::Ciphertext<Scheme> c1(context);
        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            encryptor.encrypt(c1, plaintext);
        });

        heongpu::Ciphertext<Scheme> c2(context);
        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            operators.add(c1, c1, c2);
        });

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            operators.sub(c2, c1, c2);
        });

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            operators.multiply(c2, c1, c2);
        });

        time_call(start_time, stop_time, timing.relinearize, record, [&]() {
            operators.relinearize_inplace(c2, relin_key);
        });

        time_call(start_time, stop_time, timing.rescale, record, [&]() {
            operators.rescale_inplace(c2);
        });

        heongpu::Ciphertext<Scheme> c3(context);
        encryptor.encrypt(c3, plaintext);

        time_call(start_time, stop_time, timing.rotate, record, [&]() {
            operators.rotate_rows(c3, c3, galois_key, 1);
        });

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            operators.add_plain_inplace(c3, plaintext);
        });

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            operators.sub_plain_inplace(c3, plaintext);
        });

        heongpu::Ciphertext<Scheme> c4(context);
        encryptor.encrypt(c4, plaintext);

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            operators.multiply_plain(c4, plaintext, c4);
        });

        heongpu::Plaintext<Scheme> decrypted(context);

        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            decryptor.decrypt(decrypted, c3);
        });

        heongpu::HostVector<double> decoded;
        time_call(start_time, stop_time, unreported_timing, record, [&]() {
            encoder.decode(decoded, decrypted);
        });

        cudaDeviceSynchronize();
    }

    cudaEventDestroy(start_time);
    cudaEventDestroy(stop_time);

    timing.relinearize /= repeat_count;
    timing.rotate /= repeat_count;
    timing.rescale /= repeat_count;

    return {parameter.label, parameter.poly_modulus_degree, timing};
}

std::vector<NttKernelResult> run_ntt_kernel_comparison(
    const std::vector<ParameterSet>& parameters, int repeat_count,
    bool inverse)
{
    std::vector<NttKernelResult> results;
    cudaEvent_t start_time, stop_time;
    cudaEventCreate(&start_time);
    cudaEventCreate(&stop_time);

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
        const auto roots_base = heongpu::generate_primitive_root_of_unity(parameter.poly_modulus_degree, moduli);
        const auto forward_roots = heongpu::generate_ntt_table(roots_base, moduli, n_power);
        const auto inverse_roots = heongpu::generate_intt_table(roots_base, moduli, n_power);
        const auto n_inverse = heongpu::generate_n_inverse(parameter.poly_modulus_degree, moduli);

        heongpu::DeviceVector<Modulus64> device_moduli(moduli);
        heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);
        heongpu::DeviceVector<Root64> device_inverse_roots(inverse_roots);
        heongpu::DeviceVector<Ninverse64> device_n_inverse(n_inverse);
        auto phantom_tables =
            heongpu::primitive::make_phantom_ntt_tables_from_heongpu_roots(
                moduli, forward_roots, inverse_roots, n_inverse, n_power, 0);

        for (int poly_count : kKernelPolyCounts)
        {
            const int batch_size = mod_count * poly_count;
            const std::size_t element_count = parameter.poly_modulus_degree * static_cast<std::size_t>(batch_size);
            std::size_t free_memory = 0;
            std::size_t total_memory = 0;
            cudaMemGetInfo(&free_memory, &total_memory);
            const std::size_t buffer_count = 3;
            if ((buffer_count * element_count * sizeof(Data64)) > static_cast<std::size_t>(free_memory * 0.70)){
                continue;
            }

            heongpu::DeviceVector<Data64> gpuntt_data(element_count);
            heongpu::DeviceVector<Data64> phantom_data(element_count);
            heongpu::DeviceVector<Data64> phantom_no_batch_data(element_count);

            cudaStream_t stream;
            cudaStreamCreate(&stream);
            cudaMemsetAsync(gpuntt_data.data(), 0, element_count * sizeof(Data64), stream);
            cudaMemsetAsync(phantom_data.data(), 0, element_count * sizeof(Data64), stream);
            cudaMemsetAsync(phantom_no_batch_data.data(), 0, element_count * sizeof(Data64), stream);
            cudaStreamSynchronize(stream);

            gpuntt::ntt_rns_configuration<Data64> cfg = {
                .n_power = n_power,
                .ntt_type = inverse ? gpuntt::INVERSE : gpuntt::FORWARD,
                .ntt_layout = gpuntt::PerPolynomial,
                .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
                .zero_padding = false,
                .mod_inverse = device_n_inverse.data(),
                .stream = stream};

            auto run_gpuntt = [&]() {
                if (inverse){
                    gpuntt::GPU_INTT_Inplace(
                        gpuntt_data.data(), device_inverse_roots.data(),
                        device_moduli.data(), cfg, batch_size, mod_count);
                }
                else{
                    gpuntt::GPU_NTT_Inplace(
                        gpuntt_data.data(), device_forward_roots.data(),
                        device_moduli.data(), cfg, batch_size, mod_count);
                }
            };
            auto run_phantom = [&]() {
                if (inverse){
                    heongpu::primitive::ntt_inverse_inplace(
                        phantom_data.data(), device_inverse_roots.data(),
                        device_moduli.data(), cfg, batch_size, mod_count,
                        phantom_tables);
                }
                else{
                    heongpu::primitive::ntt_forward_inplace(
                        phantom_data.data(), device_forward_roots.data(),
                        device_moduli.data(), cfg, batch_size, mod_count,
                        phantom_tables);
                }
            };

            auto run_phantom_no_batch = [&]() {
                for (int poly = 0; poly < poly_count; ++poly){
                    const std::size_t offset =
                        static_cast<std::size_t>(poly) * mod_count *
                        parameter.poly_modulus_degree;
                    if (inverse)
                    {
                        heongpu::primitive::ntt_inverse_inplace(
                            phantom_no_batch_data.data() + offset,
                            device_inverse_roots.data(), device_moduli.data(),
                            cfg, mod_count, mod_count, phantom_tables);
                    }
                    else
                    {
                        heongpu::primitive::ntt_forward_inplace(
                            phantom_no_batch_data.data() + offset,
                            device_forward_roots.data(), device_moduli.data(),
                            cfg, mod_count, mod_count, phantom_tables);
                    }
                }
            };

            float gpuntt_ms = 0.0F;
            float phantom_ms = 0.0F;
            float phantom_no_batch_ms = 0.0F;
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

            switch (results.size() % 3)
            {
                case 0:
                    gpuntt_ms = measure_gpuntt();
                    phantom_ms = measure_phantom();
                    phantom_no_batch_ms = measure_phantom_no_batch();
                    break;
                case 1:
                    phantom_ms = measure_phantom();
                    phantom_no_batch_ms = measure_phantom_no_batch();
                    gpuntt_ms = measure_gpuntt();
                    break;
                default:
                    phantom_no_batch_ms = measure_phantom_no_batch();
                    gpuntt_ms = measure_gpuntt();
                    phantom_ms = measure_phantom();
                    break;
            }
            cudaStreamDestroy(stream);

            results.push_back(
                {parameter.label, inverse ? "inverse_inplace" : "forward_inplace",
                 parameter.poly_modulus_degree, poly_count, mod_count,
                 batch_size, gpuntt_ms, phantom_ms, phantom_no_batch_ms});
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
                            const std::string& transform, int repeat_count)
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
    if (transform == "forward_inplace" || transform == "inverse_inplace"){
        std::cout << "Modes: phantom_batch - batch poly operation; phantom_no_batch - one poly at a time." << std::endl;
    }
 
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

    std::cout << std::left << std::setw(58) << "parameter / N / mod_count / mode";
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

        const int mod_count =
            row_results.empty() ? 0 : row_results.front().mod_count;
        const std::string row_prefix =
            current_label + " / N=" + std::to_string(current_n) +
            " / M=" + std::to_string(mod_count);

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

        if (transform == "forward_inplace" || transform == "inverse_inplace")
        {
            print_row("phantom_batch", false);
            print_row("phantom_no_batch", true);
        }
        else
        {
            print_row("phantom_no_batch", false);
        }
    }
}

void print_operation_summary_table(
    const std::vector<BenchmarkResult>& gpuntt_results,
    const std::vector<BenchmarkResult>& phantom_results)
{
    std::cout << "\n=================== HEOperator operation speedup summary ==================="<< std::endl;
    std::cout << "Measured HEOperator paths that call primitive::ntt_*:" << std::endl;
    std::cout << "  relin_II: HEArithmeticOperator::relinearize -> "
                 "HEOperator::relinearize_external_product_method2_inplace_ckks" << std::endl;
    std::cout << "  rotate_col_II: GPU rotate_col -> HEArithmeticOperator::rotate_rows -> "
                 "HEOperator::apply_galois_ckks_method_II" << std::endl;
    std::cout << "  rescale_leveled: HEArithmeticOperator::rescale -> HEOperator::rescale_inplace_ckks_leveled" << std::endl;
    std::cout << "Values are GPUNTT ms / PhantomNTT ms." << std::endl;
    std::cout << "Interpretation: >1.00 PhantomNTT faster, <1.00 GPUNTT faster." << std::endl;
    std::cout << std::left << std::setw(34) << "parameter / N"
              << std::right << std::setw(12) << "relin_II"
              << std::right << std::setw(16) << "rotate_col_II"
              << std::right << std::setw(16) << "rescale_leveled"
              << std::endl;
    std::cout << std::string(78, '-')
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
                  << std::defaultfloat << std::endl;
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
        run_ntt_kernel_comparison(parameters, repeat_count, false);
    print_ntt_kernel_table(forward_ntt_results, "forward_inplace", repeat_count);

    const auto inverse_ntt_results =
        run_ntt_kernel_comparison(parameters, repeat_count, true);
    print_ntt_kernel_table(inverse_ntt_results, "inverse_inplace", repeat_count);

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
            gpuntt = run_case(parameter, false, repeat_count);
            phantom = run_case(parameter, true, repeat_count);
        }
        else
        {
            phantom = run_case(parameter, true, repeat_count);
            gpuntt = run_case(parameter, false, repeat_count);
        }

        gpuntt_results.push_back(gpuntt);
        phantom_results.push_back(phantom);
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
