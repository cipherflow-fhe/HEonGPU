// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

#include "common.hpp"

#include <heongpu/kernel/switchkey.cuh>
#include <heongpu/ntt/ntt.cuh>
#include <heongpu/primitive/ntt.cuh>
#include <heongpu/primitive/switchkey.cuh>
#include <heongpu/switchkey/switchkey.cuh>
#include <heongpu/util/util.cuh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using ntt_plugin::build_key_modulus;
using ntt_plugin::make_context;
using ntt_plugin::make_bootstrapping_config;
using ntt_plugin::make_parameters;
using ntt_plugin::make_ntt_order;
using ntt_plugin::make_bsgs_input;
using ntt_plugin::ntt_surface_name;
using ntt_plugin::set_optimization_env;
using ntt_plugin::sort_parameters_by_n;
using ntt_plugin::NttSurface;
using ntt_plugin::ParameterSet;
using ntt_plugin::Scheme;
using heongpu::CudaException;

struct CheckResult
{
    std::string group;
    std::string name;
    std::string parameter;
    bool passed = false;
    double max_error = 0.0;
    std::size_t mismatches = 0;
};

enum class OperatorSurface
{
    Relin,
    Rotate,
    Bootstrap
};

const char* operator_surface_name(OperatorSurface surface)
{
    switch (surface)
    {
        case OperatorSurface::Relin:
            return "relinearize";
        case OperatorSurface::Rotate:
            return "rotate_rows";
        case OperatorSurface::Bootstrap:
            return "regular_bootstrapping_v2";
    }
    return "unknown";
}

std::vector<Data64> make_ntt_input(std::size_t n, int batch_size)
{
    std::vector<Data64> data(n * static_cast<std::size_t>(batch_size));
    for (int row = 0; row < batch_size; ++row)
    {
        for (std::size_t coeff = 0; coeff < n; ++coeff)
        {
            data[static_cast<std::size_t>(row) * n + coeff] =
                static_cast<Data64>(((coeff + 1) * 1315423911ULL +
                                     static_cast<std::size_t>(row) * 2654435761ULL) &
                                    0xFFFFFULL);
        }
    }
    return data;
}

std::size_t count_row_mismatches(const std::vector<Data64>& lhs,
                                 const std::vector<Data64>& rhs,
                                 const std::vector<Modulus64>& moduli,
                                 std::size_t n,
                                 int mod_count,
                                 int start_mod_idx = 0,
                                 const std::vector<int>* order = nullptr)
{
    std::size_t mismatches = 0;
    for (std::size_t idx = 0; idx < lhs.size(); ++idx)
    {
        const int row = static_cast<int>(idx / n);
        int mod_idx = start_mod_idx + (row % mod_count);
        if (order && !order->empty())
        {
            mod_idx = (*order)[static_cast<std::size_t>(row % mod_count)];
        }
        const Data64 modulus = moduli[static_cast<std::size_t>(mod_idx)].value;
        if ((lhs[idx] % modulus) != (rhs[idx] % modulus))
        {
            ++mismatches;
        }
    }
    return mismatches;
}

std::string route_parameter_label(const ParameterSet& parameter,
                                  bool use_phantom_ntt)
{
    return parameter.label + (use_phantom_ntt ? " plugin" : " gpuntt");
}

std::size_t count_modmajor_expanded_mismatches(
    const std::vector<Data64>& original,
    const std::vector<Data64>& modmajor,
    const std::vector<Modulus64>& moduli,
    std::size_t n,
    int decomp_count,
    int mod_count)
{
    std::size_t mismatches = 0;
    for (int group = 0; group < decomp_count; ++group)
    {
        for (int mod = 0; mod < mod_count; ++mod)
        {
            const Data64 modulus = moduli[static_cast<std::size_t>(mod)].value;
            const std::size_t original_row =
                (static_cast<std::size_t>(group) * mod_count + mod) * n;
            const std::size_t modmajor_row =
                (static_cast<std::size_t>(mod) * decomp_count + group) * n;
            for (std::size_t coeff = 0; coeff < n; ++coeff)
            {
                if ((original[original_row + coeff] % modulus) !=
                    (modmajor[modmajor_row + coeff] % modulus))
                {
                    ++mismatches;
                }
            }
        }
    }
    return mismatches;
}

CheckResult run_ntt_correctness(const ParameterSet& parameter,
                                NttSurface surface,
                                bool batch)
{
    const auto moduli = build_key_modulus(parameter);
    const int total_mod_count = static_cast<int>(moduli.size());
    const int mod_count =
        surface == NttSurface::PolyOrderedInverse ? 1 : total_mod_count;
    const int poly_count = 2;
    const int batch_size = mod_count * poly_count;
    const int start_mod_idx =
        surface == NttSurface::PolyOrderedInverse ? total_mod_count - 1 : 0;
    const int n_power =
        static_cast<int>(std::log2(parameter.poly_modulus_degree));
    const std::size_t n = parameter.poly_modulus_degree;
    const std::size_t element_count = n * static_cast<std::size_t>(batch_size);

    const auto roots_base =
        heongpu::generate_primitive_root_of_unity(n, moduli);
    const auto forward_roots =
        heongpu::generate_ntt_table(roots_base, moduli, n_power);
    const auto inverse_roots =
        heongpu::generate_intt_table(roots_base, moduli, n_power);
    const auto n_inverse = heongpu::generate_n_inverse(n, moduli);

    heongpu::DeviceVector<Modulus64> device_moduli(moduli);
    heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);
    heongpu::DeviceVector<Root64> device_inverse_roots(inverse_roots);
    heongpu::DeviceVector<Ninverse64> device_n_inverse(n_inverse);

    const std::vector<int> order =
        make_ntt_order(surface, mod_count, batch_size);
    std::vector<int> device_order_values =
        order.empty() ? std::vector<int>{0} : order;
    heongpu::DeviceVector<int> device_order(device_order_values);

    const auto tables =
        heongpu::ntt::make_phantom_ntt_tables_from_heongpu_roots(
            moduli, forward_roots, inverse_roots, n_inverse, n_power, 0,
            surface == NttSurface::ModulusOrderedForward);

    const auto input = make_ntt_input(n, batch_size);
    heongpu::DeviceVector<Data64> gpuntt_data(element_count);
    heongpu::DeviceVector<Data64> plugin_data(element_count);
    heongpu::DeviceVector<Data64> gpuntt_output(element_count);
    heongpu::DeviceVector<Data64> plugin_output(element_count);
    HEONGPU_CUDA_CHECK(cudaMemcpy(gpuntt_data.data(), input.data(),
                                  element_count * sizeof(Data64),
                                  cudaMemcpyHostToDevice));
    HEONGPU_CUDA_CHECK(cudaMemcpy(plugin_data.data(), input.data(),
                                  element_count * sizeof(Data64),
                                  cudaMemcpyHostToDevice));

    gpuntt::ntt_rns_configuration<Data64> cfg = {
        .n_power = n_power,
        .ntt_type = surface == NttSurface::ForwardInplace ||
                            surface == NttSurface::ModulusOrderedForward
                        ? gpuntt::FORWARD
                        : gpuntt::INVERSE,
        .ntt_layout = gpuntt::PerPolynomial,
        .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
        .zero_padding = false,
        .mod_inverse = device_n_inverse.data() + start_mod_idx,
        .stream = 0};

    if (surface == NttSurface::ForwardInplace)
    {
        gpuntt::GPU_NTT_Inplace(
            gpuntt_data.data(), device_forward_roots.data(),
            device_moduli.data(), cfg, batch_size, mod_count);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
        heongpu::ntt::phantom_ntt_inplace(
            plugin_data.data(), cfg, batch_size, mod_count, tables, batch);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }
    else if (surface == NttSurface::InverseInplace)
    {
        gpuntt::GPU_INTT_Inplace(
            gpuntt_data.data(), device_inverse_roots.data(),
            device_moduli.data(), cfg, batch_size, mod_count);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
        heongpu::ntt::phantom_intt_inplace(
            plugin_data.data(), cfg, batch_size, mod_count, tables, batch);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }
    else if (surface == NttSurface::InverseOutOfPlace)
    {
        gpuntt::GPU_INTT(
            gpuntt_data.data(), gpuntt_output.data(),
            device_inverse_roots.data(), device_moduli.data(), cfg,
            batch_size, mod_count);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
        heongpu::ntt::phantom_intt(
            plugin_data.data(), plugin_output.data(), cfg, batch_size,
            mod_count, tables, batch);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }
    else if (surface == NttSurface::ModulusOrderedForward)
    {
        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            gpuntt_data.data(), device_forward_roots.data(),
            device_moduli.data(), cfg, batch_size, mod_count,
            device_order.data());
        HEONGPU_CUDA_CHECK(cudaGetLastError());
        heongpu::ntt::phantom_ntt_modulus_ordered_inplace(
            plugin_data.data(), cfg, batch_size, mod_count, device_order.data(),
            tables, batch);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }
    else if (surface == NttSurface::ModulusOrderedInverse)
    {
        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            gpuntt_data.data(), device_inverse_roots.data(),
            device_moduli.data(), cfg, batch_size, mod_count,
            device_order.data());
        HEONGPU_CUDA_CHECK(cudaGetLastError());
        heongpu::ntt::phantom_intt_modulus_ordered_inplace(
            plugin_data.data(), cfg, batch_size, mod_count, device_order.data(),
            tables, batch);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }
    else
    {
        gpuntt::GPU_NTT_Poly_Ordered_Inplace(
            gpuntt_data.data(),
            device_inverse_roots.data() +
                (static_cast<std::size_t>(start_mod_idx) << n_power),
            device_moduli.data() + start_mod_idx, cfg, batch_size, mod_count,
            device_order.data());
        HEONGPU_CUDA_CHECK(cudaGetLastError());
        if (batch)
        {
            heongpu::ntt::phantom_intt_poly_ordered_inplace_batched(
                plugin_data.data(), cfg, batch_size, mod_count,
                device_order.data(), start_mod_idx, tables);
        }
        else
        {
            for (int row = 0; row < batch_size; ++row)
            {
                heongpu::ntt::phantom_intt_poly_ordered_inplace_batched(
                    plugin_data.data(), cfg, 1, 1, device_order.data() + row,
                    start_mod_idx, tables);
            }
        }
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }

    HEONGPU_CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Data64> gpuntt_host(element_count);
    std::vector<Data64> plugin_host(element_count);
    const Data64* gpuntt_source = surface == NttSurface::InverseOutOfPlace
                                      ? gpuntt_output.data()
                                      : gpuntt_data.data();
    const Data64* plugin_source = surface == NttSurface::InverseOutOfPlace
                                      ? plugin_output.data()
                                      : plugin_data.data();
    HEONGPU_CUDA_CHECK(cudaMemcpy(gpuntt_host.data(), gpuntt_source,
                                  element_count * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(plugin_host.data(), plugin_source,
                                  element_count * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaDeviceSynchronize());

    const std::size_t mismatches = count_row_mismatches(
        gpuntt_host, plugin_host, moduli, n, mod_count, start_mod_idx,
        &order);
    return {"isolated NTT",
            std::string(ntt_surface_name(surface)) +
                (batch ? " / batch" : " / no_batch"),
            parameter.label,
            mismatches == 0, static_cast<double>(mismatches), mismatches};
}

CheckResult run_mod_kswitch_correctness(const ParameterSet& parameter)
{
    const auto moduli = build_key_modulus(parameter);
    const int current_q_size =
        !parameter.q.empty() ? static_cast<int>(parameter.q.size())
                             : static_cast<int>(parameter.q_bit_sizes.size());
    const int current_qtilda_size = static_cast<int>(moduli.size());
    const int d = current_q_size;
    const int level = 0;
    const int iteration_count1 = d / 4;
    const int iteration_count2 = d % 4;
    const int n_power =
        static_cast<int>(std::log2(parameter.poly_modulus_degree));
    const std::size_t n = parameter.poly_modulus_degree;
    const std::size_t q_elements =
        n * static_cast<std::size_t>(current_q_size);
    const std::size_t expanded_elements =
        n * static_cast<std::size_t>(d) *
        static_cast<std::size_t>(current_qtilda_size);
    const std::size_t output_elements =
        n * static_cast<std::size_t>(2 * current_qtilda_size);
    const std::size_t relinkey_elements =
        n * static_cast<std::size_t>(2 * current_qtilda_size) *
        static_cast<std::size_t>(d);

    const auto roots_base =
        heongpu::generate_primitive_root_of_unity(n, moduli);
    const auto forward_roots =
        heongpu::generate_ntt_table(roots_base, moduli, n_power);
    const auto inverse_roots =
        heongpu::generate_intt_table(roots_base, moduli, n_power);
    const auto n_inverse = heongpu::generate_n_inverse(n, moduli);

    heongpu::DeviceVector<Modulus64> device_moduli(moduli);
    heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);

    const auto tables =
        heongpu::ntt::make_phantom_ntt_tables_from_heongpu_roots(
            moduli, forward_roots, inverse_roots, n_inverse, n_power, 0, true);

    std::vector<Data64> ciphertext_coeff(q_elements);
    for (int mod = 0; mod < current_q_size; ++mod)
    {
        const Data64 modulus = moduli[static_cast<std::size_t>(mod)].value;
        for (std::size_t coeff = 0; coeff < n; ++coeff)
        {
            ciphertext_coeff[static_cast<std::size_t>(mod) * n + coeff] =
                static_cast<Data64>(
                    (17ULL + coeff * 13ULL +
                     static_cast<std::size_t>(mod) * 257ULL) %
                    (modulus >> 2));
        }
    }

    std::vector<Data64> base_change_matrix(
        static_cast<std::size_t>(d) * current_qtilda_size, 1);
    std::vector<Data64> mi_inv(static_cast<std::size_t>(current_q_size), 1);
    std::vector<Data64> prod(
        static_cast<std::size_t>(d) * current_qtilda_size, 0);
    std::vector<int> group_sizes(static_cast<std::size_t>(d), 1);
    std::vector<int> group_locations(static_cast<std::size_t>(d));
    std::vector<int> mod_index(static_cast<std::size_t>(current_qtilda_size));
    std::vector<int> order(static_cast<std::size_t>(current_qtilda_size));
    for (int i = 0; i < d; ++i)
    {
        group_locations[static_cast<std::size_t>(i)] = i;
    }
    for (int i = 0; i < current_qtilda_size; ++i)
    {
        mod_index[static_cast<std::size_t>(i)] = i;
        order[static_cast<std::size_t>(i)] = i;
    }

    std::vector<Data64> relinkey(relinkey_elements);
    for (int group = 0; group < d; ++group)
    {
        for (int component = 0; component < 2; ++component)
        {
            for (int mod = 0; mod < current_qtilda_size; ++mod)
            {
                const Data64 modulus =
                    moduli[static_cast<std::size_t>(mod)].value;
                const std::size_t row =
                    (((static_cast<std::size_t>(group) * 2 + component) *
                          current_qtilda_size +
                      mod)
                     * n);
                for (std::size_t coeff = 0; coeff < n; ++coeff)
                {
                    relinkey[row + coeff] =
                        static_cast<Data64>(
                            (31ULL + coeff * 7ULL +
                             static_cast<std::size_t>(mod) * 41ULL +
                             static_cast<std::size_t>(group) * 97ULL +
                             static_cast<std::size_t>(component) * 503ULL) %
                            (modulus >> 2));
                }
            }
        }
    }

    heongpu::DeviceVector<Data64> device_ciphertext_coeff(ciphertext_coeff);
    heongpu::DeviceVector<Data64> device_ciphertext_ntt(ciphertext_coeff);
    heongpu::DeviceVector<Data64> device_original_expanded(expanded_elements);
    heongpu::DeviceVector<Data64> device_modmajor_expanded(expanded_elements);
    heongpu::DeviceVector<Data64> device_original_output(output_elements);
    heongpu::DeviceVector<Data64> device_modmajor_output(output_elements);
    heongpu::DeviceVector<Data64> device_base_change_matrix(base_change_matrix);
    heongpu::DeviceVector<Data64> device_mi_inv(mi_inv);
    heongpu::DeviceVector<Data64> device_prod(prod);
    heongpu::DeviceVector<Data64> device_relinkey(relinkey);
    heongpu::DeviceVector<int> device_group_sizes(group_sizes);
    heongpu::DeviceVector<int> device_group_locations(group_locations);
    heongpu::DeviceVector<int> device_mod_index(mod_index);
    heongpu::DeviceVector<int> device_order(order);
    HEONGPU_CUDA_CHECK(cudaMemset(device_original_expanded.data(), 0,
                                  expanded_elements * sizeof(Data64)));
    HEONGPU_CUDA_CHECK(cudaMemset(device_modmajor_expanded.data(), 0,
                                  expanded_elements * sizeof(Data64)));
    HEONGPU_CUDA_CHECK(cudaMemset(device_original_output.data(), 0,
                                  output_elements * sizeof(Data64)));
    HEONGPU_CUDA_CHECK(cudaMemset(device_modmajor_output.data(), 0,
                                  output_elements * sizeof(Data64)));

    gpuntt::ntt_rns_configuration<Data64> q_cfg_ntt = {
        .n_power = n_power,
        .ntt_type = gpuntt::FORWARD,
        .ntt_layout = gpuntt::PerPolynomial,
        .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
        .zero_padding = false,
        .stream = 0};
    gpuntt::GPU_NTT_Inplace(
        device_ciphertext_ntt.data(), device_forward_roots.data(),
        device_moduli.data(), q_cfg_ntt, current_q_size, current_q_size);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    gpuntt::ntt_rns_configuration<Data64> qp_cfg_ntt = {
        .n_power = n_power,
        .ntt_type = gpuntt::FORWARD,
        .ntt_layout = gpuntt::PerPolynomial,
        .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
        .zero_padding = false,
        .stream = 0};

    heongpu::base_conversion_DtoQtilde_relin_leveled_kernel<<<
        dim3((1 << n_power) >> 8, d, 1), 256>>>(
        device_ciphertext_coeff.data(), device_original_expanded.data(),
        device_moduli.data(), device_base_change_matrix.data(),
        device_mi_inv.data(), device_prod.data(), device_group_sizes.data(),
        device_group_locations.data(), n_power, d, current_qtilda_size,
        current_q_size, level, device_mod_index.data());
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
        device_original_expanded.data(), device_forward_roots.data(),
        device_moduli.data(), qp_cfg_ntt, d * current_qtilda_size,
        current_qtilda_size, device_order.data());
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
        dim3((1 << n_power) >> 8, current_qtilda_size, 1), 256>>>(
        device_original_expanded.data(), device_relinkey.data(),
        device_original_output.data(), device_moduli.data(),
        current_qtilda_size, current_qtilda_size, current_qtilda_size,
        iteration_count1, iteration_count2, level, n_power);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    heongpu::switchkey::base_conversion_DtoQtilde_relin_leveled_modmajor(
        device_ciphertext_coeff.data(), device_ciphertext_ntt.data(),
        device_modmajor_expanded.data(), device_moduli.data(),
        device_base_change_matrix.data(), device_mi_inv.data(),
        device_prod.data(), device_group_sizes.data(),
        device_group_locations.data(), n_power, d, current_qtilda_size,
        current_q_size, level, true, 0);
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::ntt::phantom_ntt_modmajor_inplace(
        device_modmajor_expanded.data(), qp_cfg_ntt, current_qtilda_size,
        current_q_size, current_q_size, d, device_group_sizes.data(),
        device_group_locations.data(), tables, true);
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::switchkey::keyswitch_multiply_accumulate_leveled_method_II_modmajor(
        device_modmajor_expanded.data(), device_relinkey.data(),
        device_modmajor_output.data(), device_moduli.data(),
        current_qtilda_size, current_qtilda_size, current_qtilda_size,
        iteration_count1, iteration_count2, level, n_power, 0);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    HEONGPU_CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Data64> original_expanded_host(expanded_elements);
    std::vector<Data64> modmajor_expanded_host(expanded_elements);
    std::vector<Data64> original_host(output_elements);
    std::vector<Data64> modmajor_host(output_elements);
    HEONGPU_CUDA_CHECK(cudaMemcpy(original_expanded_host.data(),
                                  device_original_expanded.data(),
                                  expanded_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(modmajor_expanded_host.data(),
                                  device_modmajor_expanded.data(),
                                  expanded_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(original_host.data(),
                                  device_original_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(modmajor_host.data(),
                                  device_modmajor_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaDeviceSynchronize());

    const std::size_t expanded_mismatches =
        count_modmajor_expanded_mismatches(
            original_expanded_host, modmajor_expanded_host, moduli, n, d,
            current_qtilda_size);
    const std::size_t output_mismatches = count_row_mismatches(
        original_host, modmajor_host, moduli, n, current_qtilda_size);
    const std::size_t mismatches = expanded_mismatches + output_mismatches;
    return {"KeySwitch_P1", "baseconv_ntt_keymul", parameter.label,
            mismatches == 0, static_cast<double>(mismatches), mismatches};
}

CheckResult run_bsgs_fusion_correctness(const ParameterSet& parameter)
{
    const auto moduli = build_key_modulus(parameter);
    const int limb_count = static_cast<int>(moduli.size());
    const int n_power =
        static_cast<int>(std::log2(parameter.poly_modulus_degree));
    const std::size_t n = parameter.poly_modulus_degree;
    const int component_count = 2;
    const std::size_t limb_elements =
        n * static_cast<std::size_t>(limb_count);
    const std::size_t ct_elements =
        limb_elements * static_cast<std::size_t>(component_count);
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
            moduli, forward_roots, inverse_roots, n_inverse, n_power, 0, true);

    const auto input =
        make_bsgs_input(n, component_count * limb_count, limb_count, moduli,
                        17);
    const auto addend =
        make_bsgs_input(n, limb_count, limb_count, moduli, 31);
    const auto accum =
        make_bsgs_input(n, component_count * limb_count, limb_count, moduli,
                        43);

    heongpu::DeviceVector<Modulus64> device_moduli(moduli);
    heongpu::DeviceVector<Data64> baby_base_input(input);
    heongpu::DeviceVector<Data64> baby_fused_input(input);
    heongpu::DeviceVector<Data64> giant_base_input(input);
    heongpu::DeviceVector<Data64> giant_fused_input(input);
    heongpu::DeviceVector<Data64> device_addend(addend);
    heongpu::DeviceVector<Data64> device_accum(accum);
    heongpu::DeviceVector<Data64> baby_base_output(ct_elements);
    heongpu::DeviceVector<Data64> baby_fused_output(ct_elements);
    heongpu::DeviceVector<Data64> giant_base_output(ct_elements);
    heongpu::DeviceVector<Data64> giant_fused_output(ct_elements);
    heongpu::DeviceVector<Data64> scratch(ct_elements);

    setenv("HEONGPU_USE_BSGS_FUSION", "0", 1);
    heongpu::primitive::bs_add_permute_fused(
        baby_base_input.data(), device_addend.data(), baby_base_output.data(),
        device_moduli.data(), galois_elt, n_power, limb_count, tables, 0);
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::primitive::gs_add_permute_acc_fused(
        giant_base_input.data(), device_addend.data(), device_accum.data(),
        scratch.data(), giant_base_output.data(), device_moduli.data(),
        galois_elt, n_power, limb_count, tables, 0);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    setenv("HEONGPU_USE_BSGS_FUSION", "1", 1);
    heongpu::primitive::bs_add_permute_fused(
        baby_fused_input.data(), device_addend.data(),
        baby_fused_output.data(), device_moduli.data(), galois_elt, n_power,
        limb_count, tables, 0);
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::primitive::gs_add_permute_acc_fused(
        giant_fused_input.data(), device_addend.data(), device_accum.data(),
        scratch.data(), giant_fused_output.data(), device_moduli.data(),
        galois_elt, n_power, limb_count, tables, 0);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    HEONGPU_CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Data64> baby_base_result(ct_elements);
    std::vector<Data64> baby_fused_result(ct_elements);
    std::vector<Data64> giant_base_result(ct_elements);
    std::vector<Data64> giant_fused_result(ct_elements);
    HEONGPU_CUDA_CHECK(cudaMemcpy(baby_base_result.data(),
                                  baby_base_output.data(),
                                  ct_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(baby_fused_result.data(),
                                  baby_fused_output.data(),
                                  ct_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(giant_base_result.data(),
                                  giant_base_output.data(),
                                  ct_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(giant_fused_result.data(),
                                  giant_fused_output.data(),
                                  ct_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));

    const std::size_t baby_mismatches = count_row_mismatches(
        baby_base_result, baby_fused_result, moduli, n, limb_count);
    const std::size_t giant_mismatches = count_row_mismatches(
        giant_base_result, giant_fused_result, moduli, n, limb_count);
    const std::size_t mismatches = baby_mismatches + giant_mismatches;
    return {"BSGS fusion", "baby_step + giant_step", parameter.label,
            mismatches == 0, static_cast<double>(mismatches), mismatches};
}

CheckResult run_keyswitch_part2_correctness(const ParameterSet& parameter)
{
    const int q_size = !parameter.q.empty()
                           ? static_cast<int>(parameter.q.size())
                           : static_cast<int>(parameter.q_bit_sizes.size());
    const int p_size = !parameter.p.empty()
                           ? static_cast<int>(parameter.p.size())
                           : static_cast<int>(parameter.p_bit_sizes.size());
    if (q_size == 0 || p_size == 0)
    {
        return {"KeySwitch_Part2", "part2", parameter.label,
                false, 1.0, 1};
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
            moduli, forward_roots, inverse_roots, n_inverse, n_power, 0, true);

    const auto input =
        make_bsgs_input(n, 2 * q_prime_size, q_prime_size, moduli, 71);
    const auto addend = make_bsgs_input(n, q_size, q_size, moduli, 97);

    heongpu::DeviceVector<Modulus64> device_moduli(moduli);
    heongpu::DeviceVector<Root64> device_forward_roots(forward_roots);
    heongpu::DeviceVector<Data64> device_half(half);
    heongpu::DeviceVector<Data64> device_half_mod(half_mod);
    heongpu::DeviceVector<Data64> device_last_q_modinv(last_q_modinv);
    heongpu::DeviceVector<Data64> device_input(input);
    heongpu::DeviceVector<Data64> device_addend(addend);

    heongpu::DeviceVector<Data64> baseline_gpuntt_div(output_elements);
    heongpu::DeviceVector<Data64> baseline_gpuntt_addend(q_elements);
    heongpu::DeviceVector<Data64> baseline_output(output_elements);
    heongpu::DeviceVector<Data64> baseline_phantom_div(output_elements);
    heongpu::DeviceVector<Data64> baseline_phantom_addend(q_elements);
    heongpu::DeviceVector<Data64> baseline_phantom_output(output_elements);
    heongpu::DeviceVector<Data64> original_wrapper_addend(q_elements);
    heongpu::DeviceVector<Data64> original_wrapper_output(output_elements);
    heongpu::DeviceVector<Data64> original_wrapper_scratch(q_elements);
    heongpu::DeviceVector<Data64> keyswitch_part2_addend(q_elements);
    heongpu::DeviceVector<Data64> keyswitch_part2_output(output_elements);
    heongpu::DeviceVector<Data64> keyswitch_part2_scratch(q_elements);
    heongpu::DeviceVector<Data64> keyswitch_part2_phantom_output(output_elements);
    heongpu::DeviceVector<Data64> keyswitch_part2_phantom_scratch(q_elements);

    gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
        .n_power = n_power,
        .ntt_type = gpuntt::FORWARD,
        .ntt_layout = gpuntt::PerPolynomial,
        .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
        .zero_padding = false,
        .stream = 0};

    HEONGPU_CUDA_CHECK(cudaMemcpy(baseline_gpuntt_addend.data(),
                                  device_addend.data(),
                                  q_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToDevice));
    heongpu::divide_round_lastq_extended_leveled_kernel<<<
        dim3(n >> 8, q_size, 2), 256>>>(
        device_input.data(), baseline_gpuntt_div.data(), device_moduli.data(),
        device_half.data(), device_half_mod.data(), device_last_q_modinv.data(),
        n_power, q_prime_size, q_size, q_prime_size, q_size, p_size);
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::primitive::NTT_inplace(
        baseline_gpuntt_div.data(), device_forward_roots.data(),
        device_moduli.data(), cfg_ntt, 2 * q_size, q_size, nullptr);
    heongpu::primitive::NTT_inplace(
        baseline_gpuntt_addend.data(), device_forward_roots.data(),
        device_moduli.data(), cfg_ntt, q_size, q_size, nullptr);
    heongpu::addition_switchkey<<<dim3(n >> 8, q_size, 2), 256>>>(
        baseline_gpuntt_div.data(), baseline_gpuntt_addend.data(),
        baseline_output.data(), device_moduli.data(), n_power);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    HEONGPU_CUDA_CHECK(cudaMemcpy(baseline_phantom_addend.data(),
                                  device_addend.data(),
                                  q_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToDevice));
    heongpu::divide_round_lastq_extended_leveled_kernel<<<
        dim3(n >> 8, q_size, 2), 256>>>(
        device_input.data(), baseline_phantom_div.data(), device_moduli.data(),
        device_half.data(), device_half_mod.data(), device_last_q_modinv.data(),
        n_power, q_prime_size, q_size, q_prime_size, q_size, p_size);
    HEONGPU_CUDA_CHECK(cudaGetLastError());
    heongpu::primitive::NTT_inplace(
        baseline_phantom_div.data(), device_forward_roots.data(),
        device_moduli.data(), cfg_ntt, 2 * q_size, q_size, tables);
    heongpu::primitive::NTT_inplace(
        baseline_phantom_addend.data(), device_forward_roots.data(),
        device_moduli.data(), cfg_ntt, q_size, q_size, tables);
    heongpu::addition_switchkey<<<dim3(n >> 8, q_size, 2), 256>>>(
        baseline_phantom_div.data(), baseline_phantom_addend.data(),
        baseline_phantom_output.data(), device_moduli.data(), n_power);
    HEONGPU_CUDA_CHECK(cudaGetLastError());

    HEONGPU_CUDA_CHECK(cudaMemcpy(original_wrapper_addend.data(),
                                  device_addend.data(),
                                  q_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToDevice));
    setenv("HEONGPU_USE_KSWITCH_P2", "0", 1);
    heongpu::primitive::keyswitch_part2_fused_moddown_ntt(
        device_input.data(), original_wrapper_addend.data(),
        original_wrapper_scratch.data(), original_wrapper_output.data(),
        device_forward_roots.data(), device_moduli.data(), cfg_ntt,
        device_half.data(), device_half_mod.data(),
        device_last_q_modinv.data(), n_power, q_prime_size, q_size,
        q_prime_size, q_size, p_size, nullptr, 0);

    HEONGPU_CUDA_CHECK(cudaMemcpy(keyswitch_part2_addend.data(),
                                  device_addend.data(),
                                  q_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToDevice));
    setenv("HEONGPU_USE_KSWITCH_P2", "1", 1);
    heongpu::primitive::keyswitch_part2_fused_moddown_ntt(
        device_input.data(), keyswitch_part2_addend.data(),
        keyswitch_part2_scratch.data(), keyswitch_part2_output.data(),
        device_forward_roots.data(),
        device_moduli.data(), cfg_ntt, device_half.data(),
        device_half_mod.data(), device_last_q_modinv.data(), n_power,
        q_prime_size, q_size, q_prime_size, q_size, p_size, nullptr, 0);

    setenv("HEONGPU_USE_KSWITCH_P2", "1", 1);
    heongpu::primitive::keyswitch_part2_fused_moddown_ntt(
        device_input.data(), device_addend.data(),
        keyswitch_part2_phantom_scratch.data(),
        keyswitch_part2_phantom_output.data(), device_forward_roots.data(),
        device_moduli.data(), cfg_ntt, device_half.data(),
        device_half_mod.data(), device_last_q_modinv.data(), n_power,
        q_prime_size, q_size, q_prime_size, q_size, p_size, tables, 0);

    HEONGPU_CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<Data64> baseline_host(output_elements);
    std::vector<Data64> baseline_phantom_host(output_elements);
    std::vector<Data64> original_wrapper_host(output_elements);
    std::vector<Data64> keyswitch_part2_host(output_elements);
    std::vector<Data64> keyswitch_part2_phantom_host(output_elements);
    HEONGPU_CUDA_CHECK(cudaMemcpy(baseline_host.data(), baseline_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(baseline_phantom_host.data(),
                                  baseline_phantom_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(original_wrapper_host.data(),
                                  original_wrapper_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(keyswitch_part2_host.data(),
                                  keyswitch_part2_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));
    HEONGPU_CUDA_CHECK(cudaMemcpy(keyswitch_part2_phantom_host.data(),
                                  keyswitch_part2_phantom_output.data(),
                                  output_elements * sizeof(Data64),
                                  cudaMemcpyDeviceToHost));

    const std::size_t phantom_mismatches = count_row_mismatches(
        baseline_host, baseline_phantom_host, moduli, n, q_size);
    const std::size_t original_wrapper_mismatches = count_row_mismatches(
        baseline_host, original_wrapper_host, moduli, n, q_size);
    const std::size_t keyswitch_part2_mismatches = count_row_mismatches(
        baseline_host, keyswitch_part2_host, moduli, n, q_size);
    const std::size_t keyswitch_part2_phantom_mismatches = count_row_mismatches(
        baseline_host, keyswitch_part2_phantom_host, moduli, n, q_size);
    const std::size_t mismatches =
        phantom_mismatches + original_wrapper_mismatches +
        keyswitch_part2_mismatches + keyswitch_part2_phantom_mismatches;

    return {"KeySwitch_Part2", "part2", parameter.label,
            mismatches == 0, static_cast<double>(mismatches), mismatches};
}

std::vector<Complex64> make_message(int slot_count)
{
    std::vector<Complex64> message(static_cast<std::size_t>(slot_count));
    for (int i = 0; i < slot_count; ++i)
    {
        const double real = 0.01 * static_cast<double>((i % 13) + 1);
        const double imag = 0.002 * static_cast<double>((i % 7) - 3);
        message[static_cast<std::size_t>(i)] = Complex64(real, imag);
    }
    return message;
}

double max_complex_error(const std::vector<Complex64>& lhs,
                         const std::vector<Complex64>& rhs)
{
    const std::size_t count = std::min(lhs.size(), rhs.size());
    double max_error = lhs.size() == rhs.size()
                           ? 0.0
                           : std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < count; ++i)
    {
        const double real_error = lhs[i].real() - rhs[i].real();
        const double imag_error = lhs[i].imag() - rhs[i].imag();
        max_error = std::max(
            max_error, std::sqrt(real_error * real_error +
                                 imag_error * imag_error));
    }
    return max_error;
}

std::vector<Complex64> square_message(const std::vector<Complex64>& message)
{
    std::vector<Complex64> expected(message.size());
    for (std::size_t i = 0; i < message.size(); ++i)
    {
        expected[i] = message[i] * message[i];
    }
    return expected;
}

std::vector<Complex64> rotate_message(const std::vector<Complex64>& message,
                                      int shift)
{
    std::vector<Complex64> expected(message.size());
    const int slots = static_cast<int>(message.size());
    for (int i = 0; i < slots; ++i)
    {
        const int source = (i + shift + slots) % slots;
        expected[static_cast<std::size_t>(i)] =
            message[static_cast<std::size_t>(source)];
    }
    return expected;
}

std::vector<Complex64> decode_cipher(
    heongpu::HEContext<Scheme> context,
    heongpu::HEEncoder<Scheme>& encoder,
    heongpu::HEDecryptor<Scheme>& decryptor,
    heongpu::Ciphertext<Scheme>& ciphertext)
{
    heongpu::Plaintext<Scheme> plaintext(context);
    decryptor.decrypt(plaintext, ciphertext);
    std::vector<Complex64> decoded;
    encoder.decode(decoded, plaintext);
    cudaDeviceSynchronize();
    return decoded;
}

CheckResult run_operator_correctness(const ParameterSet& parameter,
                                     OperatorSurface surface,
                                     bool use_phantom_ntt,
                                     bool use_bsgs_fusion)
{
    set_optimization_env(use_phantom_ntt, use_phantom_ntt, use_bsgs_fusion);
    const bool bootstrap = surface == OperatorSurface::Bootstrap;
    heongpu::HEContext<Scheme> context =
        make_context(parameter, use_phantom_ntt, bootstrap);

    const int secret_weight = 192;
    const int ephemeral_secret_weight = 32;
    heongpu::HEKeyGenerator<Scheme> keygen(context);
    heongpu::Secretkey<Scheme> secret_key =
        bootstrap ? heongpu::Secretkey<Scheme>(context, secret_weight)
                  : heongpu::Secretkey<Scheme>(context);
    if (bootstrap)
    {
        keygen.generate_secret_key_v2(secret_key);
    }
    else
    {
        keygen.generate_secret_key(secret_key);
    }

    heongpu::Publickey<Scheme> public_key(context);
    keygen.generate_public_key(public_key, secret_key);

    heongpu::Relinkey<Scheme> relin_key(context);
    keygen.generate_relin_key(relin_key, secret_key);

    std::vector<int> key_index = {1};
    heongpu::Switchkey<Scheme> swk_dense_to_sparse(context);
    heongpu::Switchkey<Scheme> swk_sparse_to_dense(context);

    heongpu::HEEncoder<Scheme> encoder(context);
    heongpu::HEEncryptor<Scheme> encryptor(context, public_key);
    heongpu::HEDecryptor<Scheme> decryptor(context, secret_key);
    heongpu::HEArithmeticOperator<Scheme> operators(context, encoder);

    if (bootstrap)
    {
        heongpu::Secretkey<Scheme> sparse_secret_key(context,
                                                     ephemeral_secret_weight);
        keygen.generate_secret_key_v2(sparse_secret_key);
        keygen.generate_switch_key(swk_dense_to_sparse, sparse_secret_key,
                                   secret_key);
        keygen.generate_switch_key(swk_sparse_to_dense, secret_key,
                                   sparse_secret_key);

        const auto boot_config =
            make_bootstrapping_config(context->get_key_modulus()[0].value);
        operators.generate_bootstrapping_params_v2(parameter.scale,
                                                   boot_config);
        key_index = operators.bootstrapping_key_indexs();
    }

    heongpu::Galoiskey<Scheme> galois_key(context, key_index);
    keygen.generate_galois_key(galois_key, secret_key);

    const int slot_count = static_cast<int>(parameter.poly_modulus_degree / 2);
    std::vector<Complex64> message =
        bootstrap ? std::vector<Complex64>(
                        static_cast<std::size_t>(slot_count),
                        Complex64(0.2, 0.4))
                  : make_message(slot_count);

    heongpu::Plaintext<Scheme> plaintext(context);
    encoder.encode(plaintext, message, parameter.scale);

    heongpu::Ciphertext<Scheme> input(context);
    encryptor.encrypt(input, plaintext);

    std::vector<Complex64> expected;
    heongpu::Ciphertext<Scheme> output(context);
    if (surface == OperatorSurface::Relin)
    {
        operators.multiply(input, input, output);
        operators.relinearize_inplace(output, relin_key);
        expected = square_message(message);
    }
    else if (surface == OperatorSurface::Rotate)
    {
        operators.rotate_rows(input, output, galois_key, 1);
        const auto rotate_left = rotate_message(message, 1);
        const auto rotate_right = rotate_message(message, -1);
        const auto decoded = decode_cipher(context, encoder, decryptor, output);
        const double left_error = max_complex_error(decoded, rotate_left);
        const double right_error = max_complex_error(decoded, rotate_right);
        const double rotation_error = std::min(left_error, right_error);
        constexpr double tolerance = 0.1;
        return {"CKKS operator", operator_surface_name(surface),
                route_parameter_label(parameter, use_phantom_ntt),
                rotation_error <= tolerance, rotation_error, 0};
    }
    else
    {
        for (int level = 1; level < static_cast<int>(parameter.q.size());
             ++level)
        {
            operators.mod_drop_inplace(input);
        }
        output = operators.regular_bootstrapping_v2(
            input, galois_key, relin_key, &swk_dense_to_sparse,
            &swk_sparse_to_dense);
        expected = message;
    }

    const auto decoded = decode_cipher(context, encoder, decryptor, output);
    const double max_error = max_complex_error(decoded, expected);
    const double tolerance = bootstrap ? 0.75 : 0.1;
    return {"CKKS operator", operator_surface_name(surface),
            route_parameter_label(parameter, use_phantom_ntt),
            max_error <= tolerance, max_error, 0};
}

void print_results(const std::vector<CheckResult>& results)
{
    std::cout << std::left << std::setw(22) << "group"
              << std::setw(38) << "test"
              << std::setw(32) << "parameter"
              << std::right << std::setw(14) << "max/error"
              << std::setw(14) << "mismatch"
              << std::setw(10) << "status" << std::endl;
    std::cout << std::string(130, '-') << std::endl;

    for (const auto& result : results)
    {
        std::cout << std::left << std::setw(22) << result.group
                  << std::setw(38) << result.name
                  << std::setw(32) << result.parameter
                  << std::right << std::setw(14) << std::scientific
                  << std::setprecision(3) << result.max_error
                  << std::setw(14) << result.mismatches
                  << std::setw(10) << (result.passed ? "PASS" : "FAIL")
                  << std::defaultfloat << std::endl;
    }
}

int main()
{
    std::vector<ParameterSet> parameters = make_parameters();
    sort_parameters_by_n(parameters);

    std::vector<CheckResult> results;
    const NttSurface ntt_surfaces[] = {
        NttSurface::ForwardInplace,
        NttSurface::InverseInplace,
        NttSurface::InverseOutOfPlace,
        NttSurface::ModulusOrderedForward,
        NttSurface::ModulusOrderedInverse,
        NttSurface::PolyOrderedInverse};

    for (const auto& parameter : parameters)
    {
        for (NttSurface surface : ntt_surfaces)
        {
            results.push_back(run_ntt_correctness(parameter, surface, true));
            results.push_back(run_ntt_correctness(parameter, surface, false));
        }
        results.push_back(run_mod_kswitch_correctness(parameter));
        results.push_back(run_bsgs_fusion_correctness(parameter));
        results.push_back(run_keyswitch_part2_correctness(parameter));
    }

    const OperatorSurface operator_surfaces[] = {
        OperatorSurface::Relin,
        OperatorSurface::Rotate};

    for (const auto& parameter : parameters)
    {
        if (parameter.poly_modulus_degree > 65536)
        {
            continue;
        }
        for (OperatorSurface surface : operator_surfaces)
        {
            results.push_back(
                run_operator_correctness(parameter, surface, false, false));
            results.push_back(
                run_operator_correctness(parameter, surface, true, true));
        }
    }

    const auto boot_parameter =
        std::find_if(parameters.begin(), parameters.end(),
                     [](const ParameterSet& parameter) {
                         return parameter.label == "N16QP1546H192H32";
                     });
    if (boot_parameter != parameters.end())
    {
        results.push_back(run_operator_correctness(
            *boot_parameter, OperatorSurface::Bootstrap, false, false));
        results.push_back(run_operator_correctness(
            *boot_parameter, OperatorSurface::Bootstrap, true, true));
    }

    print_results(results);

    const bool all_passed =
        std::all_of(results.begin(), results.end(),
                    [](const CheckResult& result) { return result.passed; });
    std::cout << "\nNTT plugin correctness: "
              << (all_passed ? "PASS" : "FAIL") << std::endl;
    return all_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
