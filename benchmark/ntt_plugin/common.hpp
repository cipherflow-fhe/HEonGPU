// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <heongpu/heongpu.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace ntt_plugin
{
    constexpr auto Scheme = heongpu::Scheme::CKKS;

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

    enum class NttSurface
    {
        ForwardInplace,
        InverseInplace,
        InverseOutOfPlace,
        ModulusOrderedForward,
        ModulusOrderedInverse,
        PolyOrderedInverse
    };

    inline const char* ntt_surface_name(NttSurface surface)
    {
        switch (surface)
        {
            case NttSurface::ForwardInplace:
                return "NTT_inplace";
            case NttSurface::InverseInplace:
                return "INTT_inplace";
            case NttSurface::InverseOutOfPlace:
                return "INTT";
            case NttSurface::ModulusOrderedForward:
                return "NTT_modulus_ordered_inplace";
            case NttSurface::ModulusOrderedInverse:
                return "INTT_modulus_ordered_inplace";
            case NttSurface::PolyOrderedInverse:
                return "INTT_poly_ordered_inplace";
        }
        return "unknown";
    }

    inline const char* ntt_surface_title(NttSurface surface)
    {
        switch (surface)
        {
            case NttSurface::ForwardInplace:
                return "NTT inplace";
            case NttSurface::InverseInplace:
                return "INTT inplace";
            case NttSurface::InverseOutOfPlace:
                return "INTT";
            case NttSurface::ModulusOrderedForward:
                return "Modulus_ordered NTT";
            case NttSurface::ModulusOrderedInverse:
                return "Modulus_ordered INTT";
            case NttSurface::PolyOrderedInverse:
                return "Poly_ordered INTT";
        }
        return "unknown NTT";
    }

    inline std::vector<int> make_ntt_order(NttSurface surface, int mod_count,
                                           int batch_size)
    {
        std::vector<int> order;
        if (surface == NttSurface::ModulusOrderedForward ||
            surface == NttSurface::ModulusOrderedInverse)
        {
            order.resize(static_cast<std::size_t>(mod_count));
            for (int i = 0; i < mod_count; ++i)
            {
                order[static_cast<std::size_t>(i)] = i;
            }
            std::rotate(order.begin(), order.begin() + 1, order.end());
        }
        else if (surface == NttSurface::PolyOrderedInverse)
        {
            order.resize(static_cast<std::size_t>(batch_size));
            for (int i = 0; i < batch_size; ++i)
            {
                order[static_cast<std::size_t>(i)] = batch_size - 1 - i;
            }
        }
        return order;
    }

    inline std::vector<ParameterSet> make_parameters()
    {
        return {
            {.label = "PN12QP109",
             .poly_modulus_degree = 4096,
             .q = {0x200000e001ULL, 0x100006001ULL},
             .p = {0x3ffffea001ULL},
             .scale = std::ldexp(1.0, 32)},
            {.label = "PN13QP218",
             .poly_modulus_degree = 8192,
             .q = {0x7fffffffeac001ULL, 0x3fffffffeb8001ULL,
                   0x3fffffffef8001ULL},
             .p = {0x7ffffffffb4001ULL},
             .scale = std::ldexp(1.0, 50)},
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
             .q = {0x4000000120001ULL, 0x10000140001ULL,
                   0xffffe80001ULL, 0x10000290001ULL, 0xffffc40001ULL,
                   0x100003e0001ULL, 0x10000470001ULL, 0x100004b0001ULL,
                   0xffffb20001ULL, 0x10000500001ULL, 0x10000650001ULL,
                   0xffff940001ULL, 0xffff8a0001ULL, 0xffff820001ULL,
                   0xffff780001ULL, 0x10000890001ULL, 0xffff750001ULL,
                   0x10000960001ULL},
             .p = {0x40000001b0001ULL, 0x3ffffffdf0001ULL,
                   0x4000000270001ULL},
             .scale = std::ldexp(1.0, 40)},
            {.label = "PN16QP1761",
             .poly_modulus_degree = 65536,
             .q = {0x80000000080001ULL, 0x2000000a0001ULL,
                   0x2000000e0001ULL, 0x1fffffc20001ULL,
                   0x200000440001ULL, 0x200000500001ULL,
                   0x200000620001ULL, 0x1fffff980001ULL,
                   0x2000006a0001ULL, 0x1fffff7e0001ULL,
                   0x200000860001ULL, 0x200000a60001ULL,
                   0x200000aa0001ULL, 0x200000b20001ULL,
                   0x200000c80001ULL, 0x1fffff360001ULL,
                   0x200000e20001ULL, 0x1fffff060001ULL,
                   0x200000fe0001ULL, 0x1ffffede0001ULL,
                   0x1ffffeca0001ULL, 0x1ffffeb40001ULL,
                   0x200001520001ULL, 0x1ffffe760001ULL,
                   0x2000019a0001ULL, 0x1ffffe640001ULL,
                   0x200001a00001ULL, 0x1ffffe520001ULL,
                   0x200001e80001ULL, 0x1ffffe0c0001ULL,
                   0x1ffffdee0001ULL, 0x200002480001ULL,
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

    inline void sort_parameters_by_n(std::vector<ParameterSet>& parameters)
    {
        std::sort(parameters.begin(), parameters.end(),
                  [](const ParameterSet& lhs, const ParameterSet& rhs) {
                      if (lhs.poly_modulus_degree !=
                          rhs.poly_modulus_degree)
                      {
                          return lhs.poly_modulus_degree <
                                 rhs.poly_modulus_degree;
                      }
                      return lhs.label < rhs.label;
                  });
    }

    inline std::vector<Modulus64> build_key_modulus(
        const ParameterSet& parameter)
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
        {
            moduli.emplace_back(value);
        }
        for (Data64 value : parameter.p)
        {
            moduli.emplace_back(value);
        }
        return moduli;
    }

    inline void set_coeff_modulus(heongpu::HEContext<Scheme>& context,
                                  const ParameterSet& parameter)
    {
        if (!parameter.q_bit_sizes.empty() || !parameter.p_bit_sizes.empty())
        {
            context->set_coeff_modulus_bit_sizes(parameter.q_bit_sizes,
                                                 parameter.p_bit_sizes);
        }
        else
        {
            context->set_coeff_modulus_values(parameter.q, parameter.p);
        }
    }

    inline heongpu::HEContext<Scheme> make_context(
        const ParameterSet& parameter, bool use_phantom_ntt,
        bool set_full_slot_count = false)
    {
        setenv("HEONGPU_USE_PHANTOM_NTT", use_phantom_ntt ? "1" : "0", 1);

        heongpu::HEContext<Scheme> context =
            heongpu::GenHEContext<Scheme>(heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(parameter.poly_modulus_degree);
        if (set_full_slot_count)
        {
            context->set_slot_count(
                static_cast<int>(parameter.poly_modulus_degree / 2));
        }
        set_coeff_modulus(context, parameter);
        context->generate();
        return context;
    }

    inline heongpu::BootstrappingConfigV2 make_bootstrapping_config(
        Data64 first_modulus)
    {
        heongpu::EvalModConfig eval_mod_config(
            first_modulus, 20, 256.0, 16, 30, 3, 0, std::ldexp(1.0, 60));

        return heongpu::BootstrappingConfigV2(
            heongpu::EncodingMatrixConfig(
                heongpu::LinearTransformType::SLOTS_TO_COEFFS, 12, 2.0, 3),
            eval_mod_config,
            heongpu::EncodingMatrixConfig(
                heongpu::LinearTransformType::COEFFS_TO_SLOTS, 24, 2.0, 4));
    }
} // namespace ntt_plugin
