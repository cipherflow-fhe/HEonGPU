// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#ifndef HEONGPU_NTT_CORE_H
#define HEONGPU_NTT_CORE_H

#include <cuda_runtime.h>
#include <memory>
#include <vector>

#include <gpuntt/ntt_merge/ntt.cuh>
#include <heongpu/util/util.cuh>

namespace heongpu
{
namespace ntt
{
    // Optional Activation: PhantomFHE NTT tables built from HEonGPU roots.
    struct PhantomNttTables;

    std::shared_ptr<PhantomNttTables> make_phantom_ntt_tables_from_heongpu_roots(
        const std::vector<Modulus64>& moduli,
        const std::vector<Root64>& forward_roots,
        const std::vector<Root64>& inverse_roots,
        const std::vector<Ninverse64>& n_inverse,
        int n_power,
        cudaStream_t stream = 0,
        bool forward_only = false);

    void phantom_ntt_inplace(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void phantom_ntt_inplace_batched(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void phantom_intt_inplace(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void phantom_intt_inplace_batched(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void phantom_intt(
        const Data64* input,
        Data64* output,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void phantom_intt_batched(
        const Data64* input,
        Data64* output,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

} // namespace ntt
} // namespace heongpu

#endif // HEONGPU_NTT_CORE_H
