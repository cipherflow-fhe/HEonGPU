// Copyright 2026 CipherFlow
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer : Chan Jia Lin

#include <heongpu/primitive/ntt.cuh>
#include <heongpu/ntt/ntt.cuh>

namespace heongpu
{
namespace primitive
{
    std::shared_ptr<PhantomNttTables> make_phantom_ntt_tables_from_heongpu_roots(
        const std::vector<Modulus64>& moduli,
        const std::vector<Root64>& forward_roots,
        const std::vector<Root64>& inverse_roots,
        const std::vector<Ninverse64>& n_inverse,
        int n_power,
        cudaStream_t stream,
        bool forward_only)
    {
        return ntt::make_phantom_ntt_tables_from_heongpu_roots(
            moduli, forward_roots, inverse_roots, n_inverse, n_power, stream,
            forward_only);
    }

    void NTT_inplace(
        Data64* data,
        Root64* roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            gpuntt::GPU_NTT_Inplace(
                data, roots, moduli, cfg, batch_size, mod_count);
            return;
        }

        if (cfg.n_power <= 16)
        {
            ntt::phantom_ntt_inplace_batched(
                data, cfg, batch_size, mod_count, tables);
            return;
        }

        ntt::phantom_ntt_inplace(data, cfg, batch_size, mod_count, tables);
    }

    void INTT_inplace(
        Data64* data,
        Root64* roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            gpuntt::GPU_INTT_Inplace(
                data, roots, moduli, cfg, batch_size, mod_count);
            return;
        }

        if (cfg.n_power <= 16)
        {
            ntt::phantom_intt_inplace_batched(
                data, cfg, batch_size, mod_count, tables);
            return;
        }

        ntt::phantom_intt_inplace(data, cfg, batch_size, mod_count, tables);
    }

    void INTT(
        Data64* input,
        Data64* output,
        Root64* roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            gpuntt::GPU_INTT(
                input, output, roots, moduli, cfg, batch_size, mod_count);
            return;
        }

        if (cfg.n_power <= 16)
        {
            ntt::phantom_intt_batched(
                input, output, cfg, batch_size, mod_count, tables);
            return;
        }

        ntt::phantom_intt(input, output, cfg, batch_size, mod_count, tables);
    }

} // namespace primitive
} // namespace heongpu
