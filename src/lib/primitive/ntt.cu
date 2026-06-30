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
namespace
{
    bool use_phantom_ntt(
        const std::shared_ptr<PhantomNttTables>& tables,
        const gpuntt::ntt_rns_configuration<Data64>& cfg,
        int batch_size, int mod_count)
    {
        return ntt::can_use_phantom_ntt(tables, cfg, batch_size, mod_count);
    }

} // namespace

    std::shared_ptr<PhantomNttTables> make_phantom_ntt_tables_from_heongpu_roots(
        const std::vector<Modulus64>& moduli,
        const std::vector<Root64>& forward_roots,
        const std::vector<Root64>& inverse_roots,
        const std::vector<Ninverse64>& n_inverse,
        int n_power,
        cudaStream_t stream)
    {
        return ntt::make_phantom_ntt_tables_from_heongpu_roots(
            moduli, forward_roots, inverse_roots, n_inverse, n_power, stream);
    }

    void ntt_forward_inplace(Data64* data, Root64* roots, Modulus64* moduli,
                             gpuntt::ntt_rns_configuration<Data64> cfg,
                             int batch_size, int mod_count,
                             const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (use_phantom_ntt(tables, cfg, batch_size, mod_count))
        {
            ntt::forward_phantom_ntt_inplace(
                data, cfg, batch_size, mod_count, tables);
            return;
        }

        gpuntt::GPU_NTT_Inplace(
            data, roots, moduli, cfg, batch_size, mod_count);
    }

    void ntt_inverse_inplace(Data64* data, Root64* roots, Modulus64* moduli,
                             gpuntt::ntt_rns_configuration<Data64> cfg,
                             int batch_size, int mod_count,
                             const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (use_phantom_ntt(tables, cfg, batch_size, mod_count))
        {
            ntt::inverse_phantom_ntt_inplace(
                data, cfg, batch_size, mod_count, tables);
            return;
        }

        gpuntt::GPU_INTT_Inplace(
            data, roots, moduli, cfg, batch_size, mod_count);
    }

} // namespace primitive
} // namespace heongpu
