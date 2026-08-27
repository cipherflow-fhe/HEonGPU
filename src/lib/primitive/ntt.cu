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
            ntt::phantom_ntt_inplace(
                data, cfg, batch_size, mod_count, tables, true);
            return;
        }

        ntt::phantom_ntt_inplace(
            data, cfg, batch_size, mod_count, tables, false);
    }

    void NTT_modulus_ordered_inplace(
        Data64* data,
        Root64* roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        int* order,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
                data, roots, moduli, cfg, batch_size, mod_count, order);
            return;
        }

        if (cfg.n_power < 16)
        {
            ntt::phantom_ntt_modulus_ordered_inplace(
                data, cfg, batch_size, mod_count, order, tables, true);
            return;
        }

        ntt::phantom_ntt_modulus_ordered_inplace(
            data, cfg, batch_size, mod_count, order, tables, false);
    }

    void NTT_modmajor_inplace(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int current_qp,
        int current_q,
        int first_q_count,
        int decomp_count,
        int* group_sizes,
        int* group_locations,
        const std::shared_ptr<PhantomNttTables>& tables,
        bool skip_excluded)
    {
        ntt::phantom_ntt_modmajor_inplace(
            data, cfg, current_qp, current_q, first_q_count, decomp_count,
            group_sizes, group_locations, tables, skip_excluded);
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
            ntt::phantom_intt_inplace(
                data, cfg, batch_size, mod_count, tables, true);
            return;
        }

        ntt::phantom_intt_inplace(
            data, cfg, batch_size, mod_count, tables, false);
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
            ntt::phantom_intt(
                input, output, cfg, batch_size, mod_count, tables, true);
            return;
        }

        ntt::phantom_intt(
            input, output, cfg, batch_size, mod_count, tables, false);
    }

    void INTT_modulus_ordered_inplace(
        Data64* data,
        Root64* roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        int* order,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
                data, roots, moduli, cfg, batch_size, mod_count, order);
            return;
        }

        if (cfg.n_power < 16)
        {
            ntt::phantom_intt_modulus_ordered_inplace(
                data, cfg, batch_size, mod_count, order, tables, true);
            return;
        }

        ntt::phantom_intt_modulus_ordered_inplace(
            data, cfg, batch_size, mod_count, order, tables, false);
    }

    void INTT_poly_ordered_inplace(
        Data64* data,
        Root64* roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        int* order,
        int start_mod_idx,
        const std::shared_ptr<PhantomNttTables>& tables)
    {
        if (!tables)
        {
            gpuntt::GPU_NTT_Poly_Ordered_Inplace(
                data, roots, moduli, cfg, batch_size, mod_count, order);
            return;
        }

        ntt::phantom_intt_poly_ordered_inplace_batched(
            data, cfg, batch_size, mod_count, order, start_mod_idx, tables);
    }
} // namespace primitive
} // namespace heongpu
