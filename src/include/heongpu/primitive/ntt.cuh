#ifndef HEONGPU_PRIMITIVE_NTT_CUH
#define HEONGPU_PRIMITIVE_NTT_CUH

#include <cuda_runtime.h>
#include <memory>
#include <vector>

#include <gpuntt/ntt_merge/ntt.cuh>
#include <heongpu/util/util.cuh>

namespace heongpu
{
namespace ntt
{
    struct PhantomNttTables;
} // namespace ntt

namespace primitive
{
    using PhantomNttTables = heongpu::ntt::PhantomNttTables;

    std::shared_ptr<PhantomNttTables> make_phantom_ntt_tables_from_heongpu_roots(
        const std::vector<Modulus64>& moduli,
        const std::vector<Root64>& forward_roots,
        const std::vector<Root64>& inverse_roots,
        const std::vector<Ninverse64>& n_inverse,
        int n_power,
        cudaStream_t stream = 0,
        bool forward_only = false);

    /////////////////////////////////// Forward NTT ////////////////////////////////
    void NTT_inplace(
        Data64* data,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void NTT_modulus_ordered_inplace(
        Data64* data,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        int* order,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void NTT_modmajor_inplace(
        Data64* data,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int current_qp,
        int current_q,
        int first_q_count,
        int decomp_count,
        int* group_sizes,
        int* group_locations,
        const std::shared_ptr<PhantomNttTables>& phantom_tables,
        bool skip_excluded);

    /////////////////////////////////// Inverse NTT ////////////////////////////////
    void INTT_inplace(
        Data64* data,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void INTT(
        Data64* input,
        Data64* output,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void INTT_modulus_ordered_inplace(
        Data64* data,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        int* order,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

    void INTT_poly_ordered_inplace(
        Data64* data,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        int* order,
        int start_mod_idx,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);
} // namespace primitive
} // namespace heongpu

#endif
