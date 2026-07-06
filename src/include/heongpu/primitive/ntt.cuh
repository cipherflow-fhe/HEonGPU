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

    void NTT_inplace(
        Data64* data,
        Root64* gpuntt_roots,
        Modulus64* moduli,
        gpuntt::ntt_rns_configuration<Data64> cfg,
        int batch_size,
        int mod_count,
        const std::shared_ptr<PhantomNttTables>& phantom_tables);

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

} // namespace primitive
} // namespace heongpu

#endif
