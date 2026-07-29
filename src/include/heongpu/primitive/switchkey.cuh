#ifndef HEONGPU_PRIMITIVE_SWITCHKEY_CUH
#define HEONGPU_PRIMITIVE_SWITCHKEY_CUH

#include <cuda_runtime.h>
#include <memory>

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

    void base_conversion_DtoQtilde_relin_leveled_ntt(
        Data64* ciphertext_coeff,
        Data64* ciphertext_ntt,
        Data64* output,
        Root64* gpuntt_roots,
        Modulus64* modulus,
        gpuntt::ntt_rns_configuration<Data64> cfg_ntt,
        Data64* base_change_matrix_D_to_Qtilda,
        Data64* Mi_inv_D_to_Qtilda,
        Data64* prod_D_to_Qtilda,
        int* I_j,
        int* I_location,
        int n_power,
        int d,
        int current_Qtilda_size,
        int current_Q_size,
        int first_Q_size,
        int level,
        int* mod_index,
        int* order,
        bool copy_excluded,
        const std::shared_ptr<PhantomNttTables>& phantom_tables,
        cudaStream_t stream);

    void keyswitch_multiply_accumulate_leveled_method_II(
        Data64* input,
        const Data64* relinkey,
        Data64* output,
        Modulus64* modulus,
        int first_rns_mod_count,
        int current_decomp_mod_count,
        int current_rns_mod_count,
        int iteration_count1,
        int iteration_count2,
        int level,
        int n_power,
        const std::shared_ptr<PhantomNttTables>& phantom_tables,
        cudaStream_t stream);

} // namespace primitive
} // namespace heongpu

#endif
