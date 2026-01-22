/**
 * @company CipherFlow
 */

#include "heongpu.hpp"
#include "ckks/precision.cuh"
#include <gtest/gtest.h>
#include <random>
#include <iomanip>
#include <algorithm>

// Minimum precision requirement (in bits)
constexpr double MIN_PRECISION = 15.0;

// Helper function to apply bit-reverse permutation
template<typename T>
void bit_reverse_inplace(std::vector<T>& vec, size_t slots)
{
    for (size_t i = 0; i < slots; i++)
    {
        size_t j = 0;
        size_t temp = i;
        size_t log_slots = 0;
        size_t n = slots;
        while (n > 1) {
            log_slots++;
            n >>= 1;
        }

        for (size_t k = 0; k < log_slots; k++)
        {
            j = (j << 1) | (temp & 1);
            temp >>= 1;
        }

        if (i < j)
        {
            std::swap(vec[i], vec[j]);
        }
    }
}


TEST(HEonGPU, CKKS_CoeffsToSlots_FullPacking)
{
    cudaSetDevice(0);

    {
        size_t poly_modulus_degree = 8192;
        heongpu::HEContext<heongpu::Scheme::CKKS> context(
            heongpu::keyswitching_type::KEYSWITCHING_METHOD_II,
            heongpu::sec_level_type::none);
        context.set_poly_modulus_degree(poly_modulus_degree);

        context.set_coeff_modulus_bit_sizes(
            {60, 40, 56, 56, 56, 56},  // Q moduli
            {60, 60});          // P moduli
        context.generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        const int slot_count = poly_modulus_degree / 2;
        double scale = pow(2.0, 40);

        // Generate bootstrapping parameters for CoeffsToSlots
        heongpu::BootstrappingConfig boot_config(
            heongpu::EncodingMatrixConfig(heongpu::Linear_Transform_Type::SlotsToCoeffs, 0),
            heongpu::EvalModConfig(0),  // Not using eval_mod in this test
            heongpu::EncodingMatrixConfig(heongpu::Linear_Transform_Type::CoeffsToSlots, 5)
        );
        operators.generate_bootstrapping_params(scale, boot_config);

        std::vector<int> key_index = operators.bootstrapping_key_indexs();
        heongpu::Galoiskey<heongpu::Scheme::CKKS> galois_key(context, key_index);
        keygen.generate_galois_key(galois_key, secret_key);

        // Generate random complex values
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        std::vector<Complex64> values(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            values[i] = Complex64(dis(gen), dis(gen));
        }

        // Split into real and imaginary parts
        std::vector<Complex64> values_real(slot_count);
        std::vector<Complex64> values_imag(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            values_real[i] = Complex64(values[i].real(), 0.0);
            values_imag[i] = Complex64(values[i].imag(), 0.0);
        }

        // Apply bit-reverse on the original complex vector
        bit_reverse_inplace(values, slot_count);

        // Map to a float vector (coefficient representation)
        // Real parts in first half, imaginary parts in second half
        std::vector<double> values_float(poly_modulus_degree, 0.0);
        for (int i = 0; i < slot_count; i++)
        {
            values_float[i] = values[i].real();
            values_float[i + slot_count] = values[i].imag();
        }

        // Encode coefficient-wise and encrypt
        heongpu::Plaintext<heongpu::Scheme::CKKS> plaintext(context);
        encoder.encode_coeff(plaintext, values_float, scale);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> ciphertext(context);
        encryptor.encrypt(ciphertext, plaintext);

        // Apply homomorphic CoeffsToSlots (DFT)
        heongpu::ExecutionOptions options;
        std::vector<heongpu::Ciphertext<heongpu::Scheme::CKKS>> result_cts =
            operators.coeff_to_slot_cf(ciphertext, galois_key, options);

        std::cout << "\n=== CoeffsToSlots returned " << result_cts.size() << " ciphertexts ===" << std::endl;

        // Decrypt and decode the first ciphertext (real part)
        heongpu::Plaintext<heongpu::Scheme::CKKS> result_pt_real(context);
        decryptor.decrypt(result_pt_real, result_cts[0]);

        std::vector<Complex64> output_real;
        encoder.decode(output_real, result_pt_real);

        // Extract first slot_count values from real part
        std::vector<Complex64> output_real_slots(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            output_real_slots[i] = output_real[i];
        }

        // Decrypt and decode the second ciphertext (imaginary part)
        heongpu::Plaintext<heongpu::Scheme::CKKS> result_pt_imag(context);
        decryptor.decrypt(result_pt_imag, result_cts[1]);

        std::vector<Complex64> output_imag;
        encoder.decode(output_imag, result_pt_imag);

        // Extract first slot_count values from imaginary part
        std::vector<Complex64> output_imag_slots(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            output_imag_slots[i] = output_imag[i];
        }

        // Combine real and imaginary parts
        std::vector<Complex64> output_combined(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            output_combined[i] = Complex64(
                output_real_slots[i].real(),
                output_imag_slots[i].real()
            );
        }

        // Print first few values for debugging
        std::cout << "\n=== CoeffsToSlots Test: First 10 values ===" << std::endl;
        for (int i = 0; i < std::min(10, slot_count); i++)
        {
            std::cout << "output[" << i << "] = " << output_combined[i].real()
                      << " + " << output_combined[i].imag() << "i, expected = "
                      << values_real[i].real() << " + " << values_imag[i].real() << "i" << std::endl;
        }

        // Compute precision statistics for real part
        std::vector<Complex64> expected_real_only(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            expected_real_only[i] = Complex64(values_real[i].real(), 0.0);
        }

        heongpu::PrecisionStats prec_stats_real =
            heongpu::get_precision_stats(expected_real_only, output_real_slots);

        std::cout << "\n=== CoeffsToSlots Real Part Precision ===" << std::endl;
        std::cout << prec_stats_real.to_string() << std::endl;

        // Verify real part precision meets minimum requirements
        EXPECT_GE(prec_stats_real.mean_precision.real, MIN_PRECISION)
            << "Mean real precision " << prec_stats_real.mean_precision.real
            << " is below minimum " << MIN_PRECISION << " bits";

        // Compute precision statistics for imaginary part
        std::vector<Complex64> expected_imag_only(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            expected_imag_only[i] = Complex64(values_imag[i].real(), 0.0);
        }

        heongpu::PrecisionStats prec_stats_imag =
            heongpu::get_precision_stats(expected_imag_only, output_imag_slots);

        std::cout << "\n=== CoeffsToSlots Imaginary Part Precision ===" << std::endl;
        std::cout << prec_stats_imag.to_string() << std::endl;

        // Verify imaginary part precision meets minimum requirements
        EXPECT_GE(prec_stats_imag.mean_precision.real, MIN_PRECISION)
            << "Mean imaginary precision " << prec_stats_imag.mean_precision.real
            << " is below minimum " << MIN_PRECISION << " bits";

        // Compute precision statistics for combined output
        std::vector<Complex64> expected_combined(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            expected_combined[i] = Complex64(values_real[i].real(), values_imag[i].real());
        }

        heongpu::PrecisionStats prec_stats_combined =
            heongpu::get_precision_stats(expected_combined, output_combined);

        std::cout << "\n=== CoeffsToSlots Combined Precision ===" << std::endl;
        std::cout << prec_stats_combined.to_string() << std::endl;

        // Verify combined precision meets minimum requirements
        EXPECT_GE(prec_stats_combined.mean_precision.real, MIN_PRECISION)
            << "Mean combined real precision " << prec_stats_combined.mean_precision.real
            << " is below minimum " << MIN_PRECISION << " bits";

        EXPECT_GE(prec_stats_combined.mean_precision.imag, MIN_PRECISION)
            << "Mean combined imaginary precision " << prec_stats_combined.mean_precision.imag
            << " is below minimum " << MIN_PRECISION << " bits";
    }

    cudaDeviceSynchronize();
}

TEST(HEonGPU, CKKS_SlotsToCoeffs_FullPacking)
{
    cudaSetDevice(0);

    {
        size_t poly_modulus_degree = 8192;
        heongpu::HEContext<heongpu::Scheme::CKKS> context(
            heongpu::keyswitching_type::KEYSWITCHING_METHOD_II,
            heongpu::sec_level_type::none);
        context.set_poly_modulus_degree(poly_modulus_degree);

        context.set_coeff_modulus_bit_sizes(
            {60, 40, 39, 39, 39},  // Q moduli
            {60, 60});          // P moduli
        context.generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        const int slot_count = poly_modulus_degree / 2;
        double scale = pow(2.0, 40);

        // Generate bootstrapping parameters for SlotsToCoeffs
        heongpu::BootstrappingConfig boot_config(
            heongpu::EncodingMatrixConfig(heongpu::Linear_Transform_Type::SlotsToCoeffs, 4),
            heongpu::EvalModConfig(0),  // Not using eval_mod in this test
            heongpu::EncodingMatrixConfig(heongpu::Linear_Transform_Type::CoeffsToSlots, 0)
        );
        operators.generate_bootstrapping_params(scale, boot_config);

        std::vector<int> key_index = operators.bootstrapping_key_indexs();
        heongpu::Galoiskey<heongpu::Scheme::CKKS> galois_key(context, key_index);
        keygen.generate_galois_key(galois_key, secret_key);

        // Generate test vectors (real values only, as encoding always results in real coefficients)
        std::vector<Complex64> values_real(slot_count);
        std::vector<Complex64> values_imag(slot_count);

        for (int i = 0; i < slot_count; i++)
        {
            values_real[i] = Complex64(static_cast<double>(i + 1) / slot_count, 0.0);
            values_imag[i] = Complex64(static_cast<double>(i + 1) / slot_count, 0.0);
        }

        // Encode and encrypt the test vectors
        heongpu::Plaintext<heongpu::Scheme::CKKS> plaintext(context);
        encoder.encode(plaintext, values_real, scale);
        heongpu::Ciphertext<heongpu::Scheme::CKKS> ct_real(context);
        encryptor.encrypt(ct_real, plaintext);

        // For full packing, we need a separate ciphertext for imaginary part
        encoder.encode(plaintext, values_imag, scale);
        heongpu::Ciphertext<heongpu::Scheme::CKKS> ct_imag(context);
        encryptor.encrypt(ct_imag, plaintext);

        // Apply homomorphic SlotsToCoeffs (inverse DFT)
        heongpu::ExecutionOptions options;
        heongpu::Ciphertext<heongpu::Scheme::CKKS> result =
            operators.slot_to_coeff_cf(ct_real, ct_imag, galois_key, options);
       
        heongpu::Plaintext<heongpu::Scheme::CKKS> result_pt(context);
        decryptor.decrypt(result_pt, result);

        std::vector<double> coeffs_float;
        encoder.decode_coeff(coeffs_float, result_pt);

        // Extract the coefficients and construct the complex vector
        std::vector<Complex64> values_test(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            values_test[i] = Complex64(
                coeffs_float[i],
                coeffs_float[i + slot_count]
            );
        }

        // // Combine real and imaginary parts for reference
        std::vector<Complex64> values_expected(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            values_expected[i] = Complex64(
                values_real[i].real(),
                values_imag[i].real()
            );
        }

        // Result is bit-reversed, so apply bit-reverse on reference
        bit_reverse_inplace(values_expected, slot_count);

        // Print first few values for debugging
        std::cout << "\n=== SlotsToCoeffs Test: First 10 values ===" << std::endl;
        for (int i = 0; i < std::min(10, slot_count); i++)
        {
            std::cout << "output[" << i << "] = " << values_test[i].real()
                      << " + " << values_test[i].imag() << "i, expected = "
                      << values_expected[i].real() << " + " << values_expected[i].imag() << "i" << std::endl;
        }

        // Compute precision statistics
        heongpu::PrecisionStats prec_stats =
            heongpu::get_precision_stats(values_expected, values_test);

        std::cout << "\n=== SlotsToCoeffs Precision Statistics ===" << std::endl;
        std::cout << prec_stats.to_string() << std::endl;

        // Verify precision meets minimum requirements
        EXPECT_GE(prec_stats.mean_precision.real, MIN_PRECISION)
            << "Mean real precision " << prec_stats.mean_precision.real
            << " is below minimum " << MIN_PRECISION << " bits";

        EXPECT_GE(prec_stats.mean_precision.imag, MIN_PRECISION)
            << "Mean imaginary precision " << prec_stats.mean_precision.imag
            << " is below minimum " << MIN_PRECISION << " bits";
    }

    cudaDeviceSynchronize();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
