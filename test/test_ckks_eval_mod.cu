/**
 * @company CipherFlow
 */

#include "heongpu.hpp"
#include "ckks/precision.cuh"
#include <gtest/gtest.h>
#include <random>
#include <iomanip>

// Minimum precision requirement (in bits)
constexpr double MIN_PRECISION = 15.0;

TEST(HEonGPU, CKKS_EvalMod)
{
    cudaSetDevice(0);

    {
        size_t poly_modulus_degree = 8192;
        heongpu::HEContext<heongpu::Scheme::CKKS> context(
            heongpu::keyswitching_type::KEYSWITCHING_METHOD_II,
            heongpu::sec_level_type::none);
        context.set_poly_modulus_degree(poly_modulus_degree);

        // Modulus configuration matching the example
        // context.set_coeff_modulus_bit_sizes(
        // {55,
        //     45,
        //     60, 60, 60, 60,
        //     60, 60, 60, 60,
        //     53
        //     },
        // {60, 60, 60});
        context.set_coeff_modulus_bit_sizes(
        {60, 40, 40, 40, 40, 40, 40, 40, 40, 40, //lv9
            39, 39, 39,  // stc
            60, 60, 60, 60, 60, 60, 60, 60, /// eval mod
            60}, // cts lv25
        {60, 60, 60, 60, 60});
        context.generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::Relinkey<heongpu::Scheme::CKKS> relin_key(context);
        keygen.generate_relin_key(relin_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context,
                                                              public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context,
                                                              secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context,
                                                                       encoder);

        // Generate bootstrapping parameters
        double scale = pow(2.0, 40);

        heongpu::EvalModConfig eval_mod_config(context.get_key_modulus()[0].value, 20, 256.0, 16, 30, 3, 0, pow(2.0, 60));

        heongpu::BootstrappingConfig boot_config(
            heongpu::EncodingMatrixConfig(heongpu::Linear_Transform_Type::CoeffsToSlots, 0),
            eval_mod_config,
            heongpu::EncodingMatrixConfig(heongpu::Linear_Transform_Type::SlotsToCoeffs, 0)
        );
        operators.generate_bootstrapping_params(scale, boot_config);

        const int slot_count = poly_modulus_degree / 2;

        // Create random values in range that will be reduced modulo Q
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        std::vector<Complex64> message(slot_count);
        std::vector<Complex64> expected(slot_count);

        for (int i = 0; i < slot_count; i++)
        {
            double val = dis(gen);
            message[i] = Complex64(val, 0.0);

            expected[i] = message[i] - Complex64(eval_mod_config.message_ratio_* eval_mod_config.q_diff_*std::round(message[i].real() / (eval_mod_config.message_ratio_/eval_mod_config.q_diff_)), 0.0);
            // expected[i] = message[i];
        }

        // Encode and encrypt
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        encoder.encode(P1, message, scale);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context);
        encryptor.encrypt(C1, P1);

        heongpu::ExecutionOptions options;

        double q0 = static_cast<double>(eval_mod_config.Q_);
    
        double scale1 = std::exp2(std::round(std::log2(q0 / eval_mod_config.message_ratio_)));
        double scale_up_factor1 = std::round(scale1/C1.scale());
        operators.scale_up(C1, scale_up_factor1, C1, options);

        double scale2 = std::round((eval_mod_config.scaling_factor_/ eval_mod_config.message_ratio_)/ C1.scale());
        operators.scale_up(C1, scale2, C1, options);
    
        Complex64 const1(
            1.0/(double(eval_mod_config.K_)*eval_mod_config.q_diff_), 0.0);
        
        std::cout << "const1 " << const1 << std::endl;
        // Perform eval_mod
        operators.mult_const(C1, const1, C1, options);
        operators.rescale_inplace(C1, options);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> result =
            operators.eval_mod(C1, relin_key, options);
        
        // Decrypt and decode
        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        decryptor.decrypt(P2, result);

        std::vector<Complex64> output;
        encoder.decode(output, P2);

        // Print first few output values for debugging
        std::cout << "\n=== First 10 output values ===" << std::endl;
        for (int i = 0; i < std::min(10, (int)output.size()); i++) {
            std::cout << "output[" << i << "] = " << output[i].real()
                    << " + " << output[i].imag() << "i, expected = "
                    << expected[i].real() << " + " << expected[i].imag() << "i" << std::endl;
        }
        
        // Compute precision statistics using get_precision_stats
        heongpu::PrecisionStats prec_stats =
            heongpu::get_precision_stats(expected, output);

        // Print statistics
        std::cout << "\n=== EvalMod Precision Statistics ===" << std::endl;
        std::cout << prec_stats.to_string() << std::endl;

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
