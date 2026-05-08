// Copyright 2024-2026 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include <heongpu/heongpu.hpp>
#include <gtest/gtest.h>

template <typename T>
bool fix_point_equal(T input1, T input2, T epsilon = static_cast<T>(1e-4))
{
    return std::fabs(input1 - input2) < epsilon;
}

template <typename T>
bool fix_point_array_check(const std::vector<T>& array1,
                           const std::vector<T>& array2,
                           T epsilon = static_cast<T>(1e-4))
{
    if (array1.size() != array2.size())
    {
        return false;
    }

    for (size_t i = 0; i < array1.size(); ++i)
    {
        if (!fix_point_equal(array1[i], array2[i], epsilon))
        {
            return false;
        }
    }

    return true;
}

TEST(HEonGPU, CKKS_Encoding_Decoding)
{
    {
        size_t poly_modulus_degree = 4096;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_coeff_modulus_bit_sizes({40, 30, 30}, {40});
        context->generate();

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        const int row_size = poly_modulus_degree / 2;
        std::vector<double> message(row_size, 0);
        for (int i = 0; i < row_size; i++)
        {
            message[i] = dis(gen);
        }

        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        double scale = pow(2.0, 30);
        encoder.encode(P1, message, scale);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P1);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(message, gpu_result), true);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        double number = static_cast<double>(dis(gen));
        encoder.encode(P2, number, scale);

        std::vector<double> gpu_result2;
        encoder.decode(gpu_result2, P2);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_equal(number, gpu_result2[0]), true);
    }

    cudaDeviceSynchronize();

    {
        size_t poly_modulus_degree = 8192;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        const int row_size = poly_modulus_degree / 2;
        std::vector<double> message(row_size, 0);
        for (int i = 0; i < row_size; i++)
        {
            message[i] = dis(gen);
        }

        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        double scale = pow(2.0, 30);
        encoder.encode(P1, message, scale);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P1);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(message, gpu_result), true);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        double number = static_cast<double>(dis(gen));
        encoder.encode(P2, number, scale);

        std::vector<double> gpu_result2;
        encoder.decode(gpu_result2, P2);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_equal(number, gpu_result2[0]), true);
    }

    cudaDeviceSynchronize();

    {
        size_t poly_modulus_degree = 16384;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_coeff_modulus_bit_sizes(
            {45, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35}, {45});
        context->generate();

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        const int row_size = poly_modulus_degree / 2;
        std::vector<double> message(row_size, 0);
        for (int i = 0; i < row_size; i++)
        {
            message[i] = dis(gen);
        }

        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        double scale = pow(2.0, 30);
        encoder.encode(P1, message, scale);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P1);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(message, gpu_result), true);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        double number = static_cast<double>(dis(gen));
        encoder.encode(P2, number, scale);

        std::vector<double> gpu_result2;
        encoder.decode(gpu_result2, P2);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_equal(number, gpu_result2[0]), true);
    }

    cudaDeviceSynchronize();

    {
        size_t poly_modulus_degree = 32768;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_coeff_modulus_bit_sizes({59, 40, 40, 40, 40, 40, 40, 40,
                                              40, 40, 40, 40, 40, 40, 40, 40,
                                              40, 40, 40},
                                             {59});
        context->generate();

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        const int row_size = poly_modulus_degree / 2;
        std::vector<double> message(row_size, 0);
        for (int i = 0; i < row_size; i++)
        {
            message[i] = dis(gen);
        }

        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        double scale = pow(2.0, 40);
        encoder.encode(P1, message, scale);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P1);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(message, gpu_result), true);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        double number = static_cast<double>(dis(gen));
        encoder.encode(P2, number, scale);

        std::vector<double> gpu_result2;
        encoder.decode(gpu_result2, P2);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_equal(number, gpu_result2[0]), true);
    }

    cudaDeviceSynchronize();

    {
        size_t poly_modulus_degree = 65536;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_coeff_modulus_bit_sizes(
            {59, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
             45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
             45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
            {59});
        context->generate();

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        const int row_size = poly_modulus_degree / 2;
        std::vector<double> message(row_size, 0);
        for (int i = 0; i < row_size; i++)
        {
            message[i] = dis(gen);
        }

        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        double scale = pow(2.0, 45);
        encoder.encode(P1, message, scale);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P1);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(message, gpu_result), true);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        double number = static_cast<double>(dis(gen));
        encoder.encode(P2, number, scale);

        std::vector<double> gpu_result2;
        encoder.decode(gpu_result2, P2);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_equal(number, gpu_result2[0]), true);
    }

    cudaDeviceSynchronize();
}

/**
 * @company CipherFlow
 */
TEST(HEonGPU, CKKS_Sparse_Encoding)
{

    // gap=2: slot_count = N/4
    {
        size_t poly_modulus_degree = 8192;
        int slot_count = poly_modulus_degree / 4;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_slot_count(slot_count);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::Relinkey<heongpu::Scheme::CKKS> relin_key(context);
        keygen.generate_relin_key(relin_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] * message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context), P2(context);
        encoder.encode(P1, message1, scale);
        encoder.encode(P2, message2, scale);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context), C2(context);
        encryptor.encrypt(C1, P1);
        encryptor.encrypt(C2, P2);

        operators.multiply_inplace(C1, C2);
        operators.relinearize_inplace(C1, relin_key);
        operators.rescale_inplace(C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();

    // gap=4: slot_count = N/8, need slot_count >= 2048 so N >= 16384
    {
        size_t poly_modulus_degree = 16384;
        int slot_count = poly_modulus_degree / 8;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_slot_count(slot_count);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::Relinkey<heongpu::Scheme::CKKS> relin_key(context);
        keygen.generate_relin_key(relin_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] * message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context), P2(context);
        encoder.encode(P1, message1, scale);
        encoder.encode(P2, message2, scale);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context), C2(context);
        encryptor.encrypt(C1, P1);
        encryptor.encrypt(C2, P2);

        operators.multiply_inplace(C1, C2);
        operators.relinearize_inplace(C1, relin_key);
        operators.rescale_inplace(C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();
}

/**
 * @company CipherFlow
 */
TEST(HEonGPU, CKKS_RingT_Encoding_FullPacking)
{

    {
        size_t poly_modulus_degree = 8192;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::Relinkey<heongpu::Scheme::CKKS> relin_key(context);
        keygen.generate_relin_key(relin_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        const int slot_count = poly_modulus_degree / 2;
        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] * message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1_ringt(context), P2_ringt(context);
        encoder.encode_ringt(P1_ringt, message1, scale);
        encoder.encode_ringt(P2_ringt, message2, scale);

        int level = 3;
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context), P2(context);
        encoder.ringt_to_pt(P1_ringt, P1, level);
        encoder.ringt_to_pt(P2_ringt, P2, level);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context), C2(context);
        encryptor.encrypt(C1, P1);
        encryptor.encrypt(C2, P2);

        operators.multiply_inplace(C1, C2);
        operators.relinearize_inplace(C1, relin_key);
        operators.rescale_inplace(C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();
}

/**
 * @company CipherFlow
 */
TEST(HEonGPU, CKKS_RingT_Encoding_SparsePacking)
{

    // gap=2: slot_count = N/4
    {
        size_t poly_modulus_degree = 8192;
        int slot_count = poly_modulus_degree / 4;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_slot_count(slot_count);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::Relinkey<heongpu::Scheme::CKKS> relin_key(context);
        keygen.generate_relin_key(relin_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] * message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1_ringt(context), P2_ringt(context);
        encoder.encode_ringt(P1_ringt, message1, scale);
        encoder.encode_ringt(P2_ringt, message2, scale);

        int level = 3;
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context), P2(context);
        encoder.ringt_to_pt(P1_ringt, P1, level);
        encoder.ringt_to_pt(P2_ringt, P2, level);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context), C2(context);
        encryptor.encrypt(C1, P1);
        encryptor.encrypt(C2, P2);

        operators.multiply_inplace(C1, C2);
        operators.relinearize_inplace(C1, relin_key);
        operators.rescale_inplace(C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();

    // gap=4: slot_count = N/8
    {
        size_t poly_modulus_degree = 16384;
        int slot_count = poly_modulus_degree / 8;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_slot_count(slot_count);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::Relinkey<heongpu::Scheme::CKKS> relin_key(context);
        keygen.generate_relin_key(relin_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] * message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1_ringt(context), P2_ringt(context);
        encoder.encode_ringt(P1_ringt, message1, scale);
        encoder.encode_ringt(P2_ringt, message2, scale);

        int level = 5;
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context), P2(context);
        encoder.ringt_to_pt(P1_ringt, P1, level);
        encoder.ringt_to_pt(P2_ringt, P2, level);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context), C2(context);
        encryptor.encrypt(C1, P1);
        encryptor.encrypt(C2, P2);

        operators.multiply_inplace(C1, C2);
        operators.relinearize_inplace(C1, relin_key);
        operators.rescale_inplace(C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();
}

/**
 * @company CipherFlow
 */
TEST(HEonGPU, CKKS_RingT_Encoding_SparsePacking_AddPlain)
{

    // gap=2: slot_count = N/4
    {
        size_t poly_modulus_degree = 8192;
        int slot_count = poly_modulus_degree / 4;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_slot_count(slot_count);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] + message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        encoder.encode(P1, message1, scale);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        encoder.encode_ringt(P2, message2, scale);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context);
        encryptor.encrypt(C1, P1);

        operators.add_plain(C1, P2, C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();

    // gap=4: slot_count = N/8
    {
        size_t poly_modulus_degree = 16384;
        int slot_count = poly_modulus_degree / 8;
        heongpu::HEContext<heongpu::Scheme::CKKS> context =
            heongpu::GenHEContext<heongpu::Scheme::CKKS>(
                heongpu::sec_level_type::none);
        context->set_poly_modulus_degree(poly_modulus_degree);
        context->set_slot_count(slot_count);
        context->set_coeff_modulus_bit_sizes({40, 30, 30, 30, 30, 30}, {40});
        context->generate();

        heongpu::HEKeyGenerator<heongpu::Scheme::CKKS> keygen(context);
        heongpu::Secretkey<heongpu::Scheme::CKKS> secret_key(context);
        keygen.generate_secret_key(secret_key);

        heongpu::Publickey<heongpu::Scheme::CKKS> public_key(context);
        keygen.generate_public_key(public_key, secret_key);

        heongpu::HEEncoder<heongpu::Scheme::CKKS> encoder(context);
        heongpu::HEEncryptor<heongpu::Scheme::CKKS> encryptor(context, public_key);
        heongpu::HEDecryptor<heongpu::Scheme::CKKS> decryptor(context, secret_key);
        heongpu::HEArithmeticOperator<heongpu::Scheme::CKKS> operators(context, encoder);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<double> message1(slot_count), message2(slot_count);
        for (int i = 0; i < slot_count; i++)
        {
            message1[i] = dis(gen);
            message2[i] = dis(gen);
        }
        std::vector<double> expected(slot_count);
        for (int i = 0; i < slot_count; i++)
            expected[i] = message1[i] + message2[i];

        double scale = pow(2.0, 30);
        heongpu::Plaintext<heongpu::Scheme::CKKS> P1(context);
        encoder.encode(P1, message1, scale);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P2(context);
        encoder.encode_ringt(P2, message2, scale);

        heongpu::Ciphertext<heongpu::Scheme::CKKS> C1(context);
        encryptor.encrypt(C1, P1);

        operators.add_plain(C1, P2, C1);

        heongpu::Plaintext<heongpu::Scheme::CKKS> P3(context);
        decryptor.decrypt(P3, C1);

        std::vector<double> gpu_result;
        encoder.decode(gpu_result, P3);

        cudaDeviceSynchronize();

        EXPECT_EQ(fix_point_array_check(expected, gpu_result), true);
    }

    cudaDeviceSynchronize();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}