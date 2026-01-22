// Copyright 2024-2025 Alişah Özcan
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Alişah Özcan

#include "bfv/operator.cuh"

namespace heongpu
{
    __host__
    HEOperator<Scheme::BFV>::HEOperator(HEContext<Scheme::BFV>& context,
                                        HEEncoder<Scheme::BFV>& encoder)
    {
        if (!context.context_generated_)
        {
            throw std::invalid_argument("HEContext is not generated!");
        }

        scheme_ = context.scheme_;

        n = context.n;

        n_power = context.n_power;

        Q_prime_size_ = context.Q_prime_size;
        Q_size_ = context.Q_size;
        P_size_ = context.P_size;

        bsk_mod_count_ = context.bsk_modulus;

        modulus_ = context.modulus_;

        ntt_table_ = context.ntt_table_;

        intt_table_ = context.intt_table_;

        n_inverse_ = context.n_inverse_;

        last_q_modinv_ = context.last_q_modinv_;

        base_Bsk_ = context.base_Bsk_;

        bsk_ntt_tables_ = context.bsk_ntt_tables_;

        bsk_intt_tables_ = context.bsk_intt_tables_;

        bsk_n_inverse_ = context.bsk_n_inverse_;

        m_tilde_ = context.m_tilde_;

        base_change_matrix_Bsk_ = context.base_change_matrix_Bsk_;

        inv_punctured_prod_mod_base_array_ =
            context.inv_punctured_prod_mod_base_array_;

        base_change_matrix_m_tilde_ = context.base_change_matrix_m_tilde_;

        inv_prod_q_mod_m_tilde_ = context.inv_prod_q_mod_m_tilde_;

        inv_m_tilde_mod_Bsk_ = context.inv_m_tilde_mod_Bsk_;

        prod_q_mod_Bsk_ = context.prod_q_mod_Bsk_;

        inv_prod_q_mod_Bsk_ = context.inv_prod_q_mod_Bsk_;

        plain_modulus_ = context.plain_modulus_;

        base_change_matrix_q_ = context.base_change_matrix_q_;

        base_change_matrix_msk_ = context.base_change_matrix_msk_;

        inv_punctured_prod_mod_B_array_ =
            context.inv_punctured_prod_mod_B_array_;

        inv_prod_B_mod_m_sk_ = context.inv_prod_B_mod_m_sk_;

        prod_B_mod_q_ = context.prod_B_mod_q_;

        q_Bsk_merge_modulus_ = context.q_Bsk_merge_modulus_;

        q_Bsk_merge_ntt_tables_ = context.q_Bsk_merge_ntt_tables_;

        q_Bsk_merge_intt_tables_ = context.q_Bsk_merge_intt_tables_;

        q_Bsk_n_inverse_ = context.q_Bsk_n_inverse_;

        half_p_ = context.half_p_;

        half_mod_ = context.half_mod_;

        upper_threshold_ = context.upper_threshold_;

        upper_halfincrement_ = context.upper_halfincrement_;

        Q_mod_t_ = context.Q_mod_t_;

        coeeff_div_plainmod_ = context.coeeff_div_plainmod_;

        //////

        d = context.d;
        d_tilda = context.d_tilda;
        r_prime = context.r_prime;

        // @company CipherFlow begin ---
        l_leveled_ = context.l_leveled;
        l_tilda_leveled_ = context.l_tilda_leveled;
        d_leveled_ = context.d_leveled;
        // @company CipherFlow end ---
        
        B_prime_ = context.B_prime_;
        B_prime_ntt_tables_ = context.B_prime_ntt_tables_;
        B_prime_intt_tables_ = context.B_prime_intt_tables_;
        B_prime_n_inverse_ = context.B_prime_n_inverse_;

        base_change_matrix_D_to_B_ = context.base_change_matrix_D_to_B_;
        base_change_matrix_B_to_D_ = context.base_change_matrix_B_to_D_;
        Mi_inv_D_to_B_ = context.Mi_inv_D_to_B_;
        Mi_inv_B_to_D_ = context.Mi_inv_B_to_D_;
        prod_D_to_B_ = context.prod_D_to_B_;
        prod_B_to_D_ = context.prod_B_to_D_;

        base_change_matrix_D_to_Q_tilda_ =
            context.base_change_matrix_D_to_Q_tilda_;
        Mi_inv_D_to_Q_tilda_ = context.Mi_inv_D_to_Q_tilda_;
        prod_D_to_Q_tilda_ = context.prod_D_to_Q_tilda_;     

        I_j_ = context.I_j_;
        I_location_ = context.I_location_;
        Sk_pair_ = context.Sk_pair_;

        // @company CipherFlow begin ---
        base_change_matrix_D_to_Qtilda_leveled_ =
            context.base_change_matrix_D_to_Qtilda_leveled;
        Mi_inv_D_to_Qtilda_leveled_ = context.Mi_inv_D_to_Qtilda_leveled;
        prod_D_to_Qtilda_leveled_ = context.prod_D_to_Qtilda_leveled;

        I_j_leveled_ = context.I_j_leveled;
        I_location_leveled_ = context.I_location_leveled;
        Sk_pair_leveled_ = context.Sk_pair_leveled;

        prime_location_leveled_ = context.prime_location_leveled;

        rescaled_last_q_modinv_ = context.rescaled_last_q_modinv_;
        rescaled_half_ = context.rescaled_half_;
        rescaled_half_mod_ = context.rescaled_half_mod_;
        // @company CipherFlow end ---

        prime_vector_ = context.prime_vector_;

        std::vector<int> prime_loc;
        std::vector<int> input_loc;

        int counter = Q_size_;
        for (int i = 0; i < Q_size_; i++) // @company CipherFlow
        {
            for (int j = 0; j < counter; j++)
            {
                prime_loc.push_back(j);
            }
            counter--;
            for (int j = 0; j < P_size_; j++)
            {
                prime_loc.push_back(Q_size_ + j);
            }
        }

        counter = Q_prime_size_;
        for (int i = 0; i < Q_prime_size_ - 1; i++)
        {
            int sum = counter - 1;
            for (int j = 0; j < 2; j++)
            {
                input_loc.push_back(sum);
                sum += counter;
            }
            counter--;
        }

        new_prime_locations_ = DeviceVector<int>(prime_loc);
        new_input_locations_ = DeviceVector<int>(input_loc);
        new_prime_locations = new_prime_locations_.data();
        new_input_locations = new_input_locations_.data();

        // @company CipherFlow begin ---
        std::vector<int> prime_bsk_loc;
        counter = Q_size_;
        for (int i = 0; i < Q_size_; i++)
        {
            for (int j = 0; j < counter; j++)
            {
                prime_bsk_loc.push_back(j);
            }
            counter--;
            for (int j = 0; j < bsk_mod_count_[i]; j++)
            {
                prime_bsk_loc.push_back(Q_size_ + j);
            }
        }

        new_prime_bsk_locations_ = DeviceVector<int>(prime_bsk_loc);
        new_prime_bsk_locations = new_prime_bsk_locations_.data();
        // @company CipherFlow end ---

        // Encode params
        slot_count_ = encoder.slot_count_;
        plain_modulus_pointer_ = encoder.plain_modulus_;
        n_plain_inverse_ = encoder.n_plain_inverse_;
        plain_intt_tables_ = encoder.plain_intt_tables_;
        encoding_location_ = encoder.encoding_location_;
    }

    __host__ void HEOperator<Scheme::BFV>::add(Ciphertext<Scheme::BFV>& input1,
                                               Ciphertext<Scheme::BFV>& input2,
                                               Ciphertext<Scheme::BFV>& output,
                                               const ExecutionOptions& options)
    {
        if (input1.depth_ != input2.depth_) // @company CipherFlow
        {
            throw std::logic_error("Ciphertexts leveled are not equal");
        }

        if (input1.relinearization_required_ !=
            input2.relinearization_required_)
        {
            throw std::invalid_argument("Ciphertexts can not be added because "
                                        "ciphertext sizes have to be equal!");
        }

        if (input1.in_ntt_domain_ != input2.in_ntt_domain_)
        {
            throw std::invalid_argument(
                "Both Ciphertexts should be in same domain");
        }

        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        if (input1.memory_size() < (cipher_size * n * current_decomp_count) ||
            input2.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        input_storage_manager(
            input1,
            [&](Ciphertext<Scheme::BFV>& input1_)
            {
                input_storage_manager(
                    input2,
                    [&](Ciphertext<Scheme::BFV>& input2_)
                    {
                        output_storage_manager(
                            output,
                            [&](Ciphertext<Scheme::BFV>& output_)
                            {
                                DeviceVector<Data64> output_memory(
                                    (cipher_size * n * current_decomp_count), // @company CipherFlow
                                    options.stream_);

                                addition<<<dim3((n >> 8), current_decomp_count, cipher_size), // @company CipherFlow
                                           256, 0, options.stream_>>>(
                                    input1_.data(), input2_.data(),
                                    output_memory.data(), modulus_->data(),
                                    n_power);
                                HEONGPU_CUDA_CHECK(cudaGetLastError());

                                output_.scheme_ = scheme_;
                                output_.ring_size_ = n;
                                output_.coeff_modulus_count_ = Q_size_;
                                output_.cipher_size_ = cipher_size;
                                output_.depth_ = input1_.depth_; // @company CipherFlow
                                output_.in_ntt_domain_ = input1_.in_ntt_domain_;
                                output_.rescale_required_ =
                                    (input1_.rescale_required_ ||
                                     input2_.rescale_required_); // @company CipherFlow
                                output_.relinearization_required_ =
                                    input1_.relinearization_required_;
                                output_.ciphertext_generated_ = true;

                                output_.memory_set(std::move(output_memory));
                            },
                            options);
                    },
                    options, (&input2 == &output));
            },
            options, (&input1 == &output));
    }

    __host__ void HEOperator<Scheme::BFV>::sub(Ciphertext<Scheme::BFV>& input1,
                                               Ciphertext<Scheme::BFV>& input2,
                                               Ciphertext<Scheme::BFV>& output,
                                               const ExecutionOptions& options)
    {
        if (input1.depth_ != input2.depth_) // @company CipherFlow
        {
            throw std::logic_error("Ciphertexts leveled are not equal");
        }

        if (input1.relinearization_required_ !=
            input2.relinearization_required_)
        {
            throw std::invalid_argument("Ciphertexts can not be added because "
                                        "ciphertext sizes have to be equal!");
        }

        if (input1.in_ntt_domain_ != input2.in_ntt_domain_)
        {
            throw std::invalid_argument(
                "Both Ciphertexts should be in same domain");
        }

        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        if (input1.memory_size() < (cipher_size * n * current_decomp_count) || 
            input2.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        input_storage_manager(
            input1,
            [&](Ciphertext<Scheme::BFV>& input1_)
            {
                input_storage_manager(
                    input2,
                    [&](Ciphertext<Scheme::BFV>& input2_)
                    {
                        output_storage_manager(
                            output,
                            [&](Ciphertext<Scheme::BFV>& output_)
                            {
                                DeviceVector<Data64> output_memory(
                                    (cipher_size * n * current_decomp_count), // @company CipherFlow
                                    options.stream_);

                                substraction<<<dim3((n >> 8), current_decomp_count, // @company CipherFlow
                                                    cipher_size),
                                               256, 0, options.stream_>>>(
                                    input1_.data(), input2_.data(),
                                    output_memory.data(), modulus_->data(),
                                    n_power);
                                HEONGPU_CUDA_CHECK(cudaGetLastError());

                                output_.scheme_ = scheme_;
                                output_.ring_size_ = n;
                                output_.coeff_modulus_count_ = Q_size_;
                                output_.cipher_size_ = cipher_size;
                                output_.depth_ = input1_.depth_; // @company CipherFlow
                                output_.in_ntt_domain_ = input1_.in_ntt_domain_;
                                output_.rescale_required_ =
                                    (input1_.rescale_required_ ||
                                     input2_.rescale_required_); // @company CipherFlow
                                output_.relinearization_required_ =
                                    input1_.relinearization_required_;
                                output_.ciphertext_generated_ = true;

                                output_.memory_set(std::move(output_memory));
                            },
                            options);
                    },
                    options, (&input2 == &output));
            },
            options, (&input1 == &output));
    }

    __host__ void
    HEOperator<Scheme::BFV>::negate(Ciphertext<Scheme::BFV>& input1,
                                    Ciphertext<Scheme::BFV>& output,
                                    const ExecutionOptions& options)
    {
        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        if (input1.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        input_storage_manager(
            input1,
            [&](Ciphertext<Scheme::BFV>& input1_)
            {
                output_storage_manager(
                    output,
                    [&](Ciphertext<Scheme::BFV>& output_)
                    {
                        DeviceVector<Data64> output_memory(
                            (cipher_size * n * current_decomp_count), options.stream_); // @company CipherFlow

                        negation<<<dim3((n >> 8), current_decomp_count, cipher_size), 256, 0, // @company CipherFlow
                                   options.stream_>>>(
                            input1_.data(), output_memory.data(),
                            modulus_->data(), n_power);
                        HEONGPU_CUDA_CHECK(cudaGetLastError());

                        output_.scheme_ = scheme_;
                        output_.ring_size_ = n;
                        output_.coeff_modulus_count_ = Q_size_;
                        output_.cipher_size_ = cipher_size;
                        output_.depth_ = input1_.depth_; // @company CipherFlow
                        output_.in_ntt_domain_ = input1_.in_ntt_domain_;
                        output_.rescale_required_ = input1_.rescale_required_; // @company CipherFlow
                        output_.relinearization_required_ =
                            input1_.relinearization_required_;
                        output_.ciphertext_generated_ = true;

                        output_.memory_set(std::move(output_memory));
                    },
                    options);
            },
            options, (&input1 == &output));
    }

    __host__ void HEOperator<Scheme::BFV>::add_plain_bfv(
        Ciphertext<Scheme::BFV>& input1, Plaintext<Scheme::BFV>& input2,
        Ciphertext<Scheme::BFV>& output, const cudaStream_t stream)
    {
        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        if (input1.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        if (input2.size() < n)
        {
            throw std::invalid_argument("Invalid Plaintext size!");
        }

        DeviceVector<Data64> output_memory((cipher_size * n * current_decomp_count), stream); // @company CipherFlow

        if (input2.is_ringt_) { // @company CipherFlow
            // @company CipherFlow begin ---
            int counter = Q_size_;
            int location = 0;
            for (int i = 0; i < input1.depth_; i++)
            {
                location += counter;
                counter--;
            }
            // @company CipherFlow end ---

            addition_plain_bfv_poly<<<dim3((n >> 8), current_decomp_count, cipher_size), 256, 0, 
                                    stream>>>(
                input1.data(), input2.data(), output_memory.data(),
                modulus_->data(), plain_modulus_, Q_mod_t_->data() + location, upper_threshold_,
                coeeff_div_plainmod_->data() + location, n_power); // @company CipherFlow
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        } else { // @company CipherFlow
            // @company CipherFlow begin ---
            addition<<<dim3((n >> 8), current_decomp_count, 1), 
                        256, 0, stream>>>(
                input1.data(), input2.data(),
                output_memory.data(), modulus_->data(),
                n_power);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
            
            global_memory_replace_kernel<<<dim3((n >> 8), current_decomp_count, 1), 256, 0,
                    stream>>>(input1.data() + n * current_decomp_count, output_memory.data() + n * current_decomp_count,
                                n_power); 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
            // @company CipherFlow end ---
        }

        output.cipher_size_ = cipher_size;

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::add_plain_bfv_inplace(
        Ciphertext<Scheme::BFV>& input1, Plaintext<Scheme::BFV>& input2,
        const cudaStream_t stream)
    {
        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        if (input1.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        if (input2.size() < n)
        {
            throw std::invalid_argument("Invalid Plaintext size!");
        }

        // @company CipherFlow begin ---
        int counter = Q_size_;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        addition_plain_bfv_poly_inplace<<<dim3((n >> 8), current_decomp_count, 1), 256, 0,
                                          stream>>>(
            input1.data(), input2.data(), input1.data(), modulus_->data(),
            plain_modulus_, Q_mod_t_->data() + location, upper_threshold_,
            coeeff_div_plainmod_->data() + location, n_power); // @company CipherFlow
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }

    __host__ void HEOperator<Scheme::BFV>::sub_plain_bfv(
        Ciphertext<Scheme::BFV>& input1, Plaintext<Scheme::BFV>& input2,
        Ciphertext<Scheme::BFV>& output, const cudaStream_t stream)
    {
        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        if (input1.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        if (input2.size() < n)
        {
            throw std::invalid_argument("Invalid Plaintext size!");
        }

        DeviceVector<Data64> output_memory((cipher_size * n * current_decomp_count), stream); // @company CipherFlow

        if (input2.is_ringt_) { // @company CipherFlow
            // @company CipherFlow begin ---
            int counter = Q_size_;
            int location = 0;
            for (int i = 0; i < input1.depth_; i++)
            {
                location += counter;
                counter--;
            }
            // @company CipherFlow end ---
            
            substraction_plain_bfv_poly<<<dim3((n >> 8), current_decomp_count, cipher_size), 256,
                                        0, stream>>>(
                input1.data(), input2.data(), output_memory.data(),
                modulus_->data(), plain_modulus_, Q_mod_t_->data() + location, upper_threshold_,
                coeeff_div_plainmod_->data() + location, n_power); // @company CipherFlow
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        } else {
            // @company CipherFlow begin ---
            substraction<<<dim3((n >> 8), current_decomp_count, 1), 
                        256, 0, stream>>>(
                input1.data(), input2.data(),
                output_memory.data(), modulus_->data(),
                n_power);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
            
            global_memory_replace_kernel<<<dim3((n >> 8), current_decomp_count, 1), 256, 0,
                    stream>>>(input1.data() + n * current_decomp_count, output_memory.data() + n * current_decomp_count,
                                n_power); 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
            // @company CipherFlow end ---
        }

        output.cipher_size_ = cipher_size;

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::sub_plain_bfv_inplace(
        Ciphertext<Scheme::BFV>& input1, Plaintext<Scheme::BFV>& input2,
        const cudaStream_t stream)
    {
        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        int cipher_size = input1.relinearization_required_ ? 3 : 2;

        if (input1.memory_size() < (cipher_size * n * current_decomp_count)) // @company CipherFlow
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        if (input2.size() < n)
        {
            throw std::invalid_argument("Invalid Plaintext size!");
        }

        // @company CipherFlow begin ---
        int counter = Q_size_;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        substraction_plain_bfv_poly_inplace<<<dim3((n >> 8), current_decomp_count, 1), 256,
                                              0, stream>>>(
            input1.data(), input2.data(), input1.data(), modulus_->data(),
            plain_modulus_, Q_mod_t_->data() + location, upper_threshold_,
            coeeff_div_plainmod_->data() + location, n_power); // @company CipherFlow
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }


    __host__ void HEOperator<Scheme::BFV>::multiply_bfv(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& input2,
        Ciphertext<Scheme::BFV>& output, const cudaStream_t stream)
    {   
        // @company CipherFlow
        if (input1.depth_ != input2.depth_)
        {
            throw std::logic_error("Ciphertexts leveled are not equal");
        }

        if ((input1.in_ntt_domain_ != false) ||
            (input2.in_ntt_domain_ != false))
        {
            throw std::invalid_argument("Ciphertexts should be in same domain");
        }

        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow

        if (input1.memory_size() < (2 * n * current_decomp_count) ||
            input2.memory_size() < (2 * n * current_decomp_count))
        {
            throw std::invalid_argument("Invalid Ciphertexts size!");
        }

        DeviceVector<Data64> output_memory((3 * n * current_decomp_count), stream); // @company CipherFlow

        DeviceVector<Data64> temp_mul((4 * n * (bsk_mod_count_[0] + Q_size_)) +
                                          (3 * n * (bsk_mod_count_[0] + Q_size_)),
                                      stream);  // @company CipherFlow
        Data64* temp1_mul = temp_mul.data();
        Data64* temp2_mul = temp1_mul + (4 * n * (bsk_mod_count_[0] + Q_size_));  // @company CipherFlow
       
        // @company CipherFlow begin ---
        int counter = Q_size_;
        int location = 0;
        int location1 = 0;
        int location2 = 0;
        int location3 = 0;
        int location4 = 0;
        int location5 = 0;
        for (int i = 0; i < input1.depth_; i++)
        {   
            location += counter;
            location1 += bsk_mod_count_[i];
            location2 += bsk_mod_count_[i] * counter;
            location3 += bsk_mod_count_[i]-1;
            location4 += (bsk_mod_count_[i]-1) * counter;
            location5 += counter + bsk_mod_count_[i];
            counter--;
        }
        // @company CipherFlow end ---

        fast_convertion<<<dim3((n >> 8), 4, 1), 256, 0, stream>>>(
            input1.data(), input2.data(), temp1_mul, modulus_->data(),
            base_Bsk_->data(), m_tilde_, inv_prod_q_mod_m_tilde_->data() + input1.depth_,
            inv_m_tilde_mod_Bsk_->data(), prod_q_mod_Bsk_->data() + location1,
            base_change_matrix_Bsk_->data() + location2,
            base_change_matrix_m_tilde_->data() + location,
            inv_punctured_prod_mod_base_array_->data() + location, n_power, current_decomp_count,
            bsk_mod_count_[input1.depth_]); // @company CipherFlow
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = q_Bsk_n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp1_mul, q_Bsk_merge_ntt_tables_->data(),
            q_Bsk_merge_modulus_->data(), cfg_ntt,
            ((bsk_mod_count_[input1.depth_] + current_decomp_count) * 4), (bsk_mod_count_[input1.depth_] + current_decomp_count),
            new_prime_bsk_locations + location5); // @company CipherFlow

        cross_multiplication<<<dim3((n >> 8), (bsk_mod_count_[input1.depth_] + current_decomp_count), 1),
                               256, 0, stream>>>(
            temp1_mul, temp1_mul + (((bsk_mod_count_[input1.depth_] + current_decomp_count) * 2) * n),
            temp2_mul, q_Bsk_merge_modulus_->data(), n_power,
            (bsk_mod_count_[input1.depth_] + current_decomp_count),
            current_decomp_count, input1.depth_); // @company CipherFlow
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp2_mul, q_Bsk_merge_intt_tables_->data(),
            q_Bsk_merge_modulus_->data(), cfg_intt,
            (3 * (bsk_mod_count_[input1.depth_] + current_decomp_count)),
            (bsk_mod_count_[input1.depth_] + current_decomp_count),
            new_prime_bsk_locations + location5); // @company CipherFlow

        fast_floor<<<dim3((n >> 8), 3, 1), 256, 0, stream>>>(
            temp2_mul, output_memory.data(), modulus_->data(),
            base_Bsk_->data(), plain_modulus_,
            inv_punctured_prod_mod_base_array_->data() + location,
            base_change_matrix_Bsk_->data() + location2, inv_prod_q_mod_Bsk_->data() + location1,
            inv_punctured_prod_mod_B_array_->data() + location3,
            base_change_matrix_q_->data() + location4, base_change_matrix_msk_->data() + location3,
            inv_prod_B_mod_m_sk_->data() +  input1.depth_, prod_B_mod_q_->data() + location, n_power, current_decomp_count,
            bsk_mod_count_[input1.depth_]); // @company CipherFlow
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::multiply_plain_bfv(
        Ciphertext<Scheme::BFV>& input1, Plaintext<Scheme::BFV>& input2,
        Ciphertext<Scheme::BFV>& output, const cudaStream_t stream)
    {   
        int current_decomp_count = Q_size_ - input1.depth_; // @company CipherFlow
        DeviceVector<Data64> output_memory((2 * n * current_decomp_count), stream); // @company CipherFlow

        if (input1.in_ntt_domain_)
        {
            cipherplain_kernel<<<dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
                input1.data(), input2.data(), output_memory.data(),
                modulus_->data(), n_power, current_decomp_count); // @company CipherFlow
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> temp_plain_mul(n * Q_size_, stream);
            Data64* temp1_plain_mul = temp_plain_mul.data();

            threshold_kernel<<<dim3((n >> 8), current_decomp_count, 1), 256, 0, stream>>>(
                input2.data(), temp1_plain_mul, modulus_->data(),
                upper_halfincrement_->data(), upper_threshold_, n_power,
                current_decomp_count); // @company CipherFlow
            HEONGPU_CUDA_CHECK(cudaGetLastError());

            gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
                .n_power = n_power,
                .ntt_type = gpuntt::FORWARD,
                .ntt_layout = gpuntt::PerPolynomial,                
                .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
                .zero_padding = false,
                .stream = stream};

            gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
                .n_power = n_power,
                .ntt_type = gpuntt::INVERSE,
                .ntt_layout = gpuntt::PerPolynomial,                
                .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
                .zero_padding = false,
                .mod_inverse = n_inverse_->data(),
                .stream = stream};

            gpuntt::GPU_NTT_Inplace(temp1_plain_mul, ntt_table_->data(),
                                    modulus_->data(), cfg_ntt, current_decomp_count,
                                    current_decomp_count); // @company CipherFlow

            gpuntt::GPU_NTT(input1.data(), output_memory.data(),
                            ntt_table_->data(), modulus_->data(), cfg_ntt,
                            2 * current_decomp_count, current_decomp_count); // @company CipherFlow

            cipherplain_kernel<<<dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
                output_memory.data(), temp1_plain_mul, output_memory.data(),
                modulus_->data(), n_power, current_decomp_count); // @company CipherFlow
            HEONGPU_CUDA_CHECK(cudaGetLastError());

            gpuntt::GPU_INTT_Inplace(output_memory.data(), intt_table_->data(),
                                    modulus_->data(), cfg_intt, 2 * current_decomp_count,
                                    current_decomp_count); // @company CipherFlow
        }

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::relinearize_seal_method_inplace(
        Ciphertext<Scheme::BFV>& input1, Relinkey<Scheme::BFV>& relin_key,
        const cudaStream_t stream)
    {
        DeviceVector<Data64> temp_relin(
            (n * Q_size_ * Q_prime_size_) + (2 * n * Q_prime_size_), stream);

        Data64* temp1_relin = temp_relin.data();
        Data64* temp2_relin = temp1_relin + (n * Q_size_ * Q_prime_size_);

        cipher_broadcast_kernel<<<dim3((n >> 8), Q_size_, 1), 256, 0, stream>>>(
            input1.data() + (Q_size_ << (n_power + 1)), temp1_relin,
            modulus_->data(), n_power, Q_prime_size_);

        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Inplace(temp1_relin, ntt_table_->data(),
                                modulus_->data(), cfg_ntt,
                                Q_size_ * Q_prime_size_, Q_prime_size_);

        int iteration_count_1 = Q_size_ / 4;
        int iteration_count_2 = Q_size_ % 4;
        // TODO: make it efficient
        if (relin_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_relin, relin_key.data(), temp2_relin, modulus_->data(),
                n_power, Q_prime_size_, iteration_count_1, iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(relin_key.host_location_, stream);
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_relin, key_location.data(), temp2_relin, modulus_->data(),
                n_power, Q_prime_size_, iteration_count_1, iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT_Inplace(temp2_relin, intt_table_->data(),
                                modulus_->data(), cfg_intt, 2 * Q_prime_size_,
                                Q_prime_size_);

        divide_round_lastq_kernel<<<dim3((n >> 8), Q_size_, 2), 256, 0,
                                    stream>>>(
            temp2_relin, input1.data(), input1.data(), modulus_->data(),
            half_p_->data(), half_mod_->data(), last_q_modinv_->data(), n_power,
            Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }

    __host__ void
    HEOperator<Scheme::BFV>::relinearize_external_product_method_inplace(
        Ciphertext<Scheme::BFV>& input1, Relinkey<Scheme::BFV>& relin_key,
        const cudaStream_t stream)
    {
        DeviceVector<Data64> temp_relin_new((n * d * r_prime) +
                                                (2 * n * d_tilda * r_prime) +
                                                (2 * n * Q_prime_size_),
                                            stream);
        Data64* temp1_relin_new = temp_relin_new.data();
        Data64* temp2_relin_new = temp1_relin_new + (n * d * r_prime);
        Data64* temp3_relin_new = temp2_relin_new + (2 * n * d_tilda * r_prime);

        base_conversion_DtoB_relin_kernel<<<dim3((n >> 8), d, 1), 256, 0,
                                            stream>>>(
            input1.data() + (Q_size_ << (n_power + 1)), temp1_relin_new,
            modulus_->data(), B_prime_->data(),
            base_change_matrix_D_to_B_->data(), Mi_inv_D_to_B_->data(),
            prod_D_to_B_->data(), I_j_->data(), I_location_->data(), n_power,
            Q_size_, d_tilda, d, r_prime);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Inplace(temp1_relin_new, B_prime_ntt_tables_->data(),
                                B_prime_->data(), cfg_ntt, d * r_prime,
                                r_prime);

        // TODO: make it efficient
        if (relin_key.storage_type_ == storage_type::DEVICE)
        {
            multiply_accumulate_extended_kernel<<<
                dim3((n >> 8), r_prime, d_tilda), 256, 0, stream>>>(
                temp1_relin_new, relin_key.data(), temp2_relin_new,
                B_prime_->data(), n_power, d_tilda, d, r_prime);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(relin_key.host_location_, stream);
            multiply_accumulate_extended_kernel<<<
                dim3((n >> 8), r_prime, d_tilda), 256, 0, stream>>>(
                temp1_relin_new, key_location.data(), temp2_relin_new,
                B_prime_->data(), n_power, d_tilda, d, r_prime);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = B_prime_n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT_Inplace(temp2_relin_new, B_prime_intt_tables_->data(),
                                B_prime_->data(), cfg_intt,
                                2 * r_prime * d_tilda, r_prime);

        base_conversion_BtoD_relin_kernel<<<dim3((n >> 8), d_tilda, 2), 256, 0,
                                            stream>>>(
            temp2_relin_new, temp3_relin_new, modulus_->data(),
            B_prime_->data(), base_change_matrix_B_to_D_->data(),
            Mi_inv_B_to_D_->data(), prod_B_to_D_->data(), I_j_->data(),
            I_location_->data(), n_power, Q_prime_size_, d_tilda, d, r_prime);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        divide_round_lastq_extended_kernel<<<dim3((n >> 8), Q_size_, 2), 256, 0,
                                             stream>>>(
            temp3_relin_new, input1.data(), input1.data(), modulus_->data(),
            half_p_->data(), half_mod_->data(), last_q_modinv_->data(), n_power,
            Q_prime_size_, Q_size_, P_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }

    __host__ void
    HEOperator<Scheme::BFV>::relinearize_external_product_method2_inplace(
        Ciphertext<Scheme::BFV>& input1, Relinkey<Scheme::BFV>& relin_key,
        const cudaStream_t stream)
    {   
        // @company CipherFlow begin ---
        int first_rns_mod_count = Q_prime_size_;
        int current_rns_mod_count = Q_prime_size_ - input1.depth_;

        int first_decomp_count = Q_size_;
        int current_decomp_count = Q_size_ - input1.depth_;
        // @company CipherFlow end ---

        DeviceVector<Data64> temp_relin(
            (n * Q_size_ * Q_prime_size_) + (2 * n * Q_prime_size_), stream);
        Data64* temp1_relin = temp_relin.data();
        Data64* temp2_relin = temp1_relin + (n * Q_size_ * Q_prime_size_);

        // @company CipherFlow begin ---
        int counter = first_rns_mod_count;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        base_conversion_DtoQtilde_relin_leveled_kernel<<<
            dim3((n >> 8), d_leveled_->operator[](input1.depth_), 1), 256, 0,
            stream>>>(
            input1.data() + (current_decomp_count << (n_power + 1)),
            temp1_relin, modulus_->data(),
            base_change_matrix_D_to_Qtilda_leveled_->operator[](input1.depth_)
                .data(),
            Mi_inv_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            prod_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            I_j_leveled_->operator[](input1.depth_).data(),
            I_location_leveled_->operator[](input1.depth_).data(), n_power,
            d_leveled_->operator[](input1.depth_), current_rns_mod_count,
            current_decomp_count, input1.depth_,
            prime_location_leveled_->data() + location); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp1_relin, ntt_table_->data(), modulus_->data(), cfg_ntt,
            d_leveled_->operator[](input1.depth_) * current_rns_mod_count,
            current_rns_mod_count, new_prime_locations + location); // @company CipherFlow 

        // TODO: make it efficient
        int iteration_count_1 = d_leveled_->operator[](input1.depth_) / 4; // @company CipherFlow 
        int iteration_count_2 = d_leveled_->operator[](input1.depth_) % 4; // @company CipherFlow 
        if (relin_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp1_relin, relin_key.data(), temp2_relin, modulus_->data(),
                first_rns_mod_count, current_decomp_count,
                current_rns_mod_count, iteration_count_1, iteration_count_2,
                input1.depth_, n_power); // @company CipherFlow 
            HEONGPU_CUDA_CHECK(cudaGetLastError()); 
        }
        else
        {
            DeviceVector<Data64> key_location(relin_key.host_location_, stream);
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp1_relin, key_location.data(), temp2_relin, modulus_->data(),
                first_rns_mod_count, current_decomp_count,
                current_rns_mod_count, iteration_count_1, iteration_count_2,
                input1.depth_, n_power); // @company CipherFlow 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,            
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp2_relin, intt_table_->data(), modulus_->data(), cfg_intt,
            2 * current_rns_mod_count, current_rns_mod_count,
            new_prime_locations + location); // @company CipherFlow 

        divide_round_lastq_extended_leveled_kernel<<<
            dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
            temp2_relin, temp1_relin, modulus_->data(), half_p_->data(),
            half_mod_->data(), last_q_modinv_->data(), n_power,
            current_rns_mod_count, current_decomp_count, first_rns_mod_count,
            first_decomp_count, P_size_); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError()); 

        addition<<<dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
            temp1_relin, input1.data(), input1.data(), modulus_->data(),
            n_power); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());  // @company CipherFlow 
    }

    /**
     * @company CipherFlow
     */
    __host__ void
    HEOperator<Scheme::BFV>::relinearize_external_product_method2(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Relinkey<Scheme::BFV>& relin_key, const cudaStream_t stream)
    {   
        int first_rns_mod_count = Q_prime_size_;
        int current_rns_mod_count = Q_prime_size_ - input1.depth_;

        int first_decomp_count = Q_size_;
        int current_decomp_count = Q_size_ - input1.depth_;

        DeviceVector<Data64> output_memory((2 * n * current_decomp_count), stream);

        DeviceVector<Data64> temp_relin(
            (n * Q_size_ * Q_prime_size_) + (2 * n * Q_prime_size_), stream);
        Data64* temp1_relin = temp_relin.data();
        Data64* temp2_relin = temp1_relin + (n * Q_size_ * Q_prime_size_);

        int counter = first_rns_mod_count;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }

        base_conversion_DtoQtilde_relin_leveled_kernel<<<
            dim3((n >> 8), d_leveled_->operator[](input1.depth_), 1), 256, 0,
            stream>>>(
            input1.data() + (current_decomp_count << (n_power + 1)),
            temp1_relin, modulus_->data(),
            base_change_matrix_D_to_Qtilda_leveled_->operator[](input1.depth_)
                .data(),
            Mi_inv_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            prod_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            I_j_leveled_->operator[](input1.depth_).data(),
            I_location_leveled_->operator[](input1.depth_).data(), n_power,
            d_leveled_->operator[](input1.depth_), current_rns_mod_count,
            current_decomp_count, input1.depth_,
            prime_location_leveled_->data() + location); 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp1_relin, ntt_table_->data(), modulus_->data(), cfg_ntt,
            d_leveled_->operator[](input1.depth_) * current_rns_mod_count,
            current_rns_mod_count, new_prime_locations + location); 

        // TODO: make it efficient
        int iteration_count_1 = d_leveled_->operator[](input1.depth_) / 4; 
        int iteration_count_2 = d_leveled_->operator[](input1.depth_) % 4; 
        if (relin_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp1_relin, relin_key.data(), temp2_relin, modulus_->data(),
                first_rns_mod_count, current_decomp_count,
                current_rns_mod_count, iteration_count_1, iteration_count_2,
                input1.depth_, n_power); 
            HEONGPU_CUDA_CHECK(cudaGetLastError()); 
        }
        else
        {
            DeviceVector<Data64> key_location(relin_key.host_location_, stream);
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp1_relin, key_location.data(), temp2_relin, modulus_->data(),
                first_rns_mod_count, current_decomp_count,
                current_rns_mod_count, iteration_count_1, iteration_count_2,
                input1.depth_, n_power); 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp2_relin, intt_table_->data(), modulus_->data(), cfg_intt,
            2 * current_rns_mod_count, current_rns_mod_count,
            new_prime_locations + location);

        divide_round_lastq_extended_leveled_kernel<<<
            dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
            temp2_relin, temp1_relin, modulus_->data(), half_p_->data(),
            half_mod_->data(), last_q_modinv_->data(), n_power,
            current_rns_mod_count, current_decomp_count, first_rns_mod_count,
            first_decomp_count, P_size_); 
        HEONGPU_CUDA_CHECK(cudaGetLastError()); 

        addition<<<dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
            temp1_relin, input1.data(), output_memory.data(), modulus_->data(),
            n_power); 
        HEONGPU_CUDA_CHECK(cudaGetLastError());  

        output.memory_set(std::move(output_memory));
    }

    /**
     * @company CipherFlow
     */
    __host__ void HEOperator<Scheme::BFV>::rescale_inplace_bfv_leveled(
        Ciphertext<Scheme::BFV>& input1, const cudaStream_t stream)
    {
        int first_decomp_count = Q_size_;
        int current_decomp_count = Q_size_ - input1.depth_;

        int counter = first_decomp_count - 1;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        
        DeviceVector<Data64> temp_rescale(
            (2 * n * Q_prime_size_) + (2 * n * Q_prime_size_), stream);
        Data64* temp1_rescale = temp_rescale.data();
        Data64* temp2_rescale = temp1_rescale + (2 * n * Q_prime_size_);

        divide_round_lastq_leveled_stage_one_kernel<<<dim3((n >> 8), 2, 1), 256,
                                                      0, stream>>>(
            input1.data(), temp1_rescale, modulus_->data(),
            rescaled_half_->data() + input1.depth_,
            rescaled_half_mod_->data() + location, n_power,
            current_decomp_count - 1, current_decomp_count - 1);
        
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        move_cipher_leveled_kernel<<<
            dim3((n >> 8), current_decomp_count - 1, 2), 256, 0, stream>>>(
            input1.data(), temp2_rescale, n_power, current_decomp_count - 1);
        
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        divide_round_lastq_rescale_kernel<<<
            dim3((n >> 8), current_decomp_count - 1, 2), 256, 0, stream>>>(
            temp1_rescale, temp2_rescale, input1.data(), modulus_->data(),
            rescaled_last_q_modinv_->data() + location, n_power,
            current_decomp_count - 1);

        HEONGPU_CUDA_CHECK(cudaGetLastError());

        input1.depth_++;
    }

    /**
     * @company CipherFlow
     */
    __host__ void HEOperator<Scheme::BFV>::rescale_bfv_leveled(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output, const cudaStream_t stream)
    {
        int first_decomp_count = Q_size_;
        int current_decomp_count = Q_size_ - input1.depth_;

        int counter = first_decomp_count - 1;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }

        DeviceVector<Data64> output_memory((2 * n * current_decomp_count), stream);
        
        DeviceVector<Data64> temp_rescale(
            (2 * n * Q_prime_size_) + (2 * n * Q_prime_size_), stream);
        Data64* temp1_rescale = temp_rescale.data();
        Data64* temp2_rescale = temp1_rescale + (2 * n * Q_prime_size_);

        divide_round_lastq_leveled_stage_one_kernel<<<dim3((n >> 8), 2, 1), 256,
                                                      0, stream>>>(
            input1.data(), temp1_rescale, modulus_->data(),
            rescaled_half_->data() + input1.depth_,
            rescaled_half_mod_->data() + location, n_power,
            current_decomp_count - 1, current_decomp_count - 1);
        
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        move_cipher_leveled_kernel<<<
            dim3((n >> 8), current_decomp_count - 1, 2), 256, 0, stream>>>(
            input1.data(), temp2_rescale, n_power, current_decomp_count - 1);
        
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        divide_round_lastq_rescale_kernel<<<
            dim3((n >> 8), current_decomp_count - 1, 2), 256, 0, stream>>>(
            temp1_rescale, temp2_rescale, output_memory.data(), modulus_->data(),
            rescaled_last_q_modinv_->data() + location, n_power,
            current_decomp_count - 1);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));

    }

    __host__ void HEOperator<Scheme::BFV>::rotate_method_I(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Galoiskey<Scheme::BFV>& galois_key, int shift,
        const cudaStream_t stream)
    {
        int galoiselt = steps_to_galois_elt(shift, n, galois_key.group_order_);
        bool key_exist = (galois_key.storage_type_ == storage_type::DEVICE)
                             ? (galois_key.device_location_.find(galoiselt) !=
                                galois_key.device_location_.end())
                             : (galois_key.host_location_.find(galoiselt) !=
                                galois_key.host_location_.end());
        if (key_exist)
        {
            apply_galois_method_I(input1, output, galois_key, galoiselt,
                                  stream);
        }
        else
        {
            std::vector<int> required_galoiselt;
            int shift_num = abs(shift);
            int negative = (shift < 0) ? (-1) : 1;
            while (shift_num != 0)
            {
                int power = int(log2(shift_num));
                int power_2 = pow(2, power);
                shift_num = shift_num - power_2;

                int index_in = power_2 * negative;

                if (!(galois_key.galois_elt.find(index_in) !=
                      galois_key.galois_elt.end()))
                {
                    throw std::logic_error("Galois key not present!");
                }
                galoiselt = galois_key.galois_elt[index_in];
                required_galoiselt.push_back(galoiselt);
            }

            Ciphertext<Scheme::BFV>& in_data = input1;
            for (auto& galois_elt : required_galoiselt)
            {
                apply_galois_method_I(in_data, output, galois_key, galois_elt,
                                      stream);
                in_data = output;
            }
        }
    }

    __host__ void HEOperator<Scheme::BFV>::rotate_method_II(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Galoiskey<Scheme::BFV>& galois_key, int shift,
        const cudaStream_t stream)
    {
        int galoiselt = steps_to_galois_elt(shift, n, galois_key.group_order_);
        bool key_exist = (galois_key.storage_type_ == storage_type::DEVICE)
                             ? (galois_key.device_location_.find(galoiselt) !=
                                galois_key.device_location_.end())
                             : (galois_key.host_location_.find(galoiselt) !=
                                galois_key.host_location_.end());
        if (key_exist)
        {
            apply_galois_method_II(input1, output, galois_key, galoiselt,
                                   stream);
        }
        else
        {
            std::vector<int> required_galoiselt;
            int shift_num = abs(shift);
            int negative = (shift < 0) ? (-1) : 1;
            while (shift_num != 0)
            {
                int power = int(log2(shift_num));
                int power_2 = pow(2, power);
                shift_num = shift_num - power_2;

                int index_in = power_2 * negative;

                if (!(galois_key.galois_elt.find(index_in) !=
                      galois_key.galois_elt.end()))
                {
                    throw std::logic_error("Galois key not present!");
                }
                galoiselt = galois_key.galois_elt[index_in];
                required_galoiselt.push_back(galoiselt);
            }

            Ciphertext<Scheme::BFV>& in_data = input1;
            for (auto& galois_elt : required_galoiselt)
            {
                apply_galois_method_II(in_data, output, galois_key, galois_elt,
                                       stream);
                in_data = output;
            }
        }
    }

    __host__ void HEOperator<Scheme::BFV>::apply_galois_method_I(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Galoiskey<Scheme::BFV>& galois_key, int galois_elt,
        const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        DeviceVector<Data64> temp_rotation((2 * n * Q_size_) +
                                               (n * Q_size_ * Q_prime_size_) +
                                               (2 * n * Q_prime_size_),
                                           stream);
        Data64* temp0_rotation = temp_rotation.data();
        Data64* temp1_rotation = temp0_rotation + (2 * n * Q_size_);
        Data64* temp2_rotation = temp1_rotation + (n * Q_size_ * Q_prime_size_);

        bfv_duplicate_kernel<<<dim3((n >> 8), Q_size_, 2), 256, 0, stream>>>(
            input1.data(), temp0_rotation, temp1_rotation, modulus_->data(),
            n_power, Q_prime_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Inplace(temp1_rotation, ntt_table_->data(),
                                modulus_->data(), cfg_ntt,
                                Q_size_ * Q_prime_size_, Q_prime_size_);

        // MultSum
        // TODO: make it efficient
        int iteration_count_1 = Q_size_ / 4;
        int iteration_count_2 = Q_size_ % 4;
        if (galois_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_rotation, galois_key.device_location_[galois_elt].data(),
                temp2_rotation, modulus_->data(), n_power, Q_prime_size_,
                iteration_count_1, iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(
                galois_key.host_location_[galois_elt], stream);
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_rotation, key_location.data(), temp2_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT_Inplace(temp2_rotation, intt_table_->data(),
                                modulus_->data(), cfg_intt, 2 * Q_prime_size_,
                                Q_prime_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        // ModDown + Permute
        divide_round_lastq_permute_bfv_kernel<<<dim3((n >> 8), Q_size_, 2), 256,
                                                0, stream>>>(
            temp2_rotation, temp0_rotation, output_memory.data(),
            modulus_->data(), half_p_->data(), half_mod_->data(),
            last_q_modinv_->data(), galois_elt, n_power, Q_prime_size_, Q_size_,
            P_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::apply_galois_method_II(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Galoiskey<Scheme::BFV>& galois_key, int galois_elt,
        const cudaStream_t stream)
    {
        // @company CipherFlow begin ---
        int first_rns_mod_count = Q_prime_size_; 
        int current_rns_mod_count = Q_prime_size_ - input1.depth_;

        int first_decomp_count = Q_size_;
        int current_decomp_count = Q_size_ - input1.depth_;
        // @company CipherFlow end ---

        DeviceVector<Data64> output_memory((2 * n * current_decomp_count), stream); // @company CipherFlow 

        DeviceVector<Data64> temp_rotation((2 * n * Q_size_) + (n * Q_size_) +
                                               (2 * n * Q_prime_size_ * d_leveled_->operator[](0)) +
                                               (2 * n * Q_prime_size_),
                                           stream); // @company CipherFlow 

        Data64* temp0_rotation = temp_rotation.data();
        Data64* temp1_rotation = temp0_rotation + (2 * n * Q_size_);
        Data64* temp2_rotation = temp1_rotation + (n * Q_size_);
        Data64* temp3_rotation = temp2_rotation + (2 * n * Q_prime_size_ * d_leveled_->operator[](0)); // @company CipherFlow 

        // TODO: make it efficient
        global_memory_replace_kernel<<<dim3((n >> 8), current_decomp_count, 1), 256, 0,
                                       stream>>>(input1.data(), temp0_rotation,
                                                 n_power); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        // @company CipherFlow begin ---
        int counter = first_rns_mod_count;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        base_conversion_DtoQtilde_relin_leveled_kernel<<<
            dim3((n >> 8), d_leveled_->operator[](input1.depth_), 1), 256, 0,
            stream>>>(
            input1.data() + (current_decomp_count << n_power), temp2_rotation,
            modulus_->data(),
            base_change_matrix_D_to_Qtilda_leveled_->operator[](input1.depth_)
                .data(),
            Mi_inv_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            prod_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            I_j_leveled_->operator[](input1.depth_).data(),
            I_location_leveled_->operator[](input1.depth_).data(), n_power,
            d_leveled_->operator[](input1.depth_), current_rns_mod_count,
            current_decomp_count, input1.depth_,
            prime_location_leveled_->data() + location); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp2_rotation, ntt_table_->data(), modulus_->data(), cfg_ntt,
            d_leveled_->operator[](input1.depth_) * current_rns_mod_count,
            current_rns_mod_count, new_prime_locations + location); // @company CipherFlow 

        // MultSum
        // TODO: make it efficient
        int iteration_count_1 = d_leveled_->operator[](input1.depth_) / 4; // @company CipherFlow 
        int iteration_count_2 = d_leveled_->operator[](input1.depth_) % 4; // @company CipherFlow 
        if (galois_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp2_rotation, galois_key.device_location_[galois_elt].data(),
                temp3_rotation, modulus_->data(), first_rns_mod_count,
                current_decomp_count, current_rns_mod_count, iteration_count_1,
                iteration_count_2, input1.depth_, n_power); // @company CipherFlow 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(
                galois_key.host_location_[galois_elt], stream);
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp2_rotation, key_location.data(), temp3_rotation,
                modulus_->data(), first_rns_mod_count, current_decomp_count,
                current_rns_mod_count, iteration_count_1, iteration_count_2,
                input1.depth_, n_power); // @company CipherFlow 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp3_rotation, intt_table_->data(), modulus_->data(), cfg_intt,
            2 * current_rns_mod_count, current_rns_mod_count,
            new_prime_locations + location); // @company CipherFlow 

        // ModDown + Permute
        divide_round_lastq_permute_bfv_kernel<<<
            dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
            temp3_rotation, temp0_rotation, output_memory.data(),
            modulus_->data(), half_p_->data(), half_mod_->data(),
            last_q_modinv_->data(), galois_elt, n_power, current_rns_mod_count,
            current_decomp_count, first_rns_mod_count, first_decomp_count,
            P_size_); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::rotate_columns_method_I(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Galoiskey<Scheme::BFV>& galois_key, const cudaStream_t stream)
    {
        int galoiselt = galois_key.galois_elt_zero;

        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        DeviceVector<Data64> temp_rotation((2 * n * Q_size_) +
                                               (n * Q_size_ * Q_prime_size_) +
                                               (2 * n * Q_prime_size_),
                                           stream);
        Data64* temp0_rotation = temp_rotation.data();
        Data64* temp1_rotation = temp0_rotation + (2 * n * Q_size_);
        Data64* temp2_rotation = temp1_rotation + (n * Q_size_ * Q_prime_size_);

        bfv_duplicate_kernel<<<dim3((n >> 8), Q_size_, 2), 256, 0, stream>>>(
            input1.data(), temp0_rotation, temp1_rotation, modulus_->data(),
            n_power, Q_prime_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Inplace(temp1_rotation, ntt_table_->data(),
                                modulus_->data(), cfg_ntt,
                                Q_size_ * Q_prime_size_, Q_prime_size_);

        // MultSum
        // TODO: make it efficient
        int iteration_count_1 = Q_size_ / 4;
        int iteration_count_2 = Q_size_ % 4;
        if (galois_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_rotation, galois_key.c_data(), temp2_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(galois_key.zero_host_location_,
                                              stream);
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_rotation, key_location.data(), temp2_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT_Inplace(temp2_rotation, intt_table_->data(),
                                modulus_->data(), cfg_intt, 2 * Q_prime_size_,
                                Q_prime_size_);

        // ModDown + Permute
        divide_round_lastq_permute_bfv_kernel<<<dim3((n >> 8), Q_size_, 2), 256,
                                                0, stream>>>(
            temp2_rotation, temp0_rotation, output_memory.data(),
            modulus_->data(), half_p_->data(), half_mod_->data(),
            last_q_modinv_->data(), galoiselt, n_power, Q_prime_size_, Q_size_,
            P_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::rotate_columns_method_II(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Galoiskey<Scheme::BFV>& galois_key, const cudaStream_t stream)
    {
        // @company CipherFlow begin ---
        int first_rns_mod_count = Q_prime_size_; 
        int current_rns_mod_count = Q_prime_size_ - input1.depth_;

        int first_decomp_count = Q_size_;
        int current_decomp_count = Q_size_ - input1.depth_;
        // @company CipherFlow end ---

        int galoiselt = galois_key.galois_elt_zero;

        DeviceVector<Data64> output_memory((2 * n * current_decomp_count), stream); // @company CipherFlow 

        DeviceVector<Data64> temp_rotation((2 * n * Q_size_) + (n * Q_size_) +
                                               (2 * n * Q_prime_size_ * d_leveled_->operator[](0)) +
                                               (2 * n * Q_prime_size_), // @company CipherFlow 
                                           stream);

        Data64* temp0_rotation = temp_rotation.data();
        Data64* temp1_rotation = temp0_rotation + (2 * n * Q_size_);
        Data64* temp2_rotation = temp1_rotation + (n * Q_size_);
        Data64* temp3_rotation = temp2_rotation + (2 * n * Q_prime_size_ * d_leveled_->operator[](0)); // @company CipherFlow 

        // TODO: make it efficient
        global_memory_replace_kernel<<<dim3((n >> 8), current_decomp_count, 1), 256, 0,
                                       stream>>>(input1.data(), temp0_rotation,
                                                 n_power); // @company CipherFlow
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        // @company CipherFlow begin ---
        int counter = first_rns_mod_count;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        base_conversion_DtoQtilde_relin_leveled_kernel<<<
            dim3((n >> 8), d_leveled_->operator[](input1.depth_), 1), 256, 0,
            stream>>>(
            input1.data() + (current_decomp_count << n_power), temp2_rotation,
            modulus_->data(),
            base_change_matrix_D_to_Qtilda_leveled_->operator[](input1.depth_)
                .data(),
            Mi_inv_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            prod_D_to_Qtilda_leveled_->operator[](input1.depth_).data(),
            I_j_leveled_->operator[](input1.depth_).data(),
            I_location_leveled_->operator[](input1.depth_).data(), n_power,
            d_leveled_->operator[](input1.depth_), current_rns_mod_count,
            current_decomp_count, input1.depth_,
            prime_location_leveled_->data() + location); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp2_rotation, ntt_table_->data(), modulus_->data(), cfg_ntt,
            d_leveled_->operator[](input1.depth_) * current_rns_mod_count,
            current_rns_mod_count, new_prime_locations + location); // @company CipherFlow 

        // MultSum
        // TODO: make it efficient
        int iteration_count_1 = d_leveled_->operator[](input1.depth_) / 4; // @company CipherFlow 
        int iteration_count_2 = d_leveled_->operator[](input1.depth_) % 4; // @company CipherFlow 
        if (galois_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp2_rotation, galois_key.c_data(),
                temp3_rotation, modulus_->data(), first_rns_mod_count,
                current_decomp_count, current_rns_mod_count, iteration_count_1,
                iteration_count_2, input1.depth_, n_power); // @company CipherFlow 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(galois_key.zero_host_location_,
                                              stream);
            keyswitch_multiply_accumulate_leveled_method_II_kernel<<<
                dim3((n >> 8), current_rns_mod_count, 1), 256, 0, stream>>>(
                temp2_rotation, key_location.data(), temp3_rotation,
                modulus_->data(), first_rns_mod_count, current_decomp_count,
                current_rns_mod_count, iteration_count_1, iteration_count_2,
                input1.depth_, n_power); // @company CipherFlow 
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_NTT_Modulus_Ordered_Inplace(
            temp3_rotation, intt_table_->data(), modulus_->data(), cfg_intt,
            2 * current_rns_mod_count, current_rns_mod_count,
            new_prime_locations + location); // @company CipherFlow 

        // ModDown + Permute
        divide_round_lastq_permute_bfv_kernel<<<
            dim3((n >> 8), current_decomp_count, 2), 256, 0, stream>>>(
            temp3_rotation, temp0_rotation, output_memory.data(),
            modulus_->data(), half_p_->data(), half_mod_->data(),
            last_q_modinv_->data(), galoiselt, n_power, current_rns_mod_count,
            current_decomp_count, first_rns_mod_count, first_decomp_count,
            P_size_); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::switchkey_method_I(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Switchkey<Scheme::BFV>& switch_key, const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        DeviceVector<Data64> temp_rotation((2 * n * Q_size_) +
                                               (n * Q_size_ * Q_prime_size_) +
                                               (2 * n * Q_prime_size_),
                                           stream);
        Data64* temp0_rotation = temp_rotation.data();
        Data64* temp1_rotation = temp0_rotation + (2 * n * Q_size_);
        Data64* temp2_rotation = temp1_rotation + (n * Q_size_ * Q_prime_size_);

        cipher_broadcast_switchkey_kernel<<<dim3((n >> 8), Q_size_, 2), 256, 0,
                                            stream>>>(
            input1.data(), temp0_rotation, temp1_rotation, modulus_->data(),
            n_power, Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Inplace(temp1_rotation, ntt_table_->data(),
                                modulus_->data(), cfg_ntt,
                                Q_size_ * Q_prime_size_, Q_prime_size_);

        // TODO: make it efficient
        int iteration_count_1 = Q_size_ / 4;
        int iteration_count_2 = Q_size_ % 4;
        if (switch_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_rotation, switch_key.data(), temp2_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(switch_key.host_location_,
                                              stream);
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp1_rotation, key_location.data(), temp2_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT_Inplace(temp2_rotation, intt_table_->data(),
                                modulus_->data(), cfg_intt, 2 * Q_prime_size_,
                                Q_prime_size_);

        divide_round_lastq_switchkey_kernel<<<dim3((n >> 8), Q_size_, 2), 256,
                                              0, stream>>>(
            temp2_rotation, temp0_rotation, output_memory.data(),
            modulus_->data(), half_p_->data(), half_mod_->data(),
            last_q_modinv_->data(), n_power, Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::switchkey_method_II(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        Switchkey<Scheme::BFV>& switch_key, const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        DeviceVector<Data64> temp_rotation((2 * n * Q_size_) + (n * Q_size_) +
                                               (2 * n * Q_prime_size_ * d) +
                                               (2 * n * Q_prime_size_),
                                           stream);

        Data64* temp0_rotation = temp_rotation.data();
        Data64* temp1_rotation = temp0_rotation + (2 * n * Q_size_);
        Data64* temp2_rotation = temp1_rotation + (n * Q_size_);
        Data64* temp3_rotation = temp2_rotation + (2 * n * Q_prime_size_ * d);

        cipher_broadcast_switchkey_method_II_kernel<<<
            dim3((n >> 8), Q_size_, 2), 256, 0, stream>>>(
            input1.data(), temp0_rotation, temp1_rotation, modulus_->data(),
            n_power, Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        base_conversion_DtoQtilde_relin_kernel<<<dim3((n >> 8), d, 1), 256, 0,
                                                 stream>>>(
            temp1_rotation, temp2_rotation, modulus_->data(),
            base_change_matrix_D_to_Q_tilda_->data(),
            Mi_inv_D_to_Q_tilda_->data(), prod_D_to_Q_tilda_->data(),
            I_j_->data(), I_location_->data(), n_power, Q_size_, Q_prime_size_,
            d);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT_Inplace(temp2_rotation, ntt_table_->data(),
                                modulus_->data(), cfg_ntt, d * Q_prime_size_,
                                Q_prime_size_);

        // TODO: make it efficient
        int iteration_count_1 = d / 4;
        int iteration_count_2 = d % 4;
        if (switch_key.storage_type_ == storage_type::DEVICE)
        {
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp2_rotation, switch_key.data(), temp3_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }
        else
        {
            DeviceVector<Data64> key_location(switch_key.host_location_,
                                              stream);
            keyswitch_multiply_accumulate_kernel<<<
                dim3((n >> 8), Q_prime_size_, 1), 256, 0, stream>>>(
                temp2_rotation, key_location.data(), temp3_rotation,
                modulus_->data(), n_power, Q_prime_size_, iteration_count_1,
                iteration_count_2);
            HEONGPU_CUDA_CHECK(cudaGetLastError());
        }

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT_Inplace(temp3_rotation, intt_table_->data(),
                                modulus_->data(), cfg_intt, 2 * Q_prime_size_,
                                Q_prime_size_);

        divide_round_lastq_extended_switchkey_kernel<<<
            dim3((n >> 8), Q_size_, 2), 256, 0, stream>>>(
            temp3_rotation, temp0_rotation, output_memory.data(),
            modulus_->data(), half_p_->data(), half_mod_->data(),
            last_q_modinv_->data(), n_power, Q_prime_size_, Q_size_, P_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::negacyclic_shift_poly_coeffmod(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        int index, const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        DeviceVector<Data64> temp(2 * n * Q_size_, stream);

        negacyclic_shift_poly_coeffmod_kernel<<<dim3((n >> 8), Q_size_, 2), 256,
                                                0, stream>>>(
            input1.data(), temp.data(), modulus_->data(), index, n_power);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        // TODO: do with efficient way!
        global_memory_replace_kernel<<<dim3((n >> 8), Q_size_, 2), 256, 0,
                                       stream>>>(temp.data(),
                                                 output_memory.data(), n_power);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::transform_to_ntt_bfv_plain(
        Plaintext<Scheme::BFV>& input1, Plaintext<Scheme::BFV>& output,
        const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory(n * Q_size_, stream);

        DeviceVector<Data64> temp_plain_mul(n * Q_size_, stream);
        Data64* temp1_plain_mul = temp_plain_mul.data();

        threshold_kernel<<<dim3((n >> 8), Q_size_, 1), 256, 0, stream>>>(
            input1.data(), temp1_plain_mul, modulus_->data(),
            upper_halfincrement_->data(), upper_threshold_, n_power, Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT(temp1_plain_mul, output_memory.data(),
                        ntt_table_->data(), modulus_->data(), cfg_ntt, Q_size_,
                        Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::transform_to_ntt_bfv_cipher(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        gpuntt::ntt_rns_configuration<Data64> cfg_ntt = {
            .n_power = n_power,
            .ntt_type = gpuntt::FORWARD,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .stream = stream};

        gpuntt::GPU_NTT(input1.data(), output_memory.data(), ntt_table_->data(),
                        modulus_->data(), cfg_ntt, 2 * Q_size_, Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    __host__ void HEOperator<Scheme::BFV>::transform_from_ntt_bfv_cipher(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        const cudaStream_t stream)
    {
        DeviceVector<Data64> output_memory((2 * n * Q_size_), stream);

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_inverse_->data(),
            .stream = stream};

        gpuntt::GPU_INTT(input1.data(), output_memory.data(),
                        intt_table_->data(), modulus_->data(), cfg_intt,
                        2 * Q_size_, Q_size_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        output.memory_set(std::move(output_memory));
    }

    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    //                       BOOTSRAPPING                         //
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////

    __host__ Ciphertext<Scheme::BFV>
    HEOperator<Scheme::BFV>::operator_from_ciphertext(
        Ciphertext<Scheme::BFV>& input, cudaStream_t stream)
    {
        Ciphertext<Scheme::BFV> cipher;

        cipher.coeff_modulus_count_ = input.coeff_modulus_count_;
        cipher.cipher_size_ = input.cipher_size_;
        cipher.ring_size_ = input.ring_size_;

        cipher.scheme_ = input.scheme_;
        cipher.in_ntt_domain_ = input.in_ntt_domain_;

        cipher.storage_type_ = storage_type::DEVICE;

        cipher.relinearization_required_ = input.relinearization_required_;
        cipher.ciphertext_generated_ = true;

        int cipher_memory_size = 2 * Q_size_ * n;

        cipher.device_locations_ =
            DeviceVector<Data64>(cipher_memory_size, stream);

        return cipher;
    }

    HEArithmeticOperator<Scheme::BFV>::HEArithmeticOperator(
        HEContext<Scheme::BFV>& context, HEEncoder<Scheme::BFV>& encoder)
        : HEOperator<Scheme::BFV>(context, encoder)
    {
    }

    HELogicOperator<Scheme::BFV>::HELogicOperator(
        HEContext<Scheme::BFV>& context, HEEncoder<Scheme::BFV>& encoder)
        : HEOperator<Scheme::BFV>(context, encoder)
    {
        // TODO: make it efficinet
        Data64 constant_1 = 1ULL;
        encoded_constant_one_ = DeviceVector<Data64>(slot_count_);
        fill_device_vector<<<dim3((n >> 8), 1, 1), 256>>>(
            encoded_constant_one_.data(), constant_1, slot_count_);
        HEONGPU_CUDA_CHECK(cudaGetLastError());

        gpuntt::ntt_rns_configuration<Data64> cfg_intt = {
            .n_power = n_power,
            .ntt_type = gpuntt::INVERSE,
            .ntt_layout = gpuntt::PerPolynomial,    
            .reduction_poly = gpuntt::ReductionPolynomial::X_N_plus,
            .zero_padding = false,
            .mod_inverse = n_plain_inverse_->data(),
            .stream = 0};

        gpuntt::GPU_INTT_Inplace(encoded_constant_one_.data(),
                                plain_intt_tables_->data(),
                                plain_modulus_pointer_->data(), cfg_intt, 1, 1);
    }

    __host__ void HELogicOperator<Scheme::BFV>::one_minus_cipher(
        Ciphertext<Scheme::BFV>& input1, Ciphertext<Scheme::BFV>& output,
        const ExecutionOptions& options)
    {
        // TODO: make it efficient
        negate_inplace(input1, options);

        // @company CipherFlow begin ---
        int counter = Q_size_;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        addition_plain_bfv_poly<<<dim3((n >> 8), Q_size_, 2), 256, 0,
                                  options.stream_>>>(
            input1.data(), encoded_constant_one_.data(), output.data(),
            modulus_->data(), plain_modulus_, Q_mod_t_->data() + location, upper_threshold_,
            coeeff_div_plainmod_->data() + location, n_power); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }

    __host__ void HELogicOperator<Scheme::BFV>::one_minus_cipher_inplace(
        Ciphertext<Scheme::BFV>& input1, const ExecutionOptions& options)
    {
        // TODO: make it efficient
        negate_inplace(input1, options);

        // @company CipherFlow begin ---
        int counter = Q_size_;
        int location = 0;
        for (int i = 0; i < input1.depth_; i++)
        {
            location += counter;
            counter--;
        }
        // @company CipherFlow end ---

        addition_plain_bfv_poly_inplace<<<dim3((n >> 8), Q_size_, 1), 256, 0,
                                          options.stream_>>>(
            input1.data(), encoded_constant_one_.data(), input1.data(),
            modulus_->data(), plain_modulus_, Q_mod_t_->data() + location, upper_threshold_,
            coeeff_div_plainmod_->data() + location, n_power); // @company CipherFlow 
        HEONGPU_CUDA_CHECK(cudaGetLastError());
    }

} // namespace heongpu
