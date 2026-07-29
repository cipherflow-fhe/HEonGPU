#pragma once

// Local minimal radix-8 NTT utility subset used by HEonGPU's NTT plugin.

#include <cstdint>
#include <stdexcept>

#include <cuda_runtime_api.h>

namespace heongpu::ntt::phantom_local
{
    inline constexpr dim3 gridDimNTT(4096);
    inline constexpr dim3 blockDimNTT(128);
    inline constexpr std::size_t per_thread_sample_size = 8;
    inline constexpr std::size_t per_block_pad = 4;

    inline std::size_t sample_size(std::size_t n)
    {
        switch (n)
        {
            case 2048:
            case 4096:
                return 64;
            case 8192:
                return 128;
            case 16384:
            case 32768:
            case 65536:
            case 131072:
                return 256;
            default:
                throw std::invalid_argument(
                    "unsupported polynomial degree when selecting sample size");
        }
    }

    template <typename T>
    class DeviceBuffer
    {
    public:
        DeviceBuffer() = default;

        DeviceBuffer(std::size_t n, cudaStream_t stream) : n_(n), stream_(stream)
        {
            if (n_ == 0)
            {
                return;
            }

            const cudaError_t err =
                cudaMallocAsync(&ptr_, n_ * sizeof(T), stream_);
            if (err != cudaSuccess)
            {
                throw std::runtime_error(cudaGetErrorString(err));
            }
        }

        DeviceBuffer(const DeviceBuffer&) = delete;
        DeviceBuffer& operator=(const DeviceBuffer&) = delete;

        DeviceBuffer(DeviceBuffer&& other) noexcept
            : ptr_(other.ptr_), n_(other.n_), stream_(other.stream_)
        {
            other.ptr_ = nullptr;
            other.n_ = 0;
            other.stream_ = nullptr;
        }

        DeviceBuffer& operator=(DeviceBuffer&& other) noexcept
        {
            if (this == &other)
            {
                return *this;
            }

            reset();
            ptr_ = other.ptr_;
            n_ = other.n_;
            stream_ = other.stream_;
            other.ptr_ = nullptr;
            other.n_ = 0;
            other.stream_ = nullptr;
            return *this;
        }

        ~DeviceBuffer()
        {
            reset();
        }

        T* get() const
        {
            return ptr_;
        }

        void reset()
        {
            if (!ptr_)
            {
                return;
            }

            cudaFreeAsync(ptr_, stream_ ? stream_ : cudaStreamPerThread);
            ptr_ = nullptr;
            n_ = 0;
            stream_ = nullptr;
        }

    private:
        T* ptr_ = nullptr;
        std::size_t n_ = 0;
        cudaStream_t stream_ = nullptr;
    };

    class DModulus
    {
    public:
        DModulus() = default;

        DModulus(std::uint64_t value, std::uint64_t ratio0,
                 std::uint64_t ratio1)
            : value_(value), const_ratio_{ratio0, ratio1}
        {
        }

        __device__ __host__ std::uint64_t value() const
        {
            return value_;
        }

    private:
        std::uint64_t value_ = 0;
        std::uint64_t const_ratio_[2] = {0, 0};
    };

    class DNTTTable
    {
    public:
        DNTTTable() = default;
        DNTTTable(const DNTTTable&) = delete;
        DNTTTable(DNTTTable&&) = delete;
        DNTTTable& operator=(const DNTTTable&) = delete;
        DNTTTable& operator=(DNTTTable&&) = delete;
        ~DNTTTable() = default;

        std::uint64_t n() const
        {
            return n_;
        }

        DModulus* modulus() const
        {
            return modulus_.get();
        }

        std::uint64_t* twiddle() const
        {
            return twiddle_.get();
        }

        std::uint64_t* twiddle_shoup() const
        {
            return twiddle_shoup_.get();
        }

        std::uint64_t* itwiddle() const
        {
            return itwiddle_.get();
        }

        std::uint64_t* itwiddle_shoup() const
        {
            return itwiddle_shoup_.get();
        }

        std::uint64_t* n_inv_mod_q() const
        {
            return n_inv_mod_q_.get();
        }

        std::uint64_t* n_inv_mod_q_shoup() const
        {
            return n_inv_mod_q_shoup_.get();
        }

        void init(std::uint64_t n, std::uint64_t size, cudaStream_t stream)
        {
            n_ = n;
            size_ = size;
            modulus_ = DeviceBuffer<DModulus>(size, stream);
            twiddle_ = DeviceBuffer<std::uint64_t>(n * size, stream);
            twiddle_shoup_ = DeviceBuffer<std::uint64_t>(n * size, stream);
            itwiddle_ = DeviceBuffer<std::uint64_t>(n * size, stream);
            itwiddle_shoup_ = DeviceBuffer<std::uint64_t>(n * size, stream);
            n_inv_mod_q_ = DeviceBuffer<std::uint64_t>(size, stream);
            n_inv_mod_q_shoup_ = DeviceBuffer<std::uint64_t>(size, stream);
        }

        void set(const DModulus* modulus, const std::uint64_t* twiddle,
                 const std::uint64_t* twiddle_shoup,
                 const std::uint64_t* itwiddle,
                 const std::uint64_t* itwiddle_shoup,
                 std::uint64_t n_inv_mod_q,
                 std::uint64_t n_inv_mod_q_shoup, std::uint64_t index,
                 cudaStream_t stream) const
        {
            cudaMemcpyAsync(modulus_.get() + index, modulus,
                            sizeof(DModulus), cudaMemcpyHostToDevice, stream);
            cudaMemcpyAsync(twiddle_.get() + index * n_, twiddle,
                            n_ * sizeof(std::uint64_t),
                            cudaMemcpyHostToDevice, stream);
            cudaMemcpyAsync(twiddle_shoup_.get() + index * n_, twiddle_shoup,
                            n_ * sizeof(std::uint64_t),
                            cudaMemcpyHostToDevice, stream);
            cudaMemcpyAsync(itwiddle_.get() + index * n_, itwiddle,
                            n_ * sizeof(std::uint64_t),
                            cudaMemcpyHostToDevice, stream);
            cudaMemcpyAsync(itwiddle_shoup_.get() + index * n_,
                            itwiddle_shoup, n_ * sizeof(std::uint64_t),
                            cudaMemcpyHostToDevice, stream);
            cudaMemcpyAsync(n_inv_mod_q_.get() + index, &n_inv_mod_q,
                            sizeof(std::uint64_t), cudaMemcpyHostToDevice,
                            stream);
            cudaMemcpyAsync(n_inv_mod_q_shoup_.get() + index,
                            &n_inv_mod_q_shoup, sizeof(std::uint64_t),
                            cudaMemcpyHostToDevice, stream);
        }

    private:
        std::uint64_t n_ = 0;
        std::uint64_t size_ = 0;
        DeviceBuffer<DModulus> modulus_;
        DeviceBuffer<std::uint64_t> twiddle_;
        DeviceBuffer<std::uint64_t> twiddle_shoup_;
        DeviceBuffer<std::uint64_t> itwiddle_;
        DeviceBuffer<std::uint64_t> itwiddle_shoup_;
        DeviceBuffer<std::uint64_t> n_inv_mod_q_;
        DeviceBuffer<std::uint64_t> n_inv_mod_q_shoup_;
    };

    __device__ __forceinline__ void csub_q(std::uint64_t& operand,
                                           const std::uint64_t& modulus)
    {
        const std::uint64_t tmp = operand - modulus;
        operand = tmp + (tmp >> 63) * modulus;
    }

    __device__ __forceinline__ std::uint64_t multiply_and_reduce_shoup_lazy(
        const std::uint64_t& operand1, const std::uint64_t& operand2,
        const std::uint64_t& operand2_shoup, const std::uint64_t& modulus)
    {
        const std::uint64_t hi = __umul64hi(operand1, operand2_shoup);
        return operand1 * operand2 - hi * modulus;
    }

    __device__ __forceinline__ void ct_butterfly(
        std::uint64_t& x, std::uint64_t& y, const std::uint64_t& tw,
        const std::uint64_t& tw_shoup, const std::uint64_t& mod)
    {
        const std::uint64_t hi = __umul64hi(y, tw_shoup);
        const std::uint64_t tw_y = y * tw - hi * mod;
        const std::uint64_t mod2 = 2 * mod;
        const std::uint64_t tmp = x - mod2;
        x = tmp + (tmp >> 63) * mod2;
        y = x + mod2 - tw_y;
        x += tw_y;
    }

    __device__ __forceinline__ void gs_butterfly(
        std::uint64_t& x, std::uint64_t& y, const std::uint64_t& tw,
        const std::uint64_t& tw_shoup, const std::uint64_t& mod)
    {
        const std::uint64_t mod2 = 2 * mod;
        const std::uint64_t t = x + mod2 - y;
        std::uint64_t s = x + y;
        csub_q(s, mod2);
        x = s;
        y = multiply_and_reduce_shoup_lazy(t, tw, tw_shoup, mod);
    }

    __device__ __forceinline__ void fntt8(
        std::uint64_t* s, const std::uint64_t* tw,
        const std::uint64_t* tw_shoup, std::uint64_t tw_idx,
        std::uint64_t mod)
    {
        ct_butterfly(s[0], s[4], tw[tw_idx], tw_shoup[tw_idx], mod);
        ct_butterfly(s[1], s[5], tw[tw_idx], tw_shoup[tw_idx], mod);
        ct_butterfly(s[2], s[6], tw[tw_idx], tw_shoup[tw_idx], mod);
        ct_butterfly(s[3], s[7], tw[tw_idx], tw_shoup[tw_idx], mod);
        ct_butterfly(s[0], s[2], tw[2 * tw_idx], tw_shoup[2 * tw_idx],
                     mod);
        ct_butterfly(s[1], s[3], tw[2 * tw_idx], tw_shoup[2 * tw_idx],
                     mod);
        ct_butterfly(s[4], s[6], tw[2 * tw_idx + 1],
                     tw_shoup[2 * tw_idx + 1], mod);
        ct_butterfly(s[5], s[7], tw[2 * tw_idx + 1],
                     tw_shoup[2 * tw_idx + 1], mod);
        ct_butterfly(s[0], s[1], tw[4 * tw_idx], tw_shoup[4 * tw_idx],
                     mod);
        ct_butterfly(s[2], s[3], tw[4 * tw_idx + 1],
                     tw_shoup[4 * tw_idx + 1], mod);
        ct_butterfly(s[4], s[5], tw[4 * tw_idx + 2],
                     tw_shoup[4 * tw_idx + 2], mod);
        ct_butterfly(s[6], s[7], tw[4 * tw_idx + 3],
                     tw_shoup[4 * tw_idx + 3], mod);
    }

    __device__ __forceinline__ void fntt4(
        std::uint64_t* s, const std::uint64_t* tw,
        const std::uint64_t* tw_shoup, std::uint64_t tw_idx,
        std::uint64_t mod)
    {
        ct_butterfly(s[0], s[2], tw[tw_idx], tw_shoup[tw_idx], mod);
        ct_butterfly(s[1], s[3], tw[tw_idx], tw_shoup[tw_idx], mod);
        ct_butterfly(s[0], s[1], tw[2 * tw_idx], tw_shoup[2 * tw_idx],
                     mod);
        ct_butterfly(s[2], s[3], tw[2 * tw_idx + 1],
                     tw_shoup[2 * tw_idx + 1], mod);
    }

    __device__ __forceinline__ void intt8(
        std::uint64_t* s, const std::uint64_t* tw,
        const std::uint64_t* tw_shoup, std::uint64_t tw_idx,
        std::uint64_t mod)
    {
        gs_butterfly(s[0], s[1], tw[4 * tw_idx], tw_shoup[4 * tw_idx],
                     mod);
        gs_butterfly(s[2], s[3], tw[4 * tw_idx + 1],
                     tw_shoup[4 * tw_idx + 1], mod);
        gs_butterfly(s[4], s[5], tw[4 * tw_idx + 2],
                     tw_shoup[4 * tw_idx + 2], mod);
        gs_butterfly(s[6], s[7], tw[4 * tw_idx + 3],
                     tw_shoup[4 * tw_idx + 3], mod);
        gs_butterfly(s[0], s[2], tw[2 * tw_idx], tw_shoup[2 * tw_idx],
                     mod);
        gs_butterfly(s[1], s[3], tw[2 * tw_idx], tw_shoup[2 * tw_idx],
                     mod);
        gs_butterfly(s[4], s[6], tw[2 * tw_idx + 1],
                     tw_shoup[2 * tw_idx + 1], mod);
        gs_butterfly(s[5], s[7], tw[2 * tw_idx + 1],
                     tw_shoup[2 * tw_idx + 1], mod);
        gs_butterfly(s[0], s[4], tw[tw_idx], tw_shoup[tw_idx], mod);
        gs_butterfly(s[1], s[5], tw[tw_idx], tw_shoup[tw_idx], mod);
        gs_butterfly(s[2], s[6], tw[tw_idx], tw_shoup[tw_idx], mod);
        gs_butterfly(s[3], s[7], tw[tw_idx], tw_shoup[tw_idx], mod);
    }

    __device__ __forceinline__ void intt4(
        std::uint64_t* s, const std::uint64_t* tw,
        const std::uint64_t* tw_shoup, std::uint64_t tw_idx,
        std::uint64_t mod)
    {
        gs_butterfly(s[0], s[2], tw[2 * tw_idx], tw_shoup[2 * tw_idx],
                     mod);
        gs_butterfly(s[4], s[6], tw[2 * tw_idx + 1],
                     tw_shoup[2 * tw_idx + 1], mod);
        gs_butterfly(s[0], s[4], tw[tw_idx], tw_shoup[tw_idx], mod);
        gs_butterfly(s[2], s[6], tw[tw_idx], tw_shoup[tw_idx], mod);
    }
} // namespace heongpu::ntt::phantom_local
