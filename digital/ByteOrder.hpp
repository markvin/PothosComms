// Copyright (c) 2020 Nicholas Corgan
// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <bit>
#include <complex>
#include <cstdint>
#include <type_traits>

namespace detail
{
    template <typename T>
    struct IsComplex : std::false_type {};

    template <typename T>
    struct IsComplex<std::complex<T>> : std::true_type {};

    template <typename T>
    struct BufferType
    {
        using Type = T;
    };

    template <>
    struct BufferType<std::int16_t>
    {
        using Type = std::uint16_t;
    };

    template <>
    struct BufferType<std::int32_t>
    {
        using Type = std::uint32_t;
    };

    template <>
    struct BufferType<std::int64_t>
    {
        using Type = std::uint64_t;
    };

    template <>
    struct BufferType<float>
    {
        using Type = std::uint32_t;
    };

    template <>
    struct BufferType<double>
    {
        using Type = std::uint64_t;
    };

    template <typename T>
    struct NoReinterpretCast : std::integral_constant<bool,
        std::is_same<T, typename BufferType<T>::Type>::value &&
        !IsComplex<T>::value> {};

    template <typename T>
    struct MustReinterpretCast : std::integral_constant<bool,
        !std::is_same<T, typename BufferType<T>::Type>::value &&
        !IsComplex<T>::value> {};

    template <typename T>
    static inline T byteswap(T val)
    {
        return std::byteswap(val);
    }
    template <typename T>
    static typename std::enable_if<NoReinterpretCast<T>::value, void>::type byteswapBuffer(const T* in, T* out, size_t numElements)
    {
        for (size_t elem = 0; elem < numElements; ++elem) out[elem] = byteswap(in[elem]);
    }

    template <typename T>
    static typename std::enable_if<MustReinterpretCast<T>::value, void>::type byteswapBuffer(const T* in, T* out, size_t numElements)
    {
        using BufferType = typename BufferType<T>::Type;
        static_assert(sizeof(T) == sizeof(BufferType), "type size mismatch");

        byteswapBuffer<BufferType>((const BufferType*)in, (BufferType*)out, numElements);
    }

    template <typename T>
    static typename std::enable_if<IsComplex<T>::value, void>::type byteswapBuffer(const T* in, T* out, size_t numElements)
    {
        using ScalarType = typename T::value_type;
        byteswapBuffer<ScalarType>((const ScalarType*)in, (ScalarType*)out, (numElements * 2));
    }
}

template <typename T>
static void byteswapBuffer(const T* in, T* out, size_t numElements)
{
    detail::byteswapBuffer(in, out, numElements);
}