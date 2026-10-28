// Copyright (c) 2020 Nicholas Corgan
// SPDX-License-Identifier: BSD-3-Clause

//One copy per SIMD architecture (PothosGenerateSIMDSources)
#include <Pothos/Util/XSIMDCopy.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <type_traits>

// Actually enforce EnableForSIMDConjugate
#ifdef _MSC_VER
#pragma warning(error: 4667) // no function template defined that matches forced instantiation
#endif

namespace PothosCommsSIMD { namespace POTHOS_SIMD_NAMESPACE {

namespace detail
{
    template <typename T>
    struct IsComplex : std::false_type {};

    template <typename T>
    struct IsComplex<std::complex<T>> : std::true_type {};

    // No int16 support due to XSIMD limitation
    template <typename T>
    struct IsSIMDConjugateSupported : std::integral_constant<bool,
        Pothos::Util::XSIMDTraits<T>::IsSupported &&
        !std::is_same<T, std::complex<std::int16_t>>::value &&
        IsComplex<T>::value> {};

    template <typename T>
    using EnableForSIMDConjugate = typename std::enable_if<IsSIMDConjugateSupported<T>::value>::type;

    template <typename T>
    using EnableForDefaultConjugate = typename std::enable_if<!IsSIMDConjugateSupported<T>::value>::type;

    //
    // XSIMD parses complex buffers by grouping all real fields and all
    // imaginary fields, as follows:
    //
    // [real, real, real, real, imag, imag, imag, imag]
    //
    // Given this, there is no suitable XSIMD functionality for what we need
    // here. Instead, we'll parse the input buffers as scalars and apply a
    // mask such that every other element is multiplied by -1.
    //
    // The mask has one lane for each scalar of the SIMD register: false for
    // the real fields, true for the imaginary fields.
    //

    template <typename T>
    static inline xsimd::batch_bool<T> getImagMask()
    {
        std::array<bool, xsimd::batch_bool<T>::size> mask{};
        for (size_t i = 1; i < mask.size(); i += 2) mask[i] = true;
        return xsimd::batch_bool<T>::load_unaligned(mask.data());
    }

    template <typename T>
    static void conjUnoptimized(const T* in, T* out, size_t len)
    {
        for (size_t elem = 0; elem < len; ++elem)
        {
            out[elem] = std::conj(in[elem]);
        }
    }

    template <typename T>
    static EnableForSIMDConjugate<T> conj(const T* in, T* out, size_t len)
    {
        using ScalarType = typename T::value_type;
        static constexpr size_t simdSize = xsimd::simd_traits<ScalarType>::size;

        const auto scalarLen = len * 2;
        const auto numSIMDFrames = scalarLen / simdSize;

        static const auto imagMask = getImagMask<ScalarType>();

        const ScalarType* scalarIn = (const ScalarType*)in;
        const ScalarType* scalarInPtr = (const ScalarType*)in;
        ScalarType* scalarOutPtr = (ScalarType*)out;

        static const auto NegOneReg = xsimd::batch<ScalarType>(-1);

        for (size_t frameIndex = 0; frameIndex < numSIMDFrames; ++frameIndex)
        {
            auto inReg = xsimd::load_unaligned(scalarInPtr);

            auto outReg = xsimd::select(imagMask, (inReg * NegOneReg), inReg);
            outReg.store_unaligned(scalarOutPtr);

            scalarInPtr += simdSize;
            scalarOutPtr += simdSize;
        }

        conjUnoptimized(
            (const T*)scalarInPtr,
            (T*)scalarOutPtr,
            ((scalarLen - (scalarInPtr - scalarIn)) / 2));
    }

    template <typename T>
    static EnableForDefaultConjugate<T> conj(const T* in, T* out, size_t len)
    {
        conjUnoptimized(in, out, len);
    }
}

template <typename T>
void conj(const T* in, T* out, size_t len)
{
    detail::conj(in, out, len);
}

#define CONJ(T) template void conj(const std::complex<T>*, std::complex<T>*, size_t);

    CONJ(std::int8_t)
    CONJ(std::int16_t)
    CONJ(std::int32_t)
    CONJ(std::int64_t)
    CONJ(float)
    CONJ(double)

}}
