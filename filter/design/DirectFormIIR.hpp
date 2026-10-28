// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace PothosComms::FilterDesign {

/*!
 * Direct form IIR filter (spuce iir_df with the parts of spuce fir that it
 * uses): a feedforward FIR part, then a feedback part.
 * T is the sample type (a real or a std::complex type). The sums use
 * double (std::complex<double> for complex T), then each part converts its
 * output to T, as spuce: integer types are truncated.
 */
template <typename T, typename Coeff = double>
class DirectFormIIR
{
    template <typename U> struct IsComplex : std::false_type {};
    template <typename U> struct IsComplex<std::complex<U>> : std::true_type {};
    static constexpr bool isComplex = IsComplex<T>::value;
    using Acc = std::conditional_t<isComplex, std::complex<double>, double>;

public:
    //spuce iir_df(order 1): b = {1, 1}, no feedback (feedback tap 0)
    DirectFormIIR(void):
        _b{Coeff(1), Coeff(1)}, _fb{Coeff(0)}, _zb(2, T(0)), _zfb(1, T(0))
    {
        return;
    }

    /*!
     * Set the taps: the first half is the feedforward taps b, the second
     * half is the feedback taps a. a[0] (= 1) is not used, the other a
     * taps are negated. With an odd count, the last tap is not used.
     * Also clears the history (spuce kept it when the sizes did not
     * change, but IIRFilter always called reset after set_taps).
     * \throws std::invalid_argument for fewer than 2 taps
     */
    void setTaps(std::span<const double> taps)
    {
        if (taps.size() < 2) throw std::invalid_argument("DirectFormIIR::setTaps(): fewer than 2 taps");
        const std::size_t size = taps.size()/2;
        _b.assign(taps.begin(), taps.begin() + size);
        _fb.resize(size - 1);
        for (std::size_t i = 0; i < size - 1; i++) _fb[i] = -taps[i + size + 1];
        _zb.assign(_b.size(), T(0));
        _zfb.assign(_fb.size(), T(0));
    }

    //! Clear the history
    void reset(void)
    {
        std::fill(_zb.begin(), _zb.end(), T(0));
        std::fill(_zfb.begin(), _zfb.end(), T(0));
    }

    //! Filter one sample
    T clock(const T in)
    {
        //feedforward: shift the inputs in, then the sum
        for (std::size_t i = _zb.size() - 1; i > 0; i--) _zb[i] = _zb[i - 1];
        _zb[0] = in;
        Acc sum(0);
        for (std::size_t i = 0; i < _b.size(); i++) sum = sum + _b[i] * toAcc(_zb[i]);
        const T ff = fromAcc(sum);

        //feedback: the sum of the old outputs, then shift the new output in
        sum = Acc(0);
        for (std::size_t i = 0; i < _fb.size(); i++) sum = sum + _fb[i] * toAcc(_zfb[i]);
        const T out = fromAcc(sum + toAcc(ff));
        if (not _zfb.empty())
        {
            for (std::size_t i = _zfb.size() - 1; i > 0; i--) _zfb[i] = _zfb[i - 1];
            _zfb[0] = out;
        }
        return out;
    }

    //! The feedforward taps b
    std::span<const Coeff> feedforward(void) const
    {
        return _b;
    }

    //! The feedback taps: -a[1], -a[2], ...
    std::span<const Coeff> feedback(void) const
    {
        return _fb;
    }

private:
    static Acc toAcc(const T &x)
    {
        if constexpr (isComplex) return Acc(double(x.real()), double(x.imag()));
        else return Acc(x);
    }

    static T fromAcc(const Acc &x)
    {
        if constexpr (isComplex)
        {
            using V = typename T::value_type;
            return T(static_cast<V>(x.real()), static_cast<V>(x.imag()));
        }
        else return static_cast<T>(x);
    }

    std::vector<Coeff> _b, _fb;
    std::vector<T> _zb, _zfb;
};

} //namespace PothosComms::FilterDesign
