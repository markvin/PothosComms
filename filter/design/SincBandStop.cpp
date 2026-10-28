// Copyright (c) 2013, Cardinal Peak, LLC. http://www.cardinalpeak.com
//               Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSD-3-Clause
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce), file
// spuce/filters/sinc_helper.cpp, which has this code from the FIR filter
// class of Mike Perkins (http://www.cardinalpeak.com/blog?p=1841) with
// the license below.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1) Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2) Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// 3) Neither the name of Cardinal Peak nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// CARDINAL PEAK, LLC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include "Detail.hpp"

#include <cmath>
#include <numbers>

namespace PothosComms::FilterDesign::Detail {

std::vector<double> sincBandStop(const std::size_t numTaps, const double fl, const double fu)
{
    constexpr double pi = std::numbers::pi;
    std::vector<double> taps(numTaps);
    const auto lambda = pi * fl / (0.5);
    const auto phi = pi * fu / (0.5);
    for (std::size_t n = 0; n < numTaps; n++)
    {
        const double mm = n - (numTaps - 1.0) / 2.0;
        if (mm == 0.0) taps[n] = 1.0 - (phi - lambda) / pi;
        else taps[n] = -(std::sin(mm * phi) - std::sin(mm * lambda)) / (mm * pi);
    }
    return taps;
}

} //namespace PothosComms::FilterDesign::Detail
