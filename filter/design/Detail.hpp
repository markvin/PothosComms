// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

//Internal functions of the filter design code.

#pragma once

#include <cstddef>
#include <vector>

namespace PothosComms::FilterDesign::Detail {

//windowed sinc band stop taps, fl and fu normalized to the sample rate (spuce sincBSF)
std::vector<double> sincBandStop(std::size_t numTaps, double fl, double fu);

/*!
 * Parks-McClellan (Remez exchange) design of a symmetric FIR filter with
 * 2 bands (spuce remez_fir::remez with remez_type::BANDPASS).
 * \param bands the band edges {0, pass edge, stop edge, 0.5}, 2 per band
 * \param des the desired response of each band
 * \param weight the error weight of each band
 * \throws std::runtime_error when the dense grid is too small or not
 * correct for the band edges, or when the algorithm does not converge
 */
std::vector<double> remez(int numTaps, const std::vector<double> &bands,
    const std::vector<double> &des, const std::vector<double> &weight);

} //namespace PothosComms::FilterDesign::Detail
