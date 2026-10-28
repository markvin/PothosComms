// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

//The filter design code of the comms blocks (internal library comms_filter_design).
//The functions give the same results as the spuce functions of the same name.
//Errors are std::runtime_error with the spuce what() text.

#pragma once

#include <complex>
#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

namespace PothosComms::FilterDesign {

/*!
 * Window taps (spuce design_window).
 * \param type rectangular, hann, hanning, hamming, blackman, bartlett,
 * flattop, kaiser or chebyshev; an unknown type gives rectangular
 * \param arg the beta of kaiser, the sidelobe attenuation (dB) of chebyshev
 */
std::vector<double> designWindow(std::string_view type, int numTaps, double arg = 4.0);

/*!
 * Real FIR taps (spuce design_fir).
 * \param type sinc, raised_cosine, root_raised_cosine, maxflat, gaussian or remez
 * \param band LOW_PASS, HIGH_PASS, BAND_PASS or BAND_STOP
 * \param fl the lower edge, normalized to the sample rate
 * \param fu the upper edge (band pass and band stop), normalized to the sample rate
 * \param alpha the alpha of the raised cosine types, the transition bandwidth of remez
 * \param weight the stop band weight of remez
 */
std::vector<double> designFIR(std::string_view type, std::string_view band,
    int numTaps, double fl, double fu, double alpha = 0.1, double weight = 100.0);

/*!
 * Complex FIR taps (spuce design_complex_fir).
 * \param band COMPLEX_BAND_PASS or COMPLEX_BAND_STOP
 * The other parameters are as in designFIR().
 */
std::vector<std::complex<double>> designComplexFIR(std::string_view type, std::string_view band,
    int numTaps, double fl, double fu, double alpha = 0.1, double weight = 100.0);

//! The number of Remez taps for a transition bandwidth, a pass band ripple (dB) and a stop band attenuation (dB)
std::size_t remezEstimateNumTaps(double transBw, double ripple, double atten);

//! The stop band attenuation (dB) that numTaps Remez taps can give, approximately
double remezEstimateAtten(std::size_t numTaps, double transBw, double ripple);

//! The transition bandwidth that numTaps Remez taps can give, approximately
double remezEstimateBw(std::size_t numTaps, double ripple, double atten);

//! The Remez stop band weight for a pass band ripple (dB) and a stop band attenuation (dB)
double remezEstimateWeight(double passRippleDb, double stopAttenDb);

//! IIR filter coefficients: the feedforward taps b and the feedback taps a (a[0] = 1)
struct IIRCoeffs
{
    std::vector<double> b, a;
};

/*!
 * IIR filter design (spuce design_iir, then iir_coeff::get_b and get_a).
 * \param type butterworth, chebyshev, chebyshev2 or elliptic
 * \param band LOW_PASS, HIGH_PASS, BAND_PASS or BAND_STOP
 * \param order the order of the low pass prototype (band pass and band stop double it)
 * \param fcd the cutoff frequency (low pass, high pass) or the half bandwidth
 * (band pass, band stop), normalized to the sample rate
 * \param ripple the pass band ripple in dB (chebyshev, elliptic)
 * \param stopAtten the stop band attenuation in dB (chebyshev2, elliptic)
 * \param centerFreq the center frequency of band pass and band stop
 * \return nullopt for an unknown type or band (spuce wrote a message and
 * gave zero taps or a low pass)
 */
std::optional<IIRCoeffs> designIIR(std::string_view type, std::string_view band,
    int order, double fcd, double ripple = 0.1, double stopAtten = 60.0,
    double centerFreq = 0.25);

} //namespace PothosComms::FilterDesign
