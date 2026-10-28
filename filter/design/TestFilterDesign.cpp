// Copyright (c) 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0

//Compare the filter design code with the reference data from spuce ad08df7.
//Tolerances are relative to max(1, the largest reference value):
// - 1e-9 for the windows and the FIR types with a closed form: the same
//   formulas give the same values, except for the order of operations;
// - 1e-6 for Remez and the Remez estimates: the iterative exchange
//   algorithm can stop at a slightly different point.

#include "FilterDesign.hpp"
#include "ReferenceData.hpp"

#include <Pothos/Testing.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <format>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

using namespace PothosComms::FilterDesign;

namespace {

bool sameValue(const double x, const double ref, const double tol)
{
    if (std::isnan(ref)) return std::isnan(x);
    if (std::isinf(ref)) return x == ref;
    return std::abs(x - ref) <= tol;
}

double largest(std::span<const double> ref)
{
    double m = 1.0;
    for (const auto r : ref) if (std::isfinite(r)) m = std::max(m, std::abs(r));
    return m;
}

double largest(std::span<const std::complex<double>> ref)
{
    double m = 1.0;
    for (const auto r : ref) if (std::isfinite(std::abs(r))) m = std::max(m, std::abs(r));
    return m;
}

//an empty text if out agrees with ref, else a description of the first difference
std::string compare(std::span<const double> out, std::span<const double> ref, const double relTol)
{
    if (out.size() != ref.size()) return std::format("size {} != {}", out.size(), ref.size());
    const auto tol = relTol*largest(ref);
    for (std::size_t i = 0; i < ref.size(); i++)
    {
        if (not sameValue(out[i], ref[i], tol)) return std::format("[{}] {} != {}", i, out[i], ref[i]);
    }
    return "";
}

std::string compare(std::span<const std::complex<double>> out, std::span<const std::complex<double>> ref, const double relTol)
{
    if (out.size() != ref.size()) return std::format("size {} != {}", out.size(), ref.size());
    const auto tol = relTol*largest(ref);
    for (std::size_t i = 0; i < ref.size(); i++)
    {
        if (not sameValue(out[i].real(), ref[i].real(), tol) or not sameValue(out[i].imag(), ref[i].imag(), tol))
        {
            return std::format("[{}] ({}, {}) != ({}, {})", i, out[i].real(), out[i].imag(), ref[i].real(), ref[i].imag());
        }
    }
    return "";
}

//what() of the std::runtime_error of fn, an empty text without an error
std::string errorOf(const std::function<void()> &fn)
{
    try
    {
        fn();
    }
    catch (const std::runtime_error &ex)
    {
        return ex.what();
    }
    return "";
}

} //namespace

POTHOS_TEST_BLOCK("/comms/tests", test_filter_design_window)
{
    for (const auto &c : Reference::windowCases)
    {
        const auto label = std::format("{} {} taps arg {}: ", c.type, c.numTaps, c.arg);
        const auto taps = designWindow(c.type, c.numTaps, c.arg);
        POTHOS_TEST_EQUAL(label + compare(taps, c.taps, 1e-9), label);
    }
}

POTHOS_TEST_BLOCK("/comms/tests", test_filter_design_fir)
{
    for (const auto &c : Reference::firCases)
    {
        const auto label = std::format("{} {} {} taps fl {} fu {} alpha {} weight {}: ",
            c.type, c.band, c.numTaps, c.fl, c.fu, c.alpha, c.weight);
        const double tol = (c.type == "remez") ? 1e-6 : 1e-9;
        std::vector<double> taps;
        std::vector<std::complex<double>> complexTaps;
        const auto error = errorOf([&]{
            if (c.isComplex) complexTaps = designComplexFIR(c.type, c.band, c.numTaps, c.fl, c.fu, c.alpha, c.weight);
            else taps = designFIR(c.type, c.band, c.numTaps, c.fl, c.fu, c.alpha, c.weight);
        });
        POTHOS_TEST_EQUAL(label + error, label + std::string(c.error));
        if (not c.error.empty()) continue;
        const auto diff = c.isComplex ? compare(complexTaps, c.complexTaps, tol) : compare(taps, c.taps, tol);
        POTHOS_TEST_EQUAL(label + diff, label);
    }

    //spuce writes outside of its buffers for these Remez band edges (the
    //reference data has no result): the own code gives an error
    POTHOS_TEST_THROWS(designFIR("remez", "BAND_PASS", 15, 0.05, 0.4, 0.35, 1.0), std::runtime_error);
    POTHOS_TEST_THROWS(designFIR("remez", "LOW_PASS", 15, 0.2, 0.2, 0.3, 1.0), std::runtime_error);
    POTHOS_TEST_THROWS(designComplexFIR("remez", "COMPLEX_BAND_PASS", 15, -0.2, 0.1, 0.35, 1.0), std::runtime_error);

    //spuce divides by 0 (int) for a maxflat bandwidth above 2
    POTHOS_TEST_THROWS(designFIR("maxflat", "LOW_PASS", 15, 3.0, 3.5), std::runtime_error);
}

POTHOS_TEST_BLOCK("/comms/tests", test_filter_design_remez_estimate)
{
    for (const auto &c : Reference::remezEstimateCases)
    {
        const auto label = std::format("transBw {} ripple {} atten {}: ", c.transBw, c.ripple, c.atten);
        const std::vector<double> out{
            remezEstimateAtten(c.numTaps, c.transBw, c.ripple),
            remezEstimateBw(c.numTaps, c.ripple, c.atten),
            remezEstimateWeight(c.ripple, c.atten)};
        const std::vector<double> ref{c.estAtten, c.estBw, c.estWeight};
        POTHOS_TEST_EQUAL(remezEstimateNumTaps(c.transBw, c.ripple, c.atten), c.estNumTaps);
        POTHOS_TEST_EQUAL(label + compare(out, ref, 1e-6), label);
    }
}
