// Copyright (c) 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0

//Compare the filter design code with the reference data from spuce ad08df7.
//Tolerances are relative to the largest finite value of the reference array
//(1 for an array of zeros), not to 1: IIR b taps can be 1e-7.
//Tolerances:
// - 1e-9 for the windows and the FIR types with a closed form: the same
//   formulas give the same values, except for the order of operations;
// - 1e-6 for Remez and the Remez estimates: the iterative exchange
//   algorithm can stop at a slightly different point;
// - 1e-6 for maxflat: its sums cancel. With 63 taps and a bandwidth of
//   0.45 (32 poles), the double result has an error of approximately 3e-8
//   (measured against long double), and it depends on the C library;
// - 1e-6 for the IIR designs: the elliptic functions and the pole/zero
//   transforms add more rounding steps;
// - 1e-9 for the IIR filter runs: the taps come from the reference data,
//   and the filter does the same operations as spuce (for the integer
//   types, the result is the same integer).

#include "FilterDesign.hpp"
#include "DirectFormIIR.hpp"
#include "ReferenceData.hpp"

#include <Pothos/Testing.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <format>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace PothosComms::FilterDesign;

namespace {

bool sameValue(const double x, const double ref, const double tol)
{
    if (std::isnan(ref)) return std::isnan(x);
    if (std::isinf(ref)) return x == ref;
    return std::abs(x - ref) <= tol;
}

//the largest finite absolute value, 1 when there is none
double largest(std::span<const double> ref)
{
    double m = 0.0;
    for (const auto r : ref) if (std::isfinite(r)) m = std::max(m, std::abs(r));
    return (m == 0.0) ? 1.0 : m;
}

double largest(std::span<const std::complex<double>> ref)
{
    double m = 0.0;
    for (const auto r : ref) if (std::isfinite(std::abs(r))) m = std::max(m, std::abs(r));
    return (m == 0.0) ? 1.0 : m;
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
        const double tol = (c.type == "remez" or c.type == "maxflat") ? 1e-6 : 1e-9;
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
        //each value relative to itself (they have different scales)
        const std::vector<double> atten{remezEstimateAtten(c.numTaps, c.transBw, c.ripple)};
        const std::vector<double> bw{remezEstimateBw(c.numTaps, c.ripple, c.atten)};
        const std::vector<double> weight{remezEstimateWeight(c.ripple, c.atten)};
        POTHOS_TEST_EQUAL(remezEstimateNumTaps(c.transBw, c.ripple, c.atten), c.estNumTaps);
        POTHOS_TEST_EQUAL(label + "atten " + compare(atten, std::vector<double>{c.estAtten}, 1e-6), label + "atten ");
        POTHOS_TEST_EQUAL(label + "bw " + compare(bw, std::vector<double>{c.estBw}, 1e-6), label + "bw ");
        POTHOS_TEST_EQUAL(label + "weight " + compare(weight, std::vector<double>{c.estWeight}, 1e-6), label + "weight ");
    }
}

POTHOS_TEST_BLOCK("/comms/tests", test_filter_design_iir)
{
    for (const auto &c : Reference::iirCases)
    {
        const auto label = std::format("{} {} order {} fcd {}: ", c.type, c.band, c.order, c.fcd);
        const auto coeffs = designIIR(c.type, c.band, c.order, c.fcd, c.ripple, c.stopAtten, c.centerFreq);
        POTHOS_TEST_EQUAL(label + (coeffs ? "taps" : "nullopt"), label + (c.unknown ? "nullopt" : "taps"));
        if (not coeffs) continue;
        POTHOS_TEST_EQUAL(label + "b " + compare(coeffs->b, c.b, 1e-6), label + "b ");
        POTHOS_TEST_EQUAL(label + "a " + compare(coeffs->a, c.a, 1e-6), label + "a ");
    }
}

namespace {

//the outputs of DirectFormIIR<T> for the input of the reference case, as double
template <typename T>
std::vector<double> runIIR(const Reference::IIRRunCase &c, const T in)
{
    const auto &design = Reference::iirCases[c.iirCase];
    std::vector<double> taps(design.b.begin(), design.b.end());
    taps.insert(taps.end(), design.a.begin(), design.a.end());
    DirectFormIIR<T> filt;
    filt.setTaps(taps);
    std::vector<double> out;
    for (int n = 0; n < 64; n++) out.push_back(double(filt.clock((n == 0 or c.input == "step") ? in : T(0))));
    return out;
}

template <typename T>
std::vector<std::complex<double>> runIIR(const Reference::IIRRunCase &c, const std::complex<T> in)
{
    const auto &design = Reference::iirCases[c.iirCase];
    std::vector<double> taps(design.b.begin(), design.b.end());
    taps.insert(taps.end(), design.a.begin(), design.a.end());
    DirectFormIIR<std::complex<T>> filt;
    filt.setTaps(taps);
    std::vector<std::complex<double>> out;
    for (int n = 0; n < 64; n++)
    {
        const auto y = filt.clock((n == 0 or c.input == "step") ? in : std::complex<T>(0, 0));
        out.emplace_back(double(y.real()), double(y.imag()));
    }
    return out;
}

//compare the real and the complex run of type T with the reference case
template <typename T>
void checkIIRRun(const Reference::IIRRunCase &c)
{
    const auto label = std::format("iirCases[{}] {} {}: ", c.iirCase, c.input, c.type);
    const T x = std::is_integral_v<T> ? T(50) : T(1);
    POTHOS_TEST_EQUAL(label + compare(runIIR<T>(c, x), c.outReal, 1e-9), label);
    POTHOS_TEST_EQUAL(label + compare(runIIR<T>(c, std::complex<T>(x, T(x/2))), c.outComplex, 1e-9), label);
}

} //namespace

POTHOS_TEST_BLOCK("/comms/tests", test_filter_design_iir_run)
{
    for (const auto &c : Reference::iirRunCases)
    {
        if (c.type == "double") checkIIRRun<double>(c);
        else if (c.type == "float") checkIIRRun<float>(c);
        else if (c.type == "int64") checkIIRRun<std::int64_t>(c);
        else if (c.type == "int32") checkIIRRun<std::int32_t>(c);
        else if (c.type == "int16") checkIIRRun<std::int16_t>(c);
        else if (c.type == "int8") checkIIRRun<std::int8_t>(c);
        else POTHOS_TEST_EQUAL(c.type, "a known type");
    }
}

POTHOS_TEST_BLOCK("/comms/tests", test_filter_design_iir_taps)
{
    //POTHOS_TEST_EQUAL can evaluate its arguments more than one time:
    //keep the outputs of clock() in variables
    DirectFormIIR<double> f;
    POTHOS_TEST_THROWS(f.setTaps(std::vector<double>{}), std::invalid_argument);
    POTHOS_TEST_THROWS(f.setTaps(std::vector<double>{1.0}), std::invalid_argument);

    //2 taps: b = {0.5}, no feedback
    f.setTaps(std::vector<double>{0.5, 1.0});
    const double y = f.clock(2.0);
    POTHOS_TEST_EQUAL(y, 1.0);
    POTHOS_TEST_EQUAL(f.feedforward().size(), 1u);
    POTHOS_TEST_EQUAL(f.feedback().size(), 0u);

    //odd count: b = {1, 0}, feedback = {0.5}, the last tap (9) is not used
    f.setTaps(std::vector<double>{1.0, 0.0, 1.0, -0.5, 9.0});
    const double y0 = f.clock(1.0);
    const double y1 = f.clock(0.0);
    const double y2 = f.clock(0.0);
    POTHOS_TEST_EQUAL(y0, 1.0);
    POTHOS_TEST_EQUAL(y1, 0.5);
    POTHOS_TEST_EQUAL(y2, 0.25);
    POTHOS_TEST_EQUAL(f.feedback().size(), 1u);
    POTHOS_TEST_EQUAL(f.feedback()[0], 0.5);

    //reset clears the history
    f.reset();
    const double y3 = f.clock(0.0);
    POTHOS_TEST_EQUAL(y3, 0.0);
}
