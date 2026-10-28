// Copyright (c) 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0

//Writes filter/design/ReferenceData.hpp to stdout.
//The data comes from spuce ad08df7, the filter design library that comms
//used before it got its own filter design code. The build of comms does not
//compile this program. Refer to README.md in this directory.

#include <spuce/filters/design_fir.h>
#include <spuce/filters/design_iir.h>
#include <spuce/filters/design_window.h>
#include <spuce/filters/iir_df.h>
#include <spuce/filters/remez_estimate.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::ostringstream body; //the data arrays and the case arrays
std::size_t nextArray = 0;

std::string number(const double x)
{
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    return std::format("{:.17g}", x);
}

std::string quoted(const std::string &s)
{
    std::string out = "\"";
    for (const char c : s)
    {
        if (c == '\n') out += "\\n";
        else if (c == '"' or c == '\\') out += std::string("\\") + c;
        else out += c;
    }
    return out + "\"";
}

//write a data array, return the expression for its span
std::string realArray(const std::vector<double> &x)
{
    if (x.empty()) return "{}";
    const auto name = std::format("d{}", nextArray++);
    body << "inline constexpr double " << name << "[] = {";
    for (std::size_t i = 0; i < x.size(); i++) body << (i ? ", " : "") << number(x[i]);
    body << "};\n";
    return "Detail::" + name;
}

std::string complexArray(const std::vector<std::complex<double>> &x)
{
    if (x.empty()) return "{}";
    const auto name = std::format("d{}", nextArray++);
    body << "inline constexpr std::complex<double> " << name << "[] = {";
    for (std::size_t i = 0; i < x.size(); i++)
    {
        body << (i ? ", " : "") << "{" << number(x[i].real()) << ", " << number(x[i].imag()) << "}";
    }
    body << "};\n";
    return "Detail::" + name;
}

//spuce writes some messages to std::cout: keep them out of the header
std::string captureCout(const std::function<void()> &fn)
{
    std::ostringstream captured;
    auto *old = std::cout.rdbuf(captured.rdbuf());
    try
    {
        fn();
    }
    catch (...)
    {
        std::cout.rdbuf(old);
        throw;
    }
    std::cout.rdbuf(old);
    return captured.str();
}

const std::vector<std::string> windowTypes{
    "rectangular", "hann", "hamming", "blackman", "bartlett", "flattop", "kaiser", "chebyshev",
    "hanning", //not in the block docs, but design_window knows it
    "foo"};    //unknown: design_window gives rectangular

std::vector<std::string> windowCases;

void addWindowCase(const std::string &type, const int numTaps, const double arg)
{
    const auto taps = spuce::design_window(type, numTaps, arg);
    windowCases.push_back(std::format("    WindowCase{{{}, {}, {}, {}}},\n",
        quoted(type), numTaps, number(arg), realArray(taps)));
}

void makeWindowCases()
{
    for (const auto &type : windowTypes)
    {
        std::vector<double> args{4.0};
        if (type == "kaiser") args = {0.0, 4.0, 8.6};
        if (type == "chebyshev") args = {40.0, 60.0, 100.0};
        for (const int numTaps : {1, 2, 7, 8, 63, 64})
        {
            for (const double arg : args) addWindowCase(type, numTaps, arg);
        }
    }
}

//the same calculation as remez_fir::remez and remez_fir::createDenseGrid
//(GRIDDENSITY 16, positive symmetry). remez allocates gridSize points, but
//createDenseGrid writes the sum of k points, and it writes grid[j - 1]
//after a band with k = 0. gridSize < 1 is safe: remez throws an error.
bool remezGridIsSafe(const int numTaps, const std::vector<double> &bands)
{
    int r = numTaps / 2;
    if (numTaps % 2 != 0) r++;
    int gridSize = 0;
    for (std::size_t i = 0; i < bands.size(); i += 2)
    {
        gridSize += (int)std::floor(0.5 + 2 * r * 16 * (bands[i + 1] - bands[i]));
    }
    if (gridSize < 1) return true;
    const double delf = 0.5 / (16 * r);
    int written = 0;
    for (std::size_t i = 0; i < bands.size(); i += 2)
    {
        const int k = (int)std::floor((bands[i + 1] - bands[i]) / delf + 0.5);
        if (k <= 0) return false;
        written += k;
    }
    return written <= gridSize;
}

std::vector<std::string> firCases;

void addFIRCase(const std::string &type, const std::string &band, const int numTaps,
    const double fl, const double fu, const double alpha, const double weight)
{
    const bool isComplex = band.starts_with("COMPLEX") or band == "complex FOO";
    const std::string spuceBand = band == "complex FOO" ? "FOO" : band;

    //remez writes outside of its buffers for some band edges (undefined
    //behavior, no reference result): do not record these cases
    if (type == "remez")
    {
        double bandwidth = fl;
        if (isComplex) bandwidth = std::abs(0.5*(fu - fl));
        else if (band == "BAND_PASS" or band == "BAND_STOP") bandwidth = 0.5*(fu - fl);
        if (not remezGridIsSafe(numTaps, {0.0, bandwidth, alpha + bandwidth, 0.5}))
        {
            std::cerr << std::format("not recorded: remez {} {} taps, fl {}, fu {}, alpha {}\n",
                band, numTaps, fl, fu, alpha);
            return;
        }
    }
    std::vector<double> taps;
    std::vector<std::complex<double>> complexTaps;
    std::string error;
    try
    {
        if (isComplex) complexTaps = spuce::design_complex_fir(type, spuceBand, numTaps, fl, fu, alpha, weight);
        else taps = spuce::design_fir(type, spuceBand, numTaps, fl, fu, alpha, weight);
    }
    catch (const std::runtime_error &ex)
    {
        error = ex.what();
        if (error.empty()) throw std::logic_error("empty what() of spuce");
    }
    firCases.push_back(std::format("    FIRCase{{{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},\n",
        quoted(type), quoted(spuceBand), isComplex, numTaps, number(fl), number(fu), number(alpha),
        number(weight), realArray(taps), complexArray(complexTaps), quoted(error)));
}

void makeFIRCases()
{
    const std::vector<std::string> types{
        "sinc", "raised_cosine", "root_raised_cosine", "maxflat", "gaussian", "remez"};
    const std::vector<std::string> bands{
        "LOW_PASS", "HIGH_PASS", "BAND_PASS", "BAND_STOP", "COMPLEX_BAND_PASS", "COMPLEX_BAND_STOP"};
    for (const auto &type : types)
    {
        //design_fir reads alpha only for these types, and weight only for remez
        const bool usesAlpha = type == "raised_cosine" or type == "root_raised_cosine" or type == "remez";
        const std::vector<double> alphas = usesAlpha ? std::vector<double>{0.1, 0.35} : std::vector<double>{0.1};
        const std::vector<double> weights = type == "remez" ? std::vector<double>{1.0, 100.0} : std::vector<double>{100.0};
        for (const auto &band : bands)
        {
            std::vector<std::pair<double, double>> edges{{0.1, 0.2}, {0.05, 0.4}};
            if (band.starts_with("COMPLEX")) edges.emplace_back(-0.2, 0.1);
            for (const int numTaps : {15, 16, 63})
                for (const auto &[fl, fu] : edges)
                    for (const double alpha : alphas)
                        for (const double weight : weights)
                            addFIRCase(type, band, numTaps, fl, fu, alpha, weight);
        }
    }

    //parameters that the blocks reject before the design, or that spuce rejects
    addFIRCase("sinc", "LOW_PASS", 0, 0.1, 0.2, 0.1, 100.0);          //0 taps: no taps, no error
    addFIRCase("sinc", "COMPLEX_BAND_PASS", 0, 0.1, 0.2, 0.1, 100.0); //0 taps: no taps, no error
    addFIRCase("sinc", "BAND_PASS", 15, 0.2, 0.1, 0.1, 100.0);        //fl >= fu: spuce uses |fu - fl|
    addFIRCase("raised_cosine", "BAND_STOP", 15, 0.2, 0.2, 0.1, 100.0);
    addFIRCase("foo", "LOW_PASS", 15, 0.1, 0.2, 0.1, 100.0);          //unknown type: error
    addFIRCase("foo", "COMPLEX_BAND_PASS", 15, 0.1, 0.2, 0.1, 100.0);
    addFIRCase("sinc", "FOO", 15, 0.1, 0.2, 0.1, 100.0);              //unknown real band: low pass taps
    addFIRCase("remez", "FOO", 15, 0.1, 0.2, 0.1, 100.0);
    addFIRCase("sinc", "complex FOO", 15, 0.1, 0.2, 0.1, 100.0);      //unknown complex band: error
}

std::vector<std::string> remezEstimateCases;

void makeRemezEstimateCases()
{
    const std::size_t numTaps = 63;
    for (const double transBw : {0.01, 0.05, 0.1})
        for (const double ripple : {0.1, 1.0})
            for (const double atten : {40.0, 80.0})
            {
                remezEstimateCases.push_back(std::format("    RemezEstimateCase{{{}, {}, {}, {}, {}, {}, {}, {}}},\n",
                    number(transBw), number(ripple), number(atten), numTaps,
                    spuce::remez_estimate_num_taps(transBw, ripple, atten),
                    number(spuce::remez_estimate_atten(numTaps, transBw, ripple)),
                    number(spuce::remez_estimate_bw(numTaps, ripple, atten)),
                    number(spuce::remez_estimate_weight(ripple, atten))));
            }
}

struct IIRDesign
{
    std::string type, band;
    int order;
    double fcd;
    std::vector<double> b, a;
};

std::vector<std::string> iirCases;
std::vector<IIRDesign> iirDesigns;

void addIIRCase(const std::string &type, const std::string &band, const int order, const double fcd)
{
    const double ripple = 0.1, stopAtten = 60.0, centerFreq = 0.25;
    std::vector<double> b, a;
    const auto messages = captureCout([&]{
        auto *filt = spuce::design_iir(type, band, order, fcd, ripple, stopAtten, centerFreq);
        b = filt->get_b();
        a = filt->get_a();
        delete filt;
    });
    //design_iir never gives nullptr: for an unknown type or band it writes a
    //message to std::cout and gives zeros (type) or low pass taps (band)
    const bool unknown = not messages.empty();
    iirCases.push_back(std::format("    IIRCase{{{}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},\n",
        quoted(type), quoted(band), order, number(fcd), number(ripple), number(stopAtten),
        number(centerFreq), realArray(b), realArray(a), unknown));
    iirDesigns.push_back(IIRDesign{type, band, order, fcd, b, a});
}

void makeIIRCases()
{
    for (const std::string type : {"butterworth", "chebyshev", "chebyshev2", "elliptic"})
        for (const std::string band : {"LOW_PASS", "HIGH_PASS", "BAND_PASS", "BAND_STOP"})
            for (int order = 1; order <= 8; order++)
                for (const double fcd : {0.05, 0.2})
                    addIIRCase(type, band, order, fcd);
    addIIRCase("foo", "LOW_PASS", 4, 0.2);
    addIIRCase("butterworth", "FOO", 4, 0.2);
}

std::vector<std::string> iirRunCases;

//the taps go into iir_df as IIRDesigner gives them to IIRFilter: b then a
template <typename T>
std::vector<T> runIIR(const IIRDesign &d, const T in, const bool step)
{
    std::vector<double> taps(d.b);
    taps.insert(taps.end(), d.a.begin(), d.a.end());
    spuce::iir_df<T, double> filt;
    filt.set_taps(taps);
    std::vector<T> out;
    for (int n = 0; n < 64; n++) out.push_back(filt.clock((n == 0 or step) ? in : T(0)));
    return out;
}

template <typename T>
std::vector<double> widen(const std::vector<T> &x)
{
    return std::vector<double>(x.begin(), x.end());
}

template <typename T>
std::vector<std::complex<double>> widen(const std::vector<std::complex<T>> &x)
{
    return std::vector<std::complex<double>>(x.begin(), x.end());
}

void makeIIRRunCases()
{
    for (std::size_t i = 0; i < iirDesigns.size(); i++)
    {
        const auto &d = iirDesigns[i];
        if (d.band != "LOW_PASS" or d.fcd != 0.2 or (d.order != 2 and d.order != 5)) continue;
        for (const bool step : {false, true})
        {
            const std::string input = step ? "step" : "impulse";
            iirRunCases.push_back(std::format("    IIRRunCase{{{}, {}, false, {}, {}}},\n",
                i, quoted(input), realArray(runIIR<double>(d, 1.0, step)),
                complexArray(runIIR<std::complex<double>>(d, {1.0, 0.5}, step))));
            iirRunCases.push_back(std::format("    IIRRunCase{{{}, {}, true, {}, {}}},\n",
                i, quoted(input), realArray(widen(runIIR<float>(d, 1.0f, step))),
                complexArray(widen(runIIR<std::complex<float>>(d, {1.0f, 0.5f}, step)))));
        }
    }
}

void writeCases(const std::string &name, const std::vector<std::string> &cases)
{
    body << "\ninline constexpr std::array " << name << "{\n";
    for (const auto &c : cases) body << c;
    body << "};\n";
}

} //namespace

int main()
{
    makeWindowCases();
    makeFIRCases();
    makeRemezEstimateCases();
    makeIIRCases();
    makeIIRRunCases();

    std::cout << R"(// Copyright (c) 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0

//Generated file: do not edit. filter/design/tools/GenerateReferenceData.cpp
//made it from spuce ad08df7 (https://github.com/hayguen/spuce).
//The tests of the filter design code compare their results with this data.

#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <limits>
#include <span>
#include <string_view>

namespace PothosComms::FilterDesign::Reference {

//design_window(type, numTaps, arg)
struct WindowCase
{
    std::string_view type;
    int numTaps;
    double arg;
    std::span<const double> taps;
};

//design_fir (isComplex false) or design_complex_fir (isComplex true)
struct FIRCase
{
    std::string_view type, band;
    bool isComplex;
    int numTaps;
    double fl, fu, alpha, weight;
    std::span<const double> taps;
    std::span<const std::complex<double>> complexTaps;
    std::string_view error; //not empty: what() of the std::runtime_error
};

//remez_estimate_num_taps(transBw, ripple, atten), remez_estimate_atten(numTaps, transBw, ripple),
//remez_estimate_bw(numTaps, ripple, atten), remez_estimate_weight(ripple, atten)
struct RemezEstimateCase
{
    double transBw, ripple, atten;
    std::size_t numTaps;
    std::size_t estNumTaps;
    double estAtten, estBw, estWeight;
};

//design_iir(type, band, order, fcd, ripple, stopAtten, centerFreq), then get_b() and get_a()
struct IIRCase
{
    std::string_view type, band;
    int order;
    double fcd, ripple, stopAtten, centerFreq;
    std::span<const double> b, a;
    bool unknown; //spuce wrote "unknown type" or "unsupported band" to std::cout
};

//iir_df<T, double> with set_taps(b then a) of iirCases[iirCase], 64 output samples.
//Input: 1 (real) and (1, 0.5) (complex) at n = 0 (impulse) or at all n (step).
//singlePrecision: T is float and std::complex<float>, else double and std::complex<double>.
struct IIRRunCase
{
    std::size_t iirCase;
    std::string_view input;
    bool singlePrecision;
    std::span<const double> outReal;
    std::span<const std::complex<double>> outComplex;
};

namespace Detail {
inline constexpr double nan = std::numeric_limits<double>::quiet_NaN();
inline constexpr double inf = std::numeric_limits<double>::infinity();
)";
    //the data arrays first, then the case arrays (outside of Detail)
    const auto writeAll = [] {
        std::ostringstream data;
        data.swap(body);
        writeCases("windowCases", windowCases);
        writeCases("firCases", firCases);
        writeCases("remezEstimateCases", remezEstimateCases);
        writeCases("iirCases", iirCases);
        writeCases("iirRunCases", iirRunCases);
        std::cout << data.str() << "} //namespace Detail\n" << body.str();
    };
    writeAll();
    std::cout << "\n} //namespace PothosComms::FilterDesign::Reference\n";
    return 0;
}
