// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

#include "FilterDesign.hpp"
#include "Detail.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

namespace PothosComms::FilterDesign {

namespace {

constexpr double pi = std::numbers::pi;

//windowless sinc low pass, normalized to a sum of 1 (spuce sinc_fir)
std::vector<double> sincFIR(const int numTaps, const double f)
{
    const std::size_t len = numTaps;
    std::vector<double> taps(len);
    double sum = 0;
    for (std::size_t i = 0; i < len; i++)
    {
        const double x = i - (len - 1.0) / 2.0;
        if (x == 0.0) taps[i] = f;
        else taps[i] = std::sin(x * f * 2.0*pi) / (x * 2.0*pi);
        sum += taps[i];
    }
    for (auto &tap : taps) tap /= sum;
    return taps;
}

//Gaussian low pass, bt = bandwidth time product (spuce gaussian_fir)
std::vector<double> gaussianFIR(const int numTaps, const double bt)
{
    const double s = 1.0 / (std::sqrt(std::log(2.0)) / (2.0 * pi * bt));
    double t0 = -0.5 * numTaps;
    double gsum = 0.0;
    std::vector<double> taps(numTaps);
    for (int i = 0; i < numTaps; i++)
    {
        t0++;
        const double ts = s * t0;
        taps[i] = std::exp(-0.5 * ts * ts);
        gsum += taps[i];
    }
    gsum = 1.0 / gsum;
    for (auto &tap : taps) tap *= gsum;
    return taps;
}

//one tap of the raised cosine impulse response (spuce raised_cosine_imp)
double raisedCosineTap(const double alpha, const double xin, const double r, const long numTaps)
{
    const double rate = 0.5*r;
    const double xindx = xin - numTaps / 2;
    const double x1 = pi * xindx / rate;
    double x2 = 1 - (4 * alpha * alpha * (xindx / rate) * (xindx / rate));
    if (x1 == 0) return 1;
    if (x2 == 0)
    {
        x2 = 8 * alpha * (xindx / rate) * (xindx / rate);
        return std::sin(x1) * std::sin(alpha * x1) / x2;
    }
    return (std::sin(x1) * std::cos(alpha * x1)) / (x1 * x2);
}

//raised cosine, normalized to a sum of 1 (spuce raised_cosine)
std::vector<double> raisedCosineFIR(const int numTaps, const double alpha, const double rate)
{
    double sum = 0;
    for (int i = 0; i < numTaps; i++) sum += raisedCosineTap(alpha, double(i), rate, numTaps);
    std::vector<double> taps(numTaps);
    for (int i = 0; i < numTaps; i++) taps[i] = (1.0/sum)*raisedCosineTap(alpha, double(i), rate, numTaps);
    return taps;
}

//one tap of the root raised cosine impulse response (spuce root_raised_cosine_imp)
double rootRaisedCosineTap(const double alpha, const double xin, const double r, const long numTaps)
{
    const double rate = 0.5*r;
    const double xindx = xin - numTaps / 2;
    double x1 = pi * xindx / rate;
    double x2 = 4 * alpha * xindx / rate;
    double x3 = x2 * x2 - 1;
    double nom, denom;
    if (x3 != 0)
    {
        if (x1 != 0)
        {
            nom = std::cos((1 + alpha) * x1);
            nom += std::sin((1 - alpha) * x1) / (4 * alpha * xindx / rate);
        }
        else
        {
            nom = std::cos((1 + alpha) * x1);
            nom += (1 - alpha) * pi / (4 * alpha);
        }
        denom = x3 * pi;
    }
    else
    {
        if (alpha == 1) return -1;
        x3 = (1 - alpha) * x1;
        x2 = (1 + alpha) * x1;
        nom = std::sin(x2) * (1 + alpha) * pi - std::cos(x3) * ((1 - alpha) * pi * rate) / (4 * alpha * xindx) +
            std::sin(x3) * rate * rate / (4 * alpha * xindx * xindx);
        denom = -32 * pi * alpha * alpha * xindx / rate;
    }
    return 4 * alpha * nom / denom;
}

//root raised cosine, normalized to a sum of 1 (spuce root_raised_cosine)
std::vector<double> rootRaisedCosineFIR(const int numTaps, const double alpha, const double rate)
{
    double gain = 1.0 / rate;
    double sum = 0;
    for (int i = 0; i < numTaps; i++) sum += gain * rootRaisedCosineTap(alpha, double(i), rate, numTaps);
    gain /= sum;
    std::vector<double> taps(numTaps);
    for (int i = 0; i < numTaps; i++) taps[i] = gain * rootRaisedCosineTap(alpha, double(i), rate, numTaps);
    return taps;
}

/*!
 * Maximally flat (Butterworth) low pass from the impulse response of the
 * analog Butterworth filter (spuce butterworth_fir). spb = the bandwidth.
 * B(p) = sum over k = 1..N/2 of lambda(k)*(1 + alpha(k)*p/wc)/(1 + 2*cos(beta(k))*p/wc + (p/wc)^2),
 * N = the even order, beta(k) = (2*k - 1)/(2*N),
 * alpha(k) = sin((N/2 - 1)*beta(k))/sin(N*beta(k)/2),
 * lambda(k) = sin(N*beta(k)/2)/sin(beta(k)) / product over m != k of 2*(cos(beta(k)) - cos(beta(m))).
 */
std::vector<double> maxflatFIR(const int numTaps, const double spb)
{
    const int spbi = (int)std::floor(1.0 / spb + 0.5);
    if (spbi == 0) throw std::runtime_error("maxflat FIR bandwidth is too large");
    int ord = (int)std::floor(numTaps / spbi + 0.5);
    if (ord % 2) ord += 1; //make even
    const int end = ord / 2;
    const double xend = end;
    std::vector<double> beta(end), alpha(end), ak(end), xk(end), yk(end), wk(end), lamda(end);
    for (int i = 0; i < end; i++)
    {
        const double x = i + 1.0;
        beta[i] = pi * (2.0 * x - 1.0) / (2.0 * ord);
        lamda[i] = std::sin(xend * beta[i]) / std::sin(beta[i]);
        alpha[i] = std::sin((xend - 1.0) * beta[i]) / std::sin(xend * beta[i]);
        ak[i] = std::cos(beta[i]);
        wk[i] = std::sqrt(1 - ak[i] * ak[i]);
        xk[i] = alpha[i];
        yk[i] = (1.0 - alpha[i] * ak[i]) / wk[i];
    }
    for (int i = 0; i < end; i++)
    {
        for (int i2 = 0; i2 < end; i2++)
        {
            if (i2 != i) lamda[i] /= (2.0 * (ak[i] - std::cos(beta[i2])));
        }
    }
    std::vector<double> taps(numTaps);
    double sum = 0;
    for (int j = 0; j < numTaps; j++)
    {
        double ht = 0.0;
        const double t = 2*pi * j / (double)spbi;
        for (int i = 0; i < end; i++)
        {
            ht += lamda[i] * std::exp(-t * ak[i]) * (xk[i] * std::cos(t * wk[i]) + yk[i] * std::sin(t * wk[i]));
        }
        taps[j] = 2.0*pi*ht;
        sum += taps[j];
    }
    for (auto &tap : taps) tap /= sum;
    return taps;
}

//move low pass taps to another band (spuce transform_fir)
std::vector<double> transformFIR(const std::string_view band, const std::vector<double> &taps, const double f)
{
    std::vector<double> out(taps.size());
    if (band == "LOW_PASS")
    {
        out = taps;
    }
    else if (band == "HIGH_PASS")
    {
        //alternate +1, -1
        double h = 1;
        for (std::size_t i = 0; i < taps.size(); i++)
        {
            out[i] = h*taps[i];
            h *= -1;
        }
    }
    else if (band == "BAND_PASS")
    {
        //modulate the taps up to the band center frequency, scale by 2 for the cosine
        for (std::size_t i = 0; i < taps.size(); i++)
        {
            out[i] = 2.0 * std::cos(2.0*pi*f*(i - (taps.size()-1)/2.0)) * taps[i];
        }
    }
    else if (band == "BAND_STOP")
    {
        //band pass, then negate and add 1 at the center tap
        out = transformFIR("BAND_PASS", taps, f);
        for (auto &tap : out) tap *= -1;
        if (not out.empty()) out[(taps.size()-1)/2] += 1;
    }
    else if (band == "INVERSE")
    {
        for (std::size_t i = 0; i < taps.size(); i++) out[i] = -taps[i];
        if (not out.empty()) out[(taps.size()-1)/2] += 1;
    }
    else
    {
        throw std::runtime_error("Unsupported band type in transform_fir :" + std::string(band));
    }
    return out;
}

//move real low pass taps to a complex band (spuce transform_complex_fir)
std::vector<std::complex<double>> transformComplexFIR(const std::string_view band, const std::vector<double> &taps, const double f)
{
    std::vector<std::complex<double>> out(taps.size());
    if (band == "COMPLEX_BAND_PASS")
    {
        //modulate the taps up to the band center frequency
        for (std::size_t i = 0; i < taps.size(); i++)
        {
            const double phase = 2.0*pi*f*(i - (taps.size()-1)/2.0);
            out[i] = taps[i]*std::polar(1.0, phase);
        }
    }
    else if (band == "COMPLEX_BAND_STOP")
    {
        out = transformComplexFIR("COMPLEX_BAND_PASS", taps, f);
        for (auto &tap : out) tap *= -1;
        if (not out.empty()) out[(taps.size()-1)/2] += 1;
    }
    else
    {
        throw std::runtime_error("Unsupported band type in transform_complex_fir :" + std::string(band));
    }
    return out;
}

std::vector<double> remezFIR(const std::string_view band, const int numTaps,
    const double fl, const double fu, const double alpha, const double weight)
{
    //for remez, the bandwidth is fl, except for band pass and band stop
    const bool isBand = band == "BAND_PASS" or band == "BAND_STOP";
    const double bw = isBand ? 0.5*(fu - fl) : fl;
    std::vector<double> des(2), w(2);
    if (band == "HIGH_PASS")
    {
        des[0] = 0.0;
        des[1] = 1.0;
        w[0] = weight;
        w[1] = 1.0;
    }
    else
    {
        des[0] = 1.0;
        des[1] = 0.0;
        w[0] = (band == "BAND_STOP") ? weight : 1.0;
        w[1] = (band == "BAND_STOP") ? 1.0 : weight;
    }
    auto taps = Detail::remez(numTaps, {0, bw, alpha + bw, 0.5}, des, w);
    if (isBand) taps = transformFIR(band, taps, 0.5*(fu + fl));
    return taps;
}

} //namespace

std::vector<double> designFIR(const std::string_view type, const std::string_view band,
    const int numTaps, const double fl, const double fu, const double alpha, const double weight)
{
    if (type == "remez") return remezFIR(band, numTaps, fl, fu, alpha, weight);

    const double centerFrequency = 0.5*(fu + fl);
    double bw = fl;
    double irate = 1.0/fl;
    if (band == "BAND_PASS" or band == "BAND_STOP")
    {
        //half of the band pass bandwidth (double sided), as an absolute value
        bw = std::abs(0.5*(fu - fl));
        irate = 1.0/std::abs(0.5*(fu - fl));
    }
    else if (band == "HIGH_PASS")
    {
        bw = 0.5 - fl;
    }

    std::vector<double> taps;
    if (type == "maxflat") taps = maxflatFIR(numTaps, std::max(bw, 0.02));
    else if (type == "gaussian") taps = gaussianFIR(numTaps, (band == "HIGH_PASS") ? fl : bw);
    else if (type == "raised_cosine") taps = raisedCosineFIR(numTaps, alpha, irate);
    else if (type == "root_raised_cosine") taps = rootRaisedCosineFIR(numTaps, alpha, irate);
    else if (type == "sinc") taps = sincFIR(numTaps, bw);
    else throw std::runtime_error("Unknown fir type : " + std::string(type));

    if (taps.empty()) return taps;
    if (band == "HIGH_PASS")
    {
        //these types do not work with the high pass transform: use the inverse
        const bool inverse = type == "raised_cosine" or type == "root_raised_cosine" or type == "gaussian";
        return transformFIR(inverse ? "INVERSE" : "HIGH_PASS", taps, centerFrequency);
    }
    if (band == "BAND_PASS") return transformFIR("BAND_PASS", taps, centerFrequency);
    if (band == "BAND_STOP")
    {
        if (type == "sinc") return Detail::sincBandStop(numTaps, fl, fu);
        if (type == "maxflat") throw std::runtime_error("maxflat FIR as BAND_STOP not supported\n");
        return transformFIR("BAND_STOP", taps, centerFrequency);
    }
    return taps;
}

std::vector<std::complex<double>> designComplexFIR(const std::string_view type, const std::string_view band,
    const int numTaps, const double fl, const double fu, const double alpha, const double weight)
{
    const double centerFrequency = 0.5*(fu + fl);
    //half of the band pass bandwidth (double sided), as an absolute value
    const double bw = std::abs(0.5*(fu - fl));

    if (type == "sinc" and band == "COMPLEX_BAND_STOP")
    {
        const auto taps = designFIR(type, "HIGH_PASS", numTaps, bw, fu, alpha, weight);
        return transformComplexFIR("COMPLEX_BAND_PASS", taps, centerFrequency);
    }
    const auto taps = designFIR(type, "LOW_PASS", numTaps, bw, fu, alpha, 1.0/weight);
    auto complexTaps = transformComplexFIR(band, taps, centerFrequency);
    if (type == "maxflat" and band == "COMPLEX_BAND_STOP")
    {
        throw std::runtime_error("maxflat FIR as COMPLEX_BAND_STOP not supported\n");
    }
    return complexTaps;
}

} //namespace PothosComms::FilterDesign
