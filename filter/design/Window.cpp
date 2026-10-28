// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

#include "FilterDesign.hpp"

#include <cmath>
#include <complex>
#include <numbers>

namespace PothosComms::FilterDesign {

namespace {

constexpr double pi = std::numbers::pi;

//y[l] = sum of x[j]*cos(2*pi*l*j/n)
void realDFT(std::vector<double> &y, const int n)
{
    const std::vector<double> x(y.begin(), y.begin() + n);
    for (int l = 0; l < n; l++)
    {
        y[l] = 0;
        for (int j = 0; j < n; j++) y[l] += x[j] * std::cos(2*pi * l * j / n);
    }
}

//y[l] = sum of x[j]*exp(i*2*pi*l*j/n)
void realDFT(std::vector<std::complex<double>> &y, const int n)
{
    const std::vector<std::complex<double>> x(y.begin(), y.begin() + n);
    for (int l = 0; l < n; l++)
    {
        y[l] = 0;
        for (int j = 0; j < n; j++)
        {
            y[l] += x[j] * std::complex<double>(std::cos(2*pi*l*j/n), std::sin(2*pi*l*j/n));
        }
    }
}

//modified Bessel function of the first kind, order 0 (for the kaiser window)
double besselI0(const double x)
{
    const double t = 1.e-08;
    const double y = 0.5 * x;
    double e = 1.0;
    double de = 1.0;
    for (int i = 1; i < 26; i++)
    {
        de *= y / double(i);
        const double sde = de * de;
        e += sde;
        if ((e * t - sde) > 0) break;
    }
    return e;
}

std::vector<double> rectangular(const int nf)
{
    return std::vector<double>(nf, 1.0);
}

std::vector<double> flattop(const int nf)
{
    const double a0 = 1.0, a1 = 1.93, a2 = 1.29, a3 = 0.388, a4 = 0.028;
    std::vector<double> w(nf);
    for (int i = 0; i < nf; i++)
    {
        w[i] = a0 - a1*std::cos((2.0*pi*i)/(nf-1))
            + a2*std::cos((4.0*pi*i)/(nf-1))
            - a3*std::cos((6.0*pi*i)/(nf-1))
            + a4*std::cos((8.0*pi*i)/(nf-1));
    }
    return w;
}

//w(n) = 0.54 - 0.46*cos(2*pi*n/(nf-1))
std::vector<double> hamming(const int nf)
{
    std::vector<double> w(nf);
    for (int i = 0; i < nf; i++) w[i] = 0.54 - 0.46*std::cos(2*pi * i / (nf-1));
    return w;
}

//w(n) = 0.5*(1 - cos(2*pi*n/(nf-1)))
std::vector<double> hann(const int nf)
{
    std::vector<double> w(nf);
    for (int i = 0; i < nf; i++) w[i] = 0.5 * (1 - std::cos(2*pi * i / (nf-1)));
    return w;
}

//w(n) = 0.5*(1 - cos(2*pi*(n+1)/(nf+1)))
std::vector<double> hanning(const int nf)
{
    std::vector<double> w(nf);
    for (int i = 0; i < nf; i++)
    {
        const double xi = i+1;
        w[i] = 0.5 * (1 - std::cos(2*pi * xi / (nf+1)));
    }
    return w;
}

//w(n) = 0.42 - 0.5*cos(2*pi*n/(nf-1)) + 0.08*cos(4*pi*n/(nf-1))
std::vector<double> blackman(const int nf)
{
    std::vector<double> w(nf);
    for (int i = 0; i < nf; i++)
    {
        w[i] = 0.42 - 0.5 * std::cos(2*pi * i / (nf-1)) + 0.08 * std::cos(2 * 2*pi * i / (nf-1));
    }
    return w;
}

std::vector<double> kaiser(const int nf, const double beta)
{
    std::vector<double> w(nf);
    const double bes = 1.0 / besselI0(beta);
    for (int i = 0; i < nf; i++)
    {
        const double xi = i - nf/2.0 + 0.5;
        const double bm = (2.0*xi/(nf - 1));
        w[i] = besselI0(beta * std::sqrt(1. - bm*bm)) * bes;
    }
    return w;
}

//Chebyshev polynomial of order m
double chebyPoly(const int m, const double a)
{
    const int odd = (m%2) ? -1 : 1;
    if (std::fabs(a) > 1.0) return odd*std::cosh(m * std::acosh(a));
    return odd*std::cos(m * std::acos(a));
}

//Dolph-Chebyshev window, rDb = the sidelobe attenuation
std::vector<double> chebyshev(const int nf, const double rDb)
{
    const double r = std::pow(10.0, rDb/20.0);
    const double x0 = std::cosh(std::acosh(r)/(nf-1.0));
    std::vector<double> w(nf);
    for (int i = 0; i < nf; i++)
    {
        const double a = std::fabs(x0*std::cos(i*pi/nf));
        w[i] = chebyPoly(nf-1, a);
        if (nf%2 == 0 and i > nf/2) w[i] *= -1;
    }
    if (nf%2 == 0)
    {
        //for even, the DFT needs complex input
        std::vector<std::complex<double>> cw(nf);
        for (int j = 0; j < nf; j++) cw[j] = w[j] * std::complex<double>(std::cos(pi*j/nf), std::sin(pi*j/nf));
        realDFT(cw, nf);
        for (int i = 0; i < nf; i++) w[i] = std::real(cw[i]);
        for (int i = 0; i < nf/2; i++) w[i+nf/2] = w[i];
        for (int i = 0; i < nf/2; i++) w[nf/2-1-i] = w[i+nf/2];
    }
    else
    {
        realDFT(w, nf);
        for (int i = 0; i < nf/2+1; i++) w[nf-1-i] = w[nf/2-i];
        for (int i = 0; i < nf/2; i++) w[nf/2-1-i] = w[nf/2+1+i];
    }
    double maxV = 0;
    for (int i = 0; i < nf; i++) if (std::fabs(w[i]) > maxV) maxV = std::fabs(w[i]);
    for (int i = 0; i < nf; i++) w[i] /= maxV;
    return w;
}

std::vector<double> bartlett(const int nf)
{
    std::vector<double> w(nf);
    for (int i = 0; i < nf / 2; i++)
    {
        const double win = 2.0 * i / (nf - 1);
        w[i] = win;
        w[nf - 1 - i] = win;
    }
    if (nf % 2 == 1) w[nf / 2] = 1.0;
    return w;
}

} //namespace

std::vector<double> designWindow(const std::string_view type, const int numTaps, const double arg)
{
    if (type == "hamming") return hamming(numTaps);
    if (type == "hanning") return hanning(numTaps);
    if (type == "hann") return hann(numTaps);
    if (type == "blackman") return blackman(numTaps);
    if (type == "kaiser") return kaiser(numTaps, arg);
    if (type == "bartlett") return bartlett(numTaps);
    if (type == "chebyshev") return chebyshev(numTaps, arg);
    if (type == "flattop") return flattop(numTaps);
    return rectangular(numTaps);
}

} //namespace PothosComms::FilterDesign
