// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

#include "FilterDesign.hpp"

#include <cfloat>
#include <cmath>
#include <complex>
#include <numbers>

namespace PothosComms::FilterDesign {

namespace {

constexpr double pi = std::numbers::pi;
using Complex = std::complex<double>;

enum class BandType {low, high, bandpass, bandstop};

//c = x convolved with y
std::vector<double> convolve(const std::vector<double> &x, const std::vector<double> &y)
{
    std::vector<double> c(x.size() + y.size() - 1);
    for (std::size_t i = 0; i < c.size(); i++)
    {
        for (std::size_t j = 0; j < x.size(); j++)
        {
            if (j <= i and i - j < y.size()) c[i] += x[j] * y[i - j];
        }
    }
    return c;
}

//the first N items of x convolved with the first M items of y
std::vector<double> partialConvolve(const std::vector<double> &x, const std::vector<double> &y, const int N, const int M)
{
    std::vector<double> c(M + N - 1);
    for (int i = 0; i < M + N - 1; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if ((i - j >= 0) and (i - j < M)) c[i] += x[j] * y[i - j];
        }
    }
    return c;
}

/*!
 * The poles and zeros of an IIR filter, and their transformation to the
 * transfer function (the part of spuce iir_coeff that design_iir uses).
 * poles and zeros keep one item of each conjugate pair; for an odd order,
 * item 0 is the real pole or zero.
 */
struct IIRCoeff
{
    IIRCoeff(const int ord, const BandType t):
        poles((ord + 1) / 2), zeros((ord + 1) / 2), aTf(ord + 1), bTf(ord + 1), type(t)
    {
        order = ord;
        n2 = (order + 1) / 2;
        odd = (order % 2);
        centerFreq = 0.25;
        c0 = std::cos(centerFreq*2*pi);
    }

    void setCenter(const double c)
    {
        centerFreq = c;
        c0 = std::cos(2*pi*c);
    }

    //keep the type and the center frequency, clear all other values
    void resize(const int ord)
    {
        order = ord;
        gain = hpfGain = 1.0;
        n2 = (order + 1) / 2;
        odd = (order % 2);
        poles.assign(n2, Complex(0.0, 0.0));
        zeros.assign(n2, Complex(0.0, 0.0));
        aTf.resize(ord + 1);
        bTf.resize(ord + 1);
    }

    //from the s plane to the z plane
    void bilinear()
    {
        hpfGain = 1.0;
        if (odd) hpfGain = 1.0 + std::real(poles[0]);
        for (int j = 0; j < n2; j++)
        {
            zeros[j] = (1.0 - zeros[j]) / (1.0 + zeros[j]);
            poles[j] = (1.0 - poles[j]) / (1.0 + poles[j]);
        }
    }

    //low pass to band pass or band stop around the center frequency c
    void makeBand(const double c)
    {
        centerFreq = c;
        c0 = std::cos(centerFreq*2*pi);
        const auto oldPoles = poles;
        const auto oldZeros = zeros;
        const int wasOdd = odd;
        resize(2 * order);

        //q is 1 for band pass, -1 for band stop
        const double q = (type == BandType::bandpass) ? 1 : -1;
        int k = 0;
        //for an odd original order, skip the first pole and zero
        for (int j = wasOdd; j < wasOdd + n2 / 2; j++)
        {
            const Complex pi_ = -oldPoles[j];
            const Complex zi = -oldZeros[j];
            const Complex p0 = c0 * (1.0 + q * pi_);
            const Complex z0 = c0 * (1.0 + q * zi);
            poles[k] = 0.5 * (p0 + std::sqrt(p0 * p0 - (4 * q) * pi_));
            poles[k + n2 / 2] = (c0 - poles[k]) / (1.0 - c0 * poles[k]);
            zeros[k] = 0.5 * (z0 + std::sqrt(z0 * z0 - (4 * q) * zi));
            zeros[k + n2 / 2] = (c0 - zeros[k]) / (1.0 - c0 * zeros[k]);
            k++;
        }

        //poles and zeros to the transfer function (without the first pole
        //and zero of an odd original order)
        convertToAb();
        if (not wasOdd) return;

        //replace the last 2 items with the second order section of the
        //first pole and zero (the polynomials are real)
        aTf.pop_back();
        aTf.pop_back();
        bTf.pop_back();
        bTf.pop_back();
        std::vector<double> p{1, 0, 0};
        //band pass zeros are the same for all center frequencies
        std::vector<double> z{1, 0, -1};
        if (type != BandType::bandpass)
        {
            z[1] = -2 * c0;
            z[2] = 1;
        }
        const Complex p0h = -oldPoles[0];
        const Complex p0hS = c0 * (1.0 + q * p0h);
        const Complex p0p = 0.5 * (p0hS + std::sqrt(p0hS * p0hS - (4 * q) * p0h));
        p[1] = -2 * std::real(p0p);
        p[2] = std::norm(p0p);
        bTf = convolve(z, bTf);
        aTf = convolve(p, aTf);
    }

    //poles and zeros to the transfer function, with the gain in b
    void convertToAb()
    {
        gain = 1.0;
        zRootToAb(zeros);
        const double zGain = gain;
        const double hpfZGain = hpfGain;
        gain = 1.0;
        hpfGain = 1.0;
        zRootToAb(poles);
        const double hpfPGain = hpfGain;
        gain /= zGain;
        if (odd) gain *= 0.5 * (1.0 - std::real(poles[0]));
        hpfGain = hpfPGain / hpfZGain;
        if (type == BandType::high) gain = hpfGain;
        aTf = p2ToPoly(poles);
        bTf = p2ToPoly(zeros);
        applyGain(gain);
    }

    //each root pair to the second order section 1 + a*z + b*z^2, packed as (a, b)
    void zRootToAb(std::vector<Complex> &z)
    {
        for (int j = odd; j < n2; j++)
        {
            gain *= (std::norm(z[j]) - 2 * std::real(z[j]) + 1.0);
            hpfGain *= (std::norm(z[j]) + 2 * std::real(z[j]) + 1.0);
            z[j] = Complex(-2 * std::real(z[j]), std::norm(z[j]));
        }
        if (gain == 0.0) gain = 1.0;
    }

    //convolve the second order sections (a, b) (and the first order section of an odd order)
    std::vector<double> p2ToPoly(const std::vector<Complex> &ab) const
    {
        std::vector<double> tf;
        std::vector<double> p2(3);
        std::vector<double> p(order + 1);
        p[0] = 1;
        p2[0] = 1;
        int m = 1;
        if (odd)
        {
            p2[1] = -std::real(ab[0]);
            p2[2] = 0;
            tf = partialConvolve(p, p2, m, 2);
            m += 1;
            for (int i = 0; i < m; i++) p[i] = tf[i];
        }
        for (int j = odd; j < n2; j++)
        {
            p2[1] = std::real(ab[j]);
            p2[2] = std::imag(ab[j]);
            tf = partialConvolve(p, p2, m, 3);
            m += 2;
            for (int i = 0; i < m; i++) p[i] = tf[i];
        }
        return tf;
    }

    void applyGain(const double g)
    {
        for (auto &b : bTf) b *= g;
    }

    //the magnitude of the frequency response at freq (radians)
    double freqzMag(const double freq) const
    {
        Complex z(1, 0);
        const Complex zInc(std::cos(freq), std::sin(freq));
        Complex nom(0), den(0);
        for (int i = 0; i < order + 1; i++)
        {
            nom += z * bTf[i];
            den += z * aTf[i];
            z *= zInc;
        }
        return std::sqrt(std::norm(nom / den));
    }

    //gain 1 at the center frequency
    void setBandpassGain()
    {
        const double w = 2.0*pi*centerFreq;
        applyGain(1.0 / freqzMag(w));
    }

    //after the prototype in the s plane: to the z plane and to the band
    void finish()
    {
        bilinear();
        if (type == BandType::bandpass or type == BandType::bandstop) makeBand(centerFreq);
        else convertToAb();
        if (type == BandType::bandpass) setBandpassGain();
    }

    std::vector<Complex> poles, zeros;
    double gain = 1.0;
    double c0 = 0.0;
    double centerFreq = 0.25;
    double hpfGain = 1.0;
    int n2 = 0;
    int odd = 0;
    int order = 0;
    std::vector<double> aTf, bTf;
    BandType type;
};

//the analog zero at infinity of spuce (FLT_MAX)
const Complex infiniteZero(FLT_MAX, 0.0);

/***********************************************************************
 * Butterworth
 **********************************************************************/
void butterworthS(IIRCoeff &filt, const double wp, const long n)
{
    int l = (n % 2 == 0) ? 1 : 0;
    const int n2 = (n + 1) / 2;
    for (int j = 0; j < n2; j++)
    {
        const double arg = -0.5 * pi * l / ((double)(n));
        if (filt.type == BandType::low)
        {
            filt.poles[j] = wp * std::exp(Complex(0, arg));
            filt.zeros[j] = infiniteZero;
        }
        else
        {
            filt.poles[j] = (1.0 / wp) * std::exp(Complex(0, arg));
            filt.zeros[j] = 0;
        }
        l += 2;
    }
}

//fcd = the cutoff frequency (1 = the sample rate), amax = the attenuation at the cutoff
void butterworthIIR(IIRCoeff &filt, const double fcd, const double amax)
{
    const double ten = 10.0;
    const auto order = filt.order;
    const double epi = std::pow((double)(std::pow(ten, amax / ten) - 1.0), (double)(1. / (2.0 * order)));
    //the pre-warped angular frequency
    const double wca = (filt.type == BandType::high) ? std::tan(pi * (0.5-fcd)) / epi : std::tan(pi*fcd) / epi;
    butterworthS(filt, wca, order);
    filt.finish();
}

/***********************************************************************
 * Chebyshev
 **********************************************************************/
void chebyshevS(IIRCoeff &filt, const double wp, const double epi, const long n)
{
    int l = 1;
    const int n2 = (n + 1) / 2;
    const double x = 1 / epi;
    const double lambda = std::pow(x*(1.0 + std::sqrt(1.0 + epi*epi)), 1.0/n);
    const double sm = 0.5*((1.0/lambda) - lambda);
    const double cm = 0.5*((1.0/lambda) + lambda);
    for (int j = 0; j < n2; j++)
    {
        const double arg = 0.5 * pi * (2*l-1) / ((double)(n));
        if (filt.type == BandType::low)
        {
            filt.poles[n2-1-j] = -wp * Complex(sm * std::sin(arg), cm * std::cos(arg));
            filt.zeros[n2-1-j] = infiniteZero;
        }
        else
        {
            filt.poles[n2-1-j] = -1.0 / (wp * Complex(sm * std::sin(arg), cm * std::cos(arg)));
            filt.zeros[n2-1-j] = 0;
        }
        l++;
    }
}

//ripple = the pass band ripple in dB
void chebyshevIIR(IIRCoeff &filt, const double fcd, const double ripple)
{
    const double ten = 10.0;
    const double rlin = std::pow(ten, ripple/ten);
    const double epi = std::sqrt(rlin - 1.0);
    const double wca = (filt.type == BandType::high) ? std::tan(pi * (0.5-fcd)) : std::tan(pi*fcd);
    chebyshevS(filt, wca, epi, filt.order);
    filt.finish();
    //an even order filter needs this gain
    if (not filt.odd) filt.applyGain(1.0/std::sqrt(rlin));
}

/***********************************************************************
 * Chebyshev type 2 (inverse Chebyshev)
 **********************************************************************/
void chebyshev2S(IIRCoeff &filt, const double wp, const double epi, const long order)
{
    int l = 1;
    const int n2 = (order + 1) / 2;
    const double x = 1 / epi;
    const double lambda = std::pow(x*(1.0 + std::sqrt(1.0 + epi*epi)), 1.0/order);
    const double sm = 0.5*((1.0/lambda) - lambda);
    const double cm = 0.5*((1.0/lambda) + lambda);
    for (int j = 0; j < n2; j++)
    {
        const double arg = pi * (2*l-1) / ((double)(2*order));
        const Complex p(sm * std::sin(arg), cm * std::cos(arg));
        if (filt.type == BandType::low)
        {
            //the regular Chebyshev pole, transformed for the inverse Chebyshev
            filt.poles[n2-1-j] = (-wp*p/std::norm(p));
            filt.zeros[n2-1-j] = Complex(0, wp/std::cos(arg));
        }
        else
        {
            filt.poles[n2-1-j] = -p/wp;
            filt.zeros[n2-1-j] = Complex(0, std::cos(arg)/wp);
        }
        l++;
    }
}

//stopband = the stop band attenuation in dB
void chebyshev2IIR(IIRCoeff &filt, const double fcd, const double stopband)
{
    const double ten = 10.0;
    const double delta = std::pow(ten, -stopband/20.0);
    const double epi = delta/std::sqrt(1 - delta*delta);
    const double wca = (filt.type == BandType::high) ? std::tan(pi * (0.5-fcd)) : std::tan(pi*fcd);
    chebyshev2S(filt, wca, epi, filt.order);
    filt.finish();
}

/***********************************************************************
 * Elliptic
 **********************************************************************/
constexpr double machep = 2.0 * FLT_MIN;

//Jacobian elliptic functions sn, cn and dn of u, parameter m
void ellpj(const double u, const double m, double &sn, double &cn, double &dn)
{
    double a[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double c[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    a[0] = 1.0;
    double b = std::sqrt(1.0 - m);
    c[0] = std::sqrt(m);
    double twon = 1.0;
    int i = 0;
    while ((std::fabs(c[i] / a[i]) > machep) and (i < 8))
    {
        const double ai = a[i];
        ++i;
        c[i] = (ai - b) / 2.0;
        const double t = std::sqrt(ai * b);
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }
    double phi = twon * a[i] * u;
    //spuce runs this loop also for i = 0 and then reads c[-1]: stop at i = 0
    do
    {
        const double t = c[i] * std::sin(phi) / a[i];
        b = phi;
        phi = (std::asin(t) + phi) / 2.0;
    } while (--i > 0);
    sn = std::sin(phi);
    const double t = std::cos(phi);
    cn = t;
    dn = t / std::cos(phi - b);
}

//incomplete elliptic integral of the first kind, phi in radians, modulus k (0..1)
double ellik(const double phi, const double k)
{
    double a0 = 1.0;
    double b0 = std::sqrt(1.0 - k * k);
    double d0 = phi;
    if (k == 1.0) return std::log((1.0 + std::sin(d0)) / std::cos(d0));
    double fac = 1.0;
    double a = 0, d = 0;
    for (int n = 1; n < 41; n++)
    {
        a = 0.5 * (a0 + b0);
        const double b = std::sqrt(a0 * b0);
        const double c = 0.5 * (a0 - b0);
        fac *= 2.0;
        d = d0 + std::atan((b0 / a0) * std::tan(d0));
        d0 = d + pi * int(d / pi + 0.5);
        a0 = a;
        b0 = b;
        if (c < 1.0e-7) break;
    }
    return d / (fac * a);
}

//complete elliptic integral of the first kind, parameter k
double ellpk(const double k)
{
    double a0 = 1.0;
    double b0 = std::sqrt(1.0 - k);
    double c = std::sqrt(k);
    if (k == 1.0) return FLT_MAX;
    double a = 0;
    for (int n = 1; n < 41; n++)
    {
        a = 0.5 * (a0 + b0);
        c = 0.5 * (a0 - b0);
        const double b = std::sqrt(a0 * b0);
        a0 = a;
        b0 = b;
        if (c < 1.0e-7) break;
    }
    return pi / (2.0 * a);
}

//the modulus from the ratio K/K' (m = (theta2/theta3)^4)
double msqrt(const double u)
{
    double dt1;
    double a = 1.0;
    double b = 1.0;
    double r = 1.0;
    const double q = std::exp(-pi * u);
    double p = q;
    do
    {
        r *= p;
        a += 2.0 * r;
        dt1 = r / a;
        r *= p;
        b += r;
        p *= q;
    } while (dt1 > 1.0e-7);
    a = b / a;
    a = 4.0 * std::sqrt(q) * a * a;
    return a;
}

//the roots in the lambda plane
double lamdaPlane(const double m, const int n, const double eps)
{
    const double Kk = ellpk(m);
    const double Kpk = ellpk(1 - m);
    double u = (n * Kpk / Kk);
    double m1 = msqrt(u);
    m1 *= m1;
    const double Kk1 = ellpk(m1);
    const double phi = std::atan(1.0 / eps);
    u = Kk * ellik(phi, 1.0 - m1) / (n * Kk1);
    return u;
}

//the s plane poles and zeros
void sPlane(IIRCoeff &filt, const int n, const double u, const double m, const double k, const double Kk, const double wc)
{
    double sn1, cn1, dn1;
    double sn, cn, dn;
    const double kon = Kk / (double)n;
    ellpj(u, 1.0 - m, sn1, cn1, dn1);
    for (int i = 0; i < (n + 1) / 2; i++)
    {
        const double b = (n - 1 - 2 * i) * kon;
        const int j = (n - 1) / 2 - i;
        ellpj(b, m, sn, cn, dn);
        double r = k * sn * sn1;
        r = 1.0 / (cn1 * cn1 + r * r);
        if (filt.type == BandType::low)
        {
            if (sn != 0) filt.zeros[j] = -Complex(0.0, wc / (k * sn));
            else filt.zeros[j] = infiniteZero;
            filt.poles[j] = -Complex(-wc * cn * dn * sn1 * cn1 * r, wc * sn * dn1 * r);
        }
        else
        {
            if (sn != 0) filt.zeros[j] = -Complex(0.0, -1.0 / (wc / (k * sn)));
            else filt.zeros[j] = 0;
            filt.poles[j] = -1.0 / Complex(-wc * cn * dn * sn1 * cn1 * r, wc * sn * dn1 * r);
        }
    }
}

//ripple = the pass band ripple in dB, stopattn = the stop band attenuation
//in dB (below 1 dB: the stop band edge)
void ellipticIIR(IIRCoeff &filt, const double fcd, const double ripple, const double stopattn)
{
    const double ten = 10.0;
    const auto order = filt.order;
    double epi = std::pow(ten, (ripple / ten));
    epi = std::sqrt(epi - 1.0);
    //the pre-warped angular frequency
    const double wca = std::tan(pi * fcd);
    const double wc = (filt.type == BandType::high) ? std::tan(pi * (0.5-fcd)) : std::tan(pi*fcd);
    double fstp;
    if (stopattn > 1.0)
    {
        const double a = std::pow(ten, (stopattn / ten));
        const double m1 = epi * epi / (a - 1.0);
        const double m2 = 1 - m1;
        const double Kk1 = ellpk(m1);
        const double Kpk1 = ellpk(m2);
        const double u = Kpk1 / (order * Kk1);
        const double k = msqrt(u);
        const double wr = 1.0 / k;
        fstp = std::atan(wca * wr) / pi;
    }
    else
    {
        fstp = stopattn;
    }
    double wr = std::tan(fstp * pi) / wca;
    if (wr < 0.0) wr = -wr;
    //the selectivity factor
    const double k = 1.0 / wr;
    const double m = k * k;
    const double Kk = ellpk(m);
    const double u = lamdaPlane(m, order, epi);
    sPlane(filt, order, u, m, k, Kk, wc);
    filt.finish();
    //an even order filter needs this gain
    if (not filt.odd) filt.applyGain(std::pow(ten, -(ripple / 20.0)));
}

} //namespace

std::optional<IIRCoeffs> designIIR(const std::string_view type, const std::string_view band,
    const int order, const double fcd, const double ripple, const double stopAtten,
    const double centerFreq)
{
    BandType bandType;
    if (band == "LOW_PASS") bandType = BandType::low;
    else if (band == "HIGH_PASS") bandType = BandType::high;
    else if (band == "BAND_PASS") bandType = BandType::bandpass;
    else if (band == "BAND_STOP") bandType = BandType::bandstop;
    else return std::nullopt;

    IIRCoeff filt(order, bandType);
    //only for band pass and band stop
    filt.setCenter(centerFreq);

    if (type == "butterworth") butterworthIIR(filt, fcd, 3.0);
    else if (type == "chebyshev") chebyshevIIR(filt, fcd, ripple);
    else if (type == "chebyshev2") chebyshev2IIR(filt, fcd, stopAtten);
    else if (type == "elliptic") ellipticIIR(filt, fcd, ripple, stopAtten);
    else return std::nullopt;

    return IIRCoeffs{filt.bTf, filt.aTf};
}

} //namespace PothosComms::FilterDesign
