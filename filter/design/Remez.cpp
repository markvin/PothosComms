// Copyright (c) Tony Kirke (spuce), 2026 Pothos contributors
// SPDX-License-Identifier: BSL-1.0
// Ported from spuce ad08df7 (https://github.com/hayguen/spuce).

#include "FilterDesign.hpp"
#include "Detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>

namespace PothosComms::FilterDesign {

namespace {

constexpr double pi = std::numbers::pi;
constexpr int maxIterations = 40;
constexpr int gridDensity = 16;

//inverse DFT of a symmetric frequency response (spuce inv_dft_symmetric)
std::vector<double> invDFTSymmetric(const std::vector<double> &A, const int N)
{
    std::vector<double> h(N);
    const double M = (N - 1.0) / 2.0;
    const int end = (N % 2 != 0) ? (int)M : (N / 2 - 1);
    for (int n = 0; n < N; n++)
    {
        double val = A[0];
        const double x = 2.0*pi * (n - M) / N;
        for (int k = 1; k <= end; k++) val += 2.0 * A[k] * std::cos(x * k);
        h[n] = val / N;
    }
    return h;
}

//the number of dense grid points of a band (spuce createDenseGrid)
int bandGridPoints(const double low, const double high, const double delf)
{
    return (int)std::floor((high - low) / delf + 0.5);
}

/*!
 * The dense frequency grid, with the desired response d and the weight w
 * on the grid. spuce does not check that the points fit into gridSize:
 * remez() checks it before.
 */
void createDenseGrid(const int r, const std::vector<double> &bands,
    const std::vector<double> &des, const std::vector<double> &weight,
    std::vector<double> &grid, std::vector<double> &d, std::vector<double> &w)
{
    const double delf = 0.5 / (gridDensity * r);
    int j = 0;
    for (std::size_t band = 0; band < bands.size()/2; band++)
    {
        double lowf = bands[2 * band];
        const double highf = bands[2 * band + 1];
        const int k = bandGridPoints(lowf, highf, delf);
        for (int i = 0; i < k; i++)
        {
            d[j] = des[band];
            w[j] = weight[band];
            grid[j] = lowf;
            lowf += delf;
            j++;
        }
        grid[j - 1] = highf;
    }
}

//the "Dee" function in Parks & Burrus
std::vector<double> calcD(const int r, const std::vector<double> &x)
{
    std::vector<double> d(r+1);
    const int ld = (r - 1) / 15 + 1;
    for (int i = 0; i <= r; i++)
    {
        double denom = 1.0;
        const double xi = x[i];
        for (int j = 0; j < ld; j++)
        {
            for (int k = j; k <= r; k += ld)
            {
                if (k != i) denom *= 2.0 * (xi - x[k]);
            }
        }
        if (std::fabs(denom) < 0.00001) denom = 0.00001;
        d[i] = 1.0 / denom;
    }
    return d;
}

std::vector<double> calcY(const int r, const std::vector<int> &ext, const std::vector<double> &d,
    const std::vector<double> &w, const std::vector<double> &ad)
{
    double numer = 0, denom = 0, sign = 1;
    for (int i = 0; i <= r; i++)
    {
        numer += ad[i] * d[ext[i]];
        denom += sign * ad[i] / w[ext[i]];
        sign = -sign;
    }
    const double dev = numer / denom;
    sign = 1;
    std::vector<double> y(r+1);
    for (int i = 0; i <= r; i++)
    {
        y[i] = d[ext[i]] - sign * dev / w[ext[i]];
        sign = -sign;
    }
    return y;
}

//the "gee" function, p301 of Parks & Burrus, Digital Filter Design
double gee(const double freq, const int r, const std::vector<double> &ad,
    const std::vector<double> &x, const std::vector<double> &y)
{
    double numer = 0;
    double denom = 0;
    const double xc = std::cos(2*pi * freq);
    for (int i = 0; i <= r; i++)
    {
        double c = xc - x[i];
        if (std::fabs(c) < 1.0e-7)
        {
            numer = y[i];
            denom = 1;
            break;
        }
        c = ad[i] / c;
        denom += c;
        numer += c * y[i];
    }
    return numer / denom;
}

/*!
 * Find the maxima and minima of the error curve. If there are more than
 * r+1 extrema, delete them as follows:
 * 1. Delete adjacent extrema that do not alternate first.
 * 2. If there is more than one extra extremum, delete the one with the
 *    smallest error. Then 1. removes the non-alternation that this makes.
 * 3. If there is one extra extremum, delete the smaller of the first and
 *    the last extremum.
 */
void search(const int r, std::vector<int> &ext, const int gridSize, const std::vector<double> &e)
{
    //spuce reads one item after the found extrema when it deletes the last one
    std::vector<int> foundExt(std::max(gridSize, r + 1) + 1);
    int k = 0;
    for (int i = 0; i <= r; i++) foundExt[i] = ext[i];

    //extremum at 0
    if (((e[0] > 0.0) and (e[0] > e[1])) or ((e[0] < 0.0) and (e[0] < e[1]))) foundExt[k++] = 0;

    //extrema inside the dense grid
    for (int i = 1; i < gridSize - 1; i++)
    {
        if (((e[i] >= e[i - 1]) and (e[i] > e[i + 1]) and (e[i] > 0.0)) or
            ((e[i] <= e[i - 1]) and (e[i] < e[i + 1]) and (e[i] < 0.0)))
            foundExt[k++] = i;
    }

    //extremum at 0.5
    int j = gridSize - 1;
    if (((e[j] > 0.0) and (e[j] > e[j - 1])) or ((e[j] < 0.0) and (e[j] < e[j - 1]))) foundExt[k++] = j;

    //remove the extra extrema
    int extra = k - (r + 1);
    while (extra > 0)
    {
        //up: the first one is a maximum
        bool up = e[foundExt[0]] > 0.0;
        int l = 0;
        bool alt = true;
        for (j = 1; j < k; j++)
        {
            if (std::fabs(e[foundExt[j]]) < std::fabs(e[foundExt[l]])) l = j; //new smallest error
            if (up and (e[foundExt[j]] < 0.0)) up = false; //a minimum
            else if (not up and (e[foundExt[j]] > 0.0)) up = true; //a maximum
            else
            {
                //two extrema that do not alternate: delete the smallest
                alt = false;
                break;
            }
        }

        //one extra extremum and all alternate: delete the smaller of the
        //first and the last. As spuce, l is a grid index here, not a
        //position in foundExt.
        if (alt and (extra == 1))
        {
            if (std::fabs(e[foundExt[k - 1]]) < std::fabs(e[foundExt[0]])) l = foundExt[k - 1];
            else l = foundExt[0];
        }

        for (j = l; j < k; j++) foundExt[j] = foundExt[j + 1];
        k--;
        extra--;
    }
    for (int i = 0; i <= r; i++) ext[i] = foundExt[i];
}

//true when the error function is small enough (the result converged)
bool isDone(const int r, const std::vector<int> &ext, const std::vector<double> &e)
{
    double min = std::fabs(e[ext[0]]);
    double max = min;
    for (int i = 1; i <= r; i++)
    {
        const double current = std::fabs(e[ext[i]]);
        if (current < min) min = current;
        if (current > max) max = current;
    }
    return ((max - min) / max) < 0.0001;
}

/*!
 * spuce converts the result of ceil() to size_t directly. A negative value
 * is undefined behavior; GCC on x86-64 converts through int64_t, which
 * gives a very large number. The reference data has such values: use the
 * same conversion here, it is defined behavior.
 */
std::size_t toSize(const double x)
{
    return static_cast<std::size_t>(static_cast<std::int64_t>(x));
}

} //namespace

std::vector<double> Detail::remez(const int numTaps, const std::vector<double> &bands,
    const std::vector<double> &des, const std::vector<double> &weight)
{
    //positive symmetry (spuce remez_type::BANDPASS)
    int r = numTaps / 2; //the number of extrema
    if (numTaps % 2 != 0) r++;

    //the dense grid size in advance, for the array sizes
    const int numBands = bands.size()/2;
    int gridSize = 0;
    for (int i = 0; i < numBands; i++)
    {
        gridSize += (int)std::floor(0.5 + 2 * r * gridDensity * (bands[2 * i + 1] - bands[2 * i]));
    }
    if (gridSize < 1) throw std::runtime_error("Gridsize issue in Remez");

    //createDenseGrid uses another formula for the number of points: it
    //can be more than gridSize, and a band without points writes before
    //the grid. spuce writes outside of its arrays then.
    int gridPoints = 0;
    for (int i = 0; i < numBands; i++)
    {
        const int k = bandGridPoints(bands[2 * i], bands[2 * i + 1], 0.5 / (gridDensity * r));
        if (k <= 0) throw std::runtime_error("Gridsize issue in Remez");
        gridPoints += k;
    }
    if (gridPoints > gridSize) throw std::runtime_error("Gridsize issue in Remez");

    std::vector<double> grid(gridSize), d(gridSize), w(gridSize), e(gridSize);
    std::vector<double> x(r + 1), y(r + 1), ad(r + 1), taps(r + 1);
    std::vector<int> ext(r + 1);

    createDenseGrid(r, bands, des, weight, grid, d, w);

    //initial guess
    for (int i = 0; i <= r; i++) ext[i] = i * (gridSize - 1) / r;

    //even number of taps: change d[] and w[] as in Parks McClellan
    if (numTaps % 2 == 0)
    {
        for (int i = 0; i < gridSize; i++)
        {
            const double c = std::cos(pi * grid[i]);
            d[i] /= c;
            w[i] *= c;
        }
    }

    //the Remez exchange algorithm
    int iter = 0;
    for (; iter < maxIterations; iter++)
    {
        for (int i = 0; i <= r; i++) x[i] = std::cos(2*pi * grid[ext[i]]);
        ad = calcD(r, x);
        y = calcY(r, ext, d, w, ad);
        for (int i = 0; i < gridSize; i++) e[i] = w[i] * (d[i] - gee(grid[i], r, ad, x, y));
        search(r, ext, gridSize, e);
        if (isDone(r, ext, e)) break;
    }
    if (iter == maxIterations) throw std::runtime_error("Reached maximum iteration count. Results may be bad");

    for (int i = 0; i <= r; i++) x[i] = std::cos(2*pi * grid[ext[i]]);
    ad = calcD(r, x);
    y = calcY(r, ext, d, w, ad);

    //the taps of the filter for the frequency sampling design; with an
    //even number of taps, change them as in Parks McClellan
    for (int i = 0; i <= numTaps / 2; i++)
    {
        const double c = (numTaps % 2 != 0) ? 1 : std::cos(pi * (double)i / numTaps);
        taps[i] = gee((double)i / numTaps, r, ad, x, y) * c;
    }
    return invDFTSymmetric(taps, numTaps);
}

double remezEstimateWeight(const double passRippleDb, const double stopAttenDb)
{
    //the pass band ripple is relative to 1, the stop band ripple is the attenuation
    const double passRipple = std::pow(10.0, passRippleDb/20.0) - 1.0;
    const double stopRipple = std::pow(10.0, -stopAttenDb/20.0);
    const double maxDev = std::max(stopRipple, passRipple);
    if (passRipple > stopRipple) return maxDev/stopRipple;
    return 1.0/(maxDev/passRipple);
}

double remezEstimateBw(const std::size_t numTaps, const double ripple, const double atten)
{
    //decrease the transition bandwidth until the estimate needs numTaps taps
    double transBw = 0.45; //must be < 0.5
    std::size_t num = 1;
    do
    {
        num = remezEstimateNumTaps(transBw, ripple, atten);
        transBw *= 0.9;
    } while (num < numTaps);
    transBw /= 0.9;
    return transBw;
}

double remezEstimateAtten(const std::size_t numTaps, const double transBw, const double ripple)
{
    //increase the stop band attenuation until the estimate needs numTaps taps
    double atten = 1.0;
    std::size_t num = 1;
    do
    {
        num = remezEstimateNumTaps(transBw, ripple, atten);
        atten += 1.0;
    } while (num < numTaps);
    atten += 1.0;
    return atten;
}

std::size_t remezEstimateNumTaps(const double transBw, const double ripple, const double atten)
{
    const double AA[9] = {
        -4.278e-01, -4.761e-01, 0,
        -5.941e-01, 7.114e-02, 0,
        -2.660e-03, 5.309e-03, 0};
    //log10 of the ripples
    const double rl = std::pow(10.0, ripple/20.0) - 1;
    const double d1 = std::log(rl)/std::log(10.0);
    const double d2 = -atten/20.0;
    const double dd[3] = {1, d2, d2*d2};
    double AD[3];
    for (std::size_t j = 0; j < 3; j++)
    {
        double sum = 0;
        for (std::size_t i = 0; i < 3; i++) sum += AA[i + 3*j]*dd[i];
        AD[j] = sum;
    }
    const double D = AD[0] + d1*AD[1] + d1*d1*AD[2];
    const double fk = 11.01217 + 0.51244*(d1 - d2);
    return toSize(std::ceil(D/transBw - fk*transBw + 1.0));
}

} //namespace PothosComms::FilterDesign
