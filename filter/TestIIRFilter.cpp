// Copyright (c) 2015-2015 Tony Kirke
// Copyright (c) 2016-2016 Josh Blum
// SPDX-License-Identifier: BSL-1.0

#include <Pothos/Testing.hpp>
#include <Pothos/Framework.hpp>
#include <Pothos/Proxy.hpp>
#include <cmath> //fabs
#include <iostream>
#include <string>
#include <utility>
#include <vector>

static double iirfilterToneGetRMS(
																	const double sampRate,
																	const double waveFreq
)
{
    auto waveSource = Pothos::BlockRegistry::make("/comms/waveform_source", "complex128");
    waveSource.call("setWaveform", "SINE");
    waveSource.call("setFrequency", waveFreq);
    waveSource.call("setSampleRate", sampRate);

    auto finiteRelease = Pothos::BlockRegistry::make("/blocks/finite_release");
    finiteRelease.call("setTotalElements", 4096);

    auto filter = Pothos::BlockRegistry::make("/comms/iir_filter", "complex128");
		filter.call("setWaitTaps", true);

    auto designer = Pothos::BlockRegistry::make("/comms/iir_designer");
    designer.call("setSampleRate", sampRate);
    designer.call("setIIRType", "butterworth");
    designer.call("setFilterType", "LOW_PASS");
    designer.call("setFrequencyLower", 0.1*sampRate);
    designer.call("setOrder", 4);

    auto probe = Pothos::BlockRegistry::make("/comms/signal_probe", "complex128");
    probe.call("setMode", "RMS");

    //run the topology
    {
        Pothos::Topology topology;
        topology.connect(designer, "tapsChanged", filter, "setTaps");
        topology.connect(waveSource, 0, finiteRelease, 0);
        topology.connect(finiteRelease, 0, filter, 0);
        topology.connect(filter, 0, probe, 0);
        topology.commit();
        POTHOS_TEST_TRUE(topology.waitInactive(0.01));
    }

    return probe.call<double>("value");
}

POTHOS_TEST_BLOCK("/comms/tests", test_iir_filter)
{

	const double rate = 1e6;
	const double freq = 30e3;
	std::cout << "freq " << freq << std::flush;
	auto rms = iirfilterToneGetRMS(rate, freq);
	std::cout << " RMS = " << rms << std::endl;
	POTHOS_TEST_TRUE(rms > 0.1);
}

/***********************************************************************
 * The text of the exception of a block call, empty without an exception
 **********************************************************************/
template <typename... Args>
static std::string callError(const Pothos::Proxy &block, const std::string &name, Args&&... args)
{
    try
    {
        block.call(name, std::forward<Args>(args)...);
    }
    catch (const Pothos::Exception &ex)
    {
        return ex.displayText();
    }
    return "";
}

/***********************************************************************
 * IIR filter taps: the first half is b, the second half is a
 **********************************************************************/
POTHOS_TEST_BLOCK("/comms/tests", test_iir_filter_taps)
{
    auto filter = Pothos::BlockRegistry::make("/comms/iir_filter", "float64");
    POTHOS_TEST_TRUE(callError(filter, "setTaps", std::vector<double>{}).find("Order cannot 0") != std::string::npos);
    POTHOS_TEST_TRUE(callError(filter, "setTaps", std::vector<double>{1.0}).find("Need at least 2 taps") != std::string::npos);
    POTHOS_TEST_EQUAL(callError(filter, "setTaps", std::vector<double>{0.5, 1.0}), "");
    POTHOS_TEST_EQUAL(callError(filter, "setTaps", std::vector<double>{1.0, 0.0, 1.0, -0.5, 9.0}), "");
}

/***********************************************************************
 * IIR designer: an unknown IIR type or band type gives an exception
 **********************************************************************/
POTHOS_TEST_BLOCK("/comms/tests", test_iir_designer_errors)
{
    auto designer = Pothos::BlockRegistry::make("/comms/iir_designer");
    designer.call("setSampleRate", 1.0);
    designer.call("setFrequencyLower", 0.1);
    designer.call("setOrder", 4);

    auto filter = Pothos::BlockRegistry::make("/comms/iir_filter", "float64");

    Pothos::Topology topology;
    topology.connect(designer, "tapsChanged", filter, "setTaps");
    topology.commit();

    POTHOS_TEST_TRUE(callError(designer, "setIIRType", "foo").find("unknown filter or band type") != std::string::npos);
    designer.call("setIIRType", "butterworth");
    POTHOS_TEST_TRUE(callError(designer, "setFilterType", "FOO").find("unknown filter or band type") != std::string::npos);
    designer.call("setFilterType", "LOW_PASS");

    topology.disconnectAll();
    topology.commit();
}
