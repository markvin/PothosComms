########################################################################
# Utility blocks module
########################################################################
pothos_comms_module(
    TARGET UtilityBlocks
    SOURCES
        SignalProbe.cpp
        Threshold.cpp
        WaveTrigger.cpp
        SplitComplex.cpp
        CombineComplex.cpp
        TestComplex.cpp
    TESTS
        test_complex_split_combine
)
