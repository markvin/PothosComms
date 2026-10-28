########################################################################
# MAC layer and packet protocol blocks module
########################################################################
pothos_comms_module(
    TARGET MacBlocks
    SOURCES
        SimpleMac.cpp
        TestSimpleMac.cpp
        SimpleLlc.cpp
        TestSimpleLlc.cpp
    TESTS
        test_simple_llc
        test_simple_llc_harsh
        test_simple_mac
)
