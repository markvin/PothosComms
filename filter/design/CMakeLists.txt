########################################################################
# Filter design code of the filter and window blocks (replaces spuce)
########################################################################
pothos_comms_static_library(comms_filter_design
    Window.cpp
    FIRDesign.cpp
    IIRDesign.cpp
    SincBandStop.cpp
    Remez.cpp)
target_include_directories(comms_filter_design PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#temporary module for the tests, until the filter module comes back
pothos_comms_module(
    TARGET FilterDesignTests
    SOURCES
        TestFilterDesign.cpp
    LIBRARIES
        comms_filter_design
    TESTS
        test_filter_design_window
        test_filter_design_fir
        test_filter_design_remez_estimate
        test_filter_design_iir
        test_filter_design_iir_run
        test_filter_design_iir_taps
)
