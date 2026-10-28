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

