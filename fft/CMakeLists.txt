########################################################################
# FFT blocks module
########################################################################
pothos_comms_module(
    TARGET FFTBlocks
    SOURCES
        FFT.cpp
        kiss_fft.c
        TestFFT.cpp
    TESTS
        test_fft_float
        test_fft_short
)
#The templated c++ implementation only handles floats.
#The C implementation handles the fixed point, but only one integer
#type at a time: complex int16 is the most useful one.
target_compile_definitions(FFTBlocks PRIVATE FIXED_POINT=16 KISS_FFT_USE_ALLOCA)
include(CheckIncludeFiles)
check_include_files(alloca.h HAS_ALLOCA_H)
if(HAS_ALLOCA_H)
    target_compile_definitions(FFTBlocks PRIVATE HAS_ALLOCA_H)
endif()
