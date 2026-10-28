########################################################################
# Math blocks module
########################################################################
#exp10 and exp10f are GNU extensions, MSVC does not have them
if(NOT MSVC)
    include(CheckSymbolExists)
    set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
    set(CMAKE_REQUIRED_LIBRARIES m)
    check_symbol_exists(exp10f "math.h" HAVE_EXP10F)
    check_symbol_exists(exp10 "math.h" HAVE_EXP10)
endif()
configure_file(Exp10.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/Exp10.hpp)

pothos_comms_module(
    TARGET MathBlocks
    SOURCES
        Abs.cpp
        TestAbs.cpp
        Angle.cpp
        TestAngle.cpp
        Arithmetic.cpp
        TestComparatorBlocks.cpp
        Comparator.cpp
        ConstComparator.cpp
        Conjugate.cpp
        TestConjugate.cpp
        Log.cpp
        TestLog.cpp
        TestArithmeticBlocks.cpp
        Scale.cpp
        TestScale.cpp
        Rotate.cpp
        TestRotate.cpp
        ConstArithmetic.cpp
        Sinc.cpp
        TestSinc.cpp
        Trigonometric.cpp
        TestTrigonometric.cpp
        Pow.cpp
        Root.cpp
        TestPowRoot.cpp
        RSqrt.cpp
        TestRSqrt.cpp
        Beta.cpp
        Gamma.cpp
        TestBeta.cpp
        TestGamma.cpp
        ErrorFunction.cpp
        TestErrorFunction.cpp
        Sigmoid.cpp
        TestSigmoid.cpp
        Exp.cpp
        TestExp.cpp
        ModF.cpp
        TestModF.cpp
    LIBRARIES
        CommsFunctions
    TESTS
        test_beta
        test_conjugate
        test_trigonometric
        test_gamma
        test_modf
        test_abs
        test_scale
        test_sinc
        test_rotate
        test_arithmetic
        test_arithmetic_feedback
        test_inline_buffer
        test_const_arithmetic
        test_log
        test_log2
        test_log10
        test_log1p
        test_angle
        test_comparator
        test_const_comparator
        test_sigmoid
        test_rsqrt
        test_exp
        test_exp2
        test_exp10
        test_expm1
        test_pow_root
        test_erf
)
target_include_directories(MathBlocks PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})
if(MSVC)
    target_compile_options(MathBlocks PRIVATE /bigobj) #arithmetic is too large in debug mode
endif()

#the SIMD kernels, chosen at run time (the dispatcher has a scalar fallback)
add_subdirectory(SIMD)
target_link_libraries(MathBlocks PRIVATE CommsMathSIMD)
