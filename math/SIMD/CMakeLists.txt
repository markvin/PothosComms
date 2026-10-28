########################################################################
## Make a static library with the SIMD implementations because MSVC
## doesn't like MathBlocksDocs.cpp depending on too many things.
########################################################################

set(SIMDInputs
    Abs.cpp
    Arithmetic.cpp
    Beta.cpp
    Comparator.cpp
    Conjugate.cpp
    ConstArithmetic.cpp
    ConstComparator.cpp
    ErrorFunction.cpp
    Gamma.cpp
    Exp.cpp
    Log.cpp
    ModF.cpp
    Pow.cpp
    Root.cpp
    RSqrt.cpp
    Sigmoid.cpp
    Sinc.cpp
    Trigonometric.cpp)

PothosGenerateSIMDSources(
    SIMDSources
    MathBlocks.json
    ${SIMDInputs})

pothos_comms_static_library(CommsMathSIMD ${SIMDSources})
target_link_libraries(CommsMathSIMD PRIVATE comms_xsimd)
#the SIMD dispatcher here, RSqrt.hpp and the generated Exp10.hpp in math/
target_include_directories(CommsMathSIMD PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/math
    ${PROJECT_BINARY_DIR}/math)
add_dependencies(CommsMathSIMD MathBlocks_SIMDDispatcher)

# This library is pure templates, so expect large object files.
if(MSVC)
    target_compile_options(CommsMathSIMD PRIVATE /bigobj)
endif()
