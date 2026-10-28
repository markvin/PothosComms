########################################################################
# Common DSP functions for the modules
########################################################################
pothos_comms_static_library(CommsFunctions fxpt_atan2.cpp)
target_include_directories(CommsFunctions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
