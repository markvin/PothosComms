########################################################################
# Demodulation blocks module
########################################################################
pothos_comms_module(
    TARGET DemodBlocks
    SOURCES
        FreqDemod.cpp
    LIBRARIES
        CommsFunctions
)
