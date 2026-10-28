########################################################################
# Window blocks module
########################################################################
pothos_comms_module(
    TARGET WindowBlocks
    SOURCES
        WindowDesigner.cpp
    LIBRARIES
        comms_filter_design
        CommsFunctions
)
