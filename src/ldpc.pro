TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += analysis_ldpc.c \
    sim_ldpc.c \
    LDPC_SSET/ldpc_cycles.c \
    decoder/ldpc_decoder.c \
    function/scm_functions.c \
    function/ldpc_functions.c \
    LDPC_SSET/ldpc_stoppingsets.c \
    decoder/ldpc_decoder_layered.c \
    decoder/main.c

HEADERS += \
    scm_types.h \
    LDPC_SSET/ldpc_cycles.h \
    function/ldpc_types.h \
    decoder/ldpc_decoder.h \
    function/scm_functions.h \
    function/ldpc_functions.h \
    LDPC_SSET/ldpc_stoppingsets.h \
    function/scm_types.h \
    decoder/ldpc_decoder_layered.h

