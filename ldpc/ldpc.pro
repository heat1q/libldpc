TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    ldpc_cycles.c \
    ldpc_main.c \
    ldpc_decoder.c \
    scm_functions.c \
    functions.c

HEADERS += \
    scm_types.h \
    ldpc_cycles.h \
    ldpc_types.h \
    ldpc_decoder.h \
    scm_functions.h \
    functions.h

