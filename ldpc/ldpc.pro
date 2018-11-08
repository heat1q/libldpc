TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    functions.c \
    ldpc_cycles.c \
    ldpc_main.c

HEADERS += \
    scm_types.h \
    functions.h \
    ldpc_cycles.h \
    ldpc_types.h

