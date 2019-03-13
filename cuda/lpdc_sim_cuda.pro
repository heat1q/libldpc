TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    sim_cuda.cpp \
    ldpc/ldpc.cpp \
    simulation.cpp \
    ldpc/decoder_layered.cpp

HEADERS += \
    ldpc/ldpc.h \
    simulation.h
