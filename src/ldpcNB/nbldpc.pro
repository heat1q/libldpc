TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ldpc.cpp

HEADERS += \
    ldpc.h \
    gf_math.h

