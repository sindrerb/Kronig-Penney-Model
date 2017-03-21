TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    vec3.cpp \
    kronigpenney.cpp \
    wavestate.cpp

HEADERS += \
    vec3.h \
    kronigpenney.h \
    wavestate.h
