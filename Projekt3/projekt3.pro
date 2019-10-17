TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        data_to_file.cpp \
        g_legendre.cpp \
        main.cpp
        GaussLegendre.cpp

HEADERS += \
        GaussLegendre.h \
        data_to_file.h \
        g_legendre.h

