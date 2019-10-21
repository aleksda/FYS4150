TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        data_to_file.cpp \
        g_laguerre.cpp \
        g_legendre.cpp \
        main.cpp \
        mesh_and_weights.cpp

HEADERS += \
        data_to_file.h \
        g_laguerre.h \
        g_legendre.h \
        mesh_and_weights.h

