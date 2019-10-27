TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += D:\Armadillo\armadillo-9.200-RC1\include
DEPENDPATH += D:\Armadillo\armadillo-9.200-RC1\include

LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
    isingmodel.cpp

HEADERS += \
    isingmodel.hpp
