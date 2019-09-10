TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += D:\Armadillo\armadillo-9.200-RC1\include
DEPENDPATH += D:\Armadillo\armadillo-9.200-RC1\include

LIBS += -larmadillo -lblas -llapack




SOURCES += \
    main.cpp \
    test_proj1.cpp
    project1b.cpp
    project1c.cpp

HEADERS += \
    project1b.h \
    project1c.h \
    project1e.h \
    project1d.h \
    test_proj1.h
