TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += D:\Armadillo\armadillo-9.200-RC1\include
DEPENDPATH += D:\Armadillo\armadillo-9.200-RC1\include

LIBS += -larmadillo -lblas -llapack

SOURCES += \
    jacobi.cpp \ 
    tridiagtoplitz.cpp \
    #tridiagtest.cpp \
    #test.cpp \
    #maintest.cpp \
    main.cpp

HEADERS += \
    jacobi.h \
    tridiagtoplitz.h \
    #tridiagtest.h \
    #catch.hpp

#DISTFILES += \
