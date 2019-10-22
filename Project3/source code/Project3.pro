TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


#INCLUDEPATH += D:\Armadillo\armadillo-9.200-RC1\include
#DEPENDPATH += D:\Armadillo\armadillo-9.200-RC1\include

#LIBS += -larmadillo -lblas -llapack
QMAKE_CXXFLAGS += -fopenmp -O2
QMAKE_LFLAGS += -fopenmp

SOURCES += \
    main.cpp \
    bruteForceMC.cpp \
    samplingMC.cpp \
    samplingMC_omp.cpp \
    lib.cpp \
    laguerre.cpp \
    legendre.cpp

HEADERS += \
    bruteForceMC.h \
    samplingMC.h \
    samplingMC_omp.h \
    lib.h \
    bruteForceMC.h \
    laguerre.h \
    legendre.h

#DISTFILES += \
