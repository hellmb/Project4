TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    functions_mpi.cpp

HEADERS += \
    functions_mpi.h

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -lblas -llapack -larmadillo

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib/openmpi
LIBS += -L/usr/local/opt/libevent/lib
LIBS += -L/usr/local/Cellar/open-mpi/2.0.1/lib -lmpi
