TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    functions.cpp

HEADERS += \
    functions.h

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -lblas -llapack -larmadillo




