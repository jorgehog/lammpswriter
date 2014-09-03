TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = mpicxx
QMAKE_LINK = $$QMAKE_CXX

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK  #--enable-shared --enable-static

QMAKE_CXXFLAGS += -std=c++11

SOURCES += example.cpp
HEADERS += lammpswriter.h datahandler.h

