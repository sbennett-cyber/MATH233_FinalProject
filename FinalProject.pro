TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ../../Important_Files/conjugate_gradient.cpp \
        ../../Important_Files/fvsolver.cpp \
        ../../Important_Files/grid2d.cpp \
        ../../Important_Files/sparsematrix_crs.cpp \
        main.cpp \
        networkfinitedifference.cpp

## Stuff for running in parallel
macx:
{
QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -I/usr/local/include
}
macx:
{
QMAKE_LFLAGS += -lomp
}
macx:
{
LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib
}

HEADERS += \
    ../../Important_Files/cf_2.h \
    ../../Important_Files/conjugate_gradient.h \
    ../../Important_Files/fvsolver.h \
    ../../Important_Files/grid2d.h \
    ../../Important_Files/matrix.h \
    ../../Important_Files/sparsematrix_crs.h \
    networkfinitedifference.h
