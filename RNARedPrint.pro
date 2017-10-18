
CONFIG += c++11

TARGET = RNARedPrint
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

SOURCES += src/DP.cpp \
    src/Nucleotide.cpp \
    src/RNARedPrint.cpp \
    src/RNAStructure.cpp \
    src/TreeDecomposition.cpp \
    src/Utils.cpp \
    src/EnergyModels.cpp

HEADERS += \
    src/DP.hpp \
    src/Nucleotide.hpp \
    src/RNAStructure.hpp \
    src/TreeDecomposition.hpp \
    src/Utils.hpp \
    src/EnergyModels.hpp

DISTFILES +=
