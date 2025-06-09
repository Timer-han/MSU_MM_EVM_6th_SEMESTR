QT += widgets

QMAKE_CXXFLAGS += -std=c++14 -fsanitize=address

HEADERS       = window.h \
                algorithm.h \
                reduce.h
SOURCES       = main.cpp \
                window.cpp \
                algorithm.cpp \
                reduce.cpp
