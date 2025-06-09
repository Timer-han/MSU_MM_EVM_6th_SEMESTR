QT += widgets

QMAKE_CXXFLAGS += -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format 

HEADERS       = window.h \
                thread_work.h \
                solve.h
SOURCES       = main.cpp \
                window.cpp \
                thread_work.cpp \
                solve.cpp
