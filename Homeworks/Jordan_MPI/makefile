CXXFLAGS = -isystem /usr/local/intel/impi/5.1.3.223/intel64/include \
           -O3 -mfpmath=sse -fstack-protector-all -g \
           -W -Wall -Wextra -Wunused -Wcast-align \
           -Werror -pedantic -pedantic-errors -Wfloat-equal \
           -Wpointer-arith -Wformat-security -Wmissing-format-attribute \
           -Wformat=1 -Wwrite-strings -Wno-long-long \
           -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual \
           -Wno-suggest-attribute=format

TARGET = a.out

SRCS = main.cpp functions.cpp
HDRS = functions.h
OBJS = $($SRCS:.cpp=.o)


all:
	mpicxx $(CXXFLAGS) $(SRCS) -o $(TARGET)
    # g++ -g $(SRCS) -o $(TARGET)
