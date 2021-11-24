CC = gcc
CXX = g++
CXXFLAGS = -g -fPIC -Wall -O2 -lgslcblas -lgsl -lm -std=c++14
LIB = /usr/local/lib/
LDFLAGS = -lgsl -lgslcblas -lm

TARGET_LIB = libchi2.so

TARGET = solar

SRCS = numerical.cpp probability.cpp event.cpp read_files.cpp wrapper.cpp

OBJS = numerical.o probability.o event.o read_files.o wrapper.o

.PHONY: all

all: $(TARGET)

$(TARGET): $(TARGET).o $(OBJS)
	$(CXX) $(CXXFLAGS) -o  $(TARGET) $(TARGET).o $(OBJS) $(LDFLAGS)

$(TARGET).o: $(TARGET).cpp $(SRCS)
	$(CXX)  $(CXXFLAGS) -c $(TARGET).cpp -L$(LIB)




numerical.o: numerical.cpp numerical.hpp
	$(CXX) $(CXXFLAGS)  -c numerical.cpp -L$(LIB)


probability.o: probability.cpp probability.hpp
	$(CXX) $(CXXFLAGS) -c probability.cpp


read_files.o: read_files.cpp read_files.hpp
	$(CXX) $(CXXFLAGS) -c read_files.cpp


wrapper.o: wrapper.cpp wrapper.hpp
	$(CXX) $(CXXFLAGS) -c wrapper.cpp


event.o: event.cpp event.hpp
	$(CXX) $(CXXFLAGS) -c event.cpp -L$(LIB)








clear_obj:
	rm -f *.o


clean:
	rm -f *.o *.out $(TARGET) $(TARGET1)
