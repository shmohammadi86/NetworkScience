#CXXFLAGS=-Wall -O0 -g -Wno-deprecated -lgomp -fopenmp
CXXFLAGS=-Wall -O3 -g -Wno-deprecated
LIBS=-lboost_program_options-mt
#CXXFLAGS=-Wall -O0 -g -ffloat-store
#CXXFLAGS=-Wall -O2 -g -ffloat-store
BIN=msgsteiner


all: $(BIN)
msgsteiner: msgsteiner.cpp mes.o proba.hpp
	g++ ${CXXFLAGS} msgsteiner.cpp ${LIBS} mes.o -o msgsteiner
clean:
	rm -f $(BIN) *.o
mes.o: mes.cpp mes.hpp proba.hpp
	g++ ${CXXFLAGS} -c mes.cpp -o mes.o

